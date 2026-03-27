"""
Flat-Band Hotspot Hunter — Corrected CRT-Chirality Bridge
==========================================================
Fixes the CRITICAL BUG in crt_chirality_bridge.py:
  - OLD (WRONG): |A| = sqrtm(A^T A)  ← matrix absolute value
  - NEW (CORRECT): |A|_ew = np.abs(A)  ← element-wise absolute value

Two paths:
  PATH A: Verify Identity 2.5 with element-wise |A|, test CRT block-diagonality
  PATH B: Flat-band autopsy — rank eigenvectors by IPR, measure Block Participation Ratio

Tests m=2310 (even, 30×16) and m=1155 (odd, 60×8).
Runtime: ~30s on consumer PC (no sqrtm needed!)
"""

import numpy as np
from math import gcd, sqrt


def get_coprime_residues(m):
    return sorted([r for r in range(1, m) if gcd(r, m) == 1])


def build_D_sym(residues, m):
    r = np.array(residues)
    diff = np.abs(r[:, None] - r[None, :])
    return np.minimum(diff, m - diff).astype(float)


def build_directed_D(residues, m):
    r = np.array(residues)
    return ((r[None, :] - r[:, None]) % m).astype(float)


def build_P_tau(residues, m):
    n = len(residues)
    res_to_idx = {r: i for i, r in enumerate(residues)}
    perm = np.array([res_to_idx[pow(r, -1, m)] for r in residues])
    P = np.zeros((n, n))
    P[np.arange(n), perm] = 1.0
    return P


def get_crt_blocks(residues, block_primes):
    blocks = {}
    for i, r in enumerate(residues):
        sig = tuple(round(np.cos(2 * np.pi * r / p), 8) for p in block_primes)
        blocks.setdefault(sig, []).append(i)
    return blocks


def off_diagonal_block_norm(M, blocks):
    block_list = list(blocks.values())
    total_off = 0.0
    total_on = 0.0
    for i, bi in enumerate(block_list):
        for j, bj in enumerate(block_list):
            sub = M[np.ix_(bi, bj)]
            norm_sq = np.sum(np.abs(sub) ** 2)
            if i == j:
                total_on += norm_sq
            else:
                total_off += norm_sq
    total = total_on + total_off
    return total_off, total_on, total_off / total if total > 0 else 0.0


def build_V_diag(residues, p):
    return np.array([np.cos(2 * np.pi * r / p) for r in residues])


def ipr(v):
    """Inverse participation ratio for a normalized state."""
    p = np.abs(v) ** 2
    return np.sum(p ** 2)


def block_participation_ratio(v, blocks):
    """How many CRT blocks does state v live in?
    BPR = 1/sum(w_k^2) where w_k = sum_{i in block_k} |v_i|^2.
    BPR = 1 means single block; BPR = N_blocks means uniform."""
    block_list = list(blocks.values())
    weights = []
    for indices in block_list:
        w = np.sum(np.abs(v[indices]) ** 2)
        weights.append(w)
    weights = np.array(weights)
    total = weights.sum()
    if total < 1e-15:
        return 0.0
    weights /= total
    return 1.0 / np.sum(weights ** 2)


# ==============================================================
# MAIN EXPERIMENT
# ==============================================================

def run_experiment(m, label):
    print(f"\n{'='*70}")
    print(f"  {label}: m = {m}")
    print(f"{'='*70}")

    residues = get_coprime_residues(m)
    phi = len(residues)
    print(f"  φ = {phi}")

    # Build operators
    D_sym = build_D_sym(residues, m)
    D_dir = build_directed_D(residues, m)
    P_tau = build_P_tau(residues, m)

    # Antisymmetric part
    A = (D_dir - D_dir.T) / 2.0

    # ============================================================
    # PATH A: Identity 2.5 with ELEMENT-WISE |A|
    # ============================================================
    print(f"\n  ======= PATH A: Identity 2.5 verification =======")

    n = phi
    J = np.ones((n, n))
    I_mat = np.eye(n)

    # Element-wise absolute value (THE CORRECT ONE)
    A_abs_ew = np.abs(A)

    # Identity 2.5: D_sym = m/2 * (J - I) - |A|_ew
    rhs = (m / 2.0) * (J - I_mat) - A_abs_ew
    identity_err = np.linalg.norm(D_sym - rhs, 'fro')
    identity_rel = identity_err / np.linalg.norm(D_sym, 'fro')
    print(f"  D_sym = m/2*(J-I) - |A|_ew")
    print(f"    Absolute error:  {identity_err:.2e}")
    print(f"    Relative error:  {identity_rel:.2e}")
    if identity_rel < 1e-13:
        print(f"    ✓ MACHINE ZERO — Identity 2.5 CONFIRMED with element-wise |A|")
    else:
        print(f"    ✗ FAILED (relative error {identity_rel:.2e})")

    # Also verify Identity 2: [D_sym, P_tau] = -[|A|_ew, P_tau]
    comm_D = D_sym @ P_tau - P_tau @ D_sym
    comm_A = A_abs_ew @ P_tau - P_tau @ A_abs_ew
    id2_err = np.linalg.norm(comm_D + comm_A, 'fro')
    id2_rel = id2_err / np.linalg.norm(comm_D, 'fro') if np.linalg.norm(comm_D, 'fro') > 0 else id2_err
    print(f"\n  [D_sym, P_tau] = -[|A|_ew, P_tau]")
    print(f"    Relative error:  {id2_rel:.2e}")
    if id2_rel < 1e-13:
        print(f"    ✓ MACHINE ZERO — Identity 2 CONFIRMED")

    # ============================================================
    # CRT BLOCK CONSTRUCTION
    # ============================================================
    print(f"\n  ======= CRT block structure =======")

    # Determine block primes
    is_even = (m % 2 == 0)
    if is_even:
        block_primes = [5, 7, 11]
    else:
        v2 = build_V_diag(residues, 2)
        if v2.max() - v2.min() > 0.1:
            block_primes = [2, 5, 7, 11]
        else:
            block_primes = [5, 7, 11]
    print(f"  Block primes: {block_primes}")

    blocks = get_crt_blocks(residues, block_primes)
    n_blocks = len(blocks)
    dims = [len(v) for v in blocks.values()]
    block_dim = dims[0] if len(set(dims)) == 1 else "mixed"
    print(f"  {n_blocks} blocks × {block_dim}D = {sum(dims)} total")

    # ============================================================
    # PATH A (cont): Element-wise |A| block-diagonality
    # ============================================================
    print(f"\n  ======= PATH A: |A|_ew CRT block-diagonality =======")

    off, on, ratio = off_diagonal_block_norm(A_abs_ew, blocks)
    print(f"  ||A_abs_ew||² on-diagonal blocks:  {on:.4f}")
    print(f"  ||A_abs_ew||² off-diagonal blocks: {off:.4f}")
    print(f"  Off-diagonal fraction: {ratio:.6e}")

    if ratio < 1e-10:
        print(f"  ✓ |A|_ew EXACTLY block-diagonalizes in CRT basis!")
    elif ratio < 0.01:
        print(f"  ~ |A|_ew approximately block-diagonal ({ratio*100:.4f}% leakage)")
    else:
        print(f"  ✗ |A|_ew does NOT block-diagonalize ({ratio*100:.2f}% off-diagonal)")
        # Decompose: how much leaks between specific block pairs?
        block_list = list(blocks.values())
        # Find the hottest off-diagonal block pair
        max_off = 0
        max_pair = (0, 0)
        for i in range(n_blocks):
            for j in range(n_blocks):
                if i != j:
                    sub = A_abs_ew[np.ix_(block_list[i], block_list[j])]
                    nrm = np.sum(sub ** 2)
                    if nrm > max_off:
                        max_off = nrm
                        max_pair = (i, j)
        print(f"    Hottest off-diagonal pair: blocks ({max_pair[0]}, {max_pair[1]})")
        print(f"    That pair's fraction of total off-diagonal: {max_off/off*100:.1f}%")

    # Also test D_sym block-diagonality for reference
    off_d, on_d, ratio_d = off_diagonal_block_norm(D_sym, blocks)
    print(f"  D_sym off-diagonal fraction: {ratio_d:.6e}")

    # Commutator [|A|_ew, P_tau] block analysis
    print(f"\n  ======= Commutator [|A|_ew, P_tau] block structure =======")
    comm_Aew = A_abs_ew @ P_tau - P_tau @ A_abs_ew
    off_c, on_c, ratio_c = off_diagonal_block_norm(comm_Aew, blocks)
    print(f"  Off-diagonal fraction of [|A|_ew, P_tau]: {ratio_c:.4f}")
    # Fill ratio (fraction of entries > 1% of max)
    comm_max = np.abs(comm_Aew).max()
    if comm_max > 0:
        fill = np.sum(np.abs(comm_Aew) > 0.01 * comm_max) / comm_Aew.size
        print(f"  Fill ratio (>1% of max): {fill*100:.1f}%")

    # ============================================================
    # PATH B: FLAT-BAND AUTOPSY
    # ============================================================
    print(f"\n  ======= PATH B: H(α_c) eigendecomposition =======")
    alpha_c = sqrt(135 / 88)
    lambda_P = np.max(np.linalg.eigvalsh(D_sym))
    C = D_sym @ P_tau - P_tau @ D_sym
    H = D_sym / lambda_P + 1j * alpha_c * C / lambda_P
    evals, evecs = np.linalg.eigh(H)

    idx_sorted = np.argsort(evals)
    perron_eval = evals[idx_sorted[-1]]
    anti_eval = evals[idx_sorted[0]]
    print(f"  Perron:      {perron_eval:.6f}")
    print(f"  Anti-Perron: {anti_eval:.6f}")
    print(f"  Gap:         {perron_eval - anti_eval:.6f}")

    # Flat band: eigenvalues near zero
    flat_mask = np.abs(evals) < 0.05
    n_flat = np.sum(flat_mask)
    print(f"  Flat band (|λ|<0.05): {n_flat}/{phi}")

    # ---- Compute IPR for every eigenvector ----
    print(f"\n  ======= PATH B: IPR ranking =======")
    ipr_vals = np.array([ipr(evecs[:, j]) * phi for j in range(phi)])
    # Sort by IPR descending
    ipr_order = np.argsort(ipr_vals)[::-1]

    print(f"  Top 10 most localized eigenstates (IPR × φ):")
    print(f"  {'Rank':>4} {'EigIdx':>6} {'λ':>10} {'IPR×φ':>10} {'|λ|<0.05':>9} {'BPR':>6}")
    for rank in range(min(10, phi)):
        j = ipr_order[rank]
        is_flat = "FLAT" if np.abs(evals[j]) < 0.05 else ""
        bpr = block_participation_ratio(evecs[:, j], blocks)
        print(f"  {rank+1:4d} {j:6d} {evals[j]:10.6f} {ipr_vals[j]:10.4f} {is_flat:>9} {bpr:6.2f}")

    print(f"\n  Bottom 5 most delocalized (closest to 1.0):")
    for rank in range(max(0, phi - 5), phi):
        j = ipr_order[rank]
        is_flat = "FLAT" if np.abs(evals[j]) < 0.05 else ""
        bpr = block_participation_ratio(evecs[:, j], blocks)
        print(f"  {rank+1:4d} {j:6d} {evals[j]:10.6f} {ipr_vals[j]:10.4f} {is_flat:>9} {bpr:6.2f}")

    # ---- Flat-band IPR statistics ----
    flat_ipr = ipr_vals[flat_mask]
    nonflat_ipr = ipr_vals[~flat_mask]
    print(f"\n  Flat-band IPR×φ:    mean={flat_ipr.mean():.4f}, max={flat_ipr.max():.4f}, "
          f"std={flat_ipr.std():.4f}")
    print(f"  Non-flat-band IPR×φ: mean={nonflat_ipr.mean():.4f}, max={nonflat_ipr.max():.4f}")

    # Kurtosis of IPR distribution
    from scipy.stats import kurtosis as scipy_kurtosis
    kurt_flat = scipy_kurtosis(flat_ipr, fisher=True)
    kurt_all = scipy_kurtosis(ipr_vals, fisher=True)
    print(f"  Kurtosis (flat band only): {kurt_flat:.1f}")
    print(f"  Kurtosis (all eigenvectors): {kurt_all:.1f}")

    # ---- Block Participation Ratio analysis ----
    print(f"\n  ======= PATH B: Block Participation Ratio =======")

    # For the most localized flat-band states, where do they live?
    flat_indices = np.where(flat_mask)[0]
    flat_ipr_with_idx = [(ipr_vals[j], j) for j in flat_indices]
    flat_ipr_with_idx.sort(reverse=True)

    print(f"\n  Top 10 most localized FLAT-BAND states:")
    print(f"  {'Rank':>4} {'λ':>10} {'IPR×φ':>10} {'BPR':>6} {'BPR/N_blk':>9} {'Peak blk':>8} {'Peak%':>7}")
    for rank in range(min(10, len(flat_ipr_with_idx))):
        ipr_val, j = flat_ipr_with_idx[rank]
        bpr = block_participation_ratio(evecs[:, j], blocks)
        # Find which block has max weight
        block_list = list(blocks.values())
        block_weights = []
        for indices in block_list:
            w = np.sum(np.abs(evecs[indices, j]) ** 2)
            block_weights.append(w)
        block_weights = np.array(block_weights)
        peak_blk = np.argmax(block_weights)
        peak_pct = block_weights[peak_blk] * 100
        print(f"  {rank+1:4d} {evals[j]:10.6f} {ipr_val:10.4f} {bpr:6.2f} "
              f"{bpr/n_blocks:9.4f} {peak_blk:8d} {peak_pct:6.1f}%")

    # Average BPR for top-20 most localized flat-band states
    n_top = min(20, len(flat_ipr_with_idx))
    top_bprs = []
    for rank in range(n_top):
        _, j = flat_ipr_with_idx[rank]
        top_bprs.append(block_participation_ratio(evecs[:, j], blocks))
    avg_bpr = np.mean(top_bprs)
    print(f"\n  Avg BPR of top-{n_top} localized flat-band states: {avg_bpr:.2f} / {n_blocks}")
    print(f"  Ratio BPR/N_blocks: {avg_bpr/n_blocks:.4f}")

    # Outcome classification
    if avg_bpr < 3:
        print(f"  → OUTCOME 1: Extreme-IPR modes TRAPPED in single CRT block!")
        print(f"    Spatial defect aligns with CRT block boundary")
    elif avg_bpr / n_blocks > 0.7:
        print(f"  → OUTCOME 2: Extreme-IPR modes SPREAD across all CRT blocks")
        print(f"    Position-localized but CRT-spread — additive ≠ multiplicative")
    else:
        print(f"  → INTERMEDIATE: BPR = {avg_bpr:.1f} (neither trapped nor fully spread)")

    # ---- Position-space analysis of extreme modes ----
    print(f"\n  ======= PATH B: Position-space hotspot =======")
    # For the single most localized flat-band state, which residues carry the weight?
    if len(flat_ipr_with_idx) > 0:
        _, j_top = flat_ipr_with_idx[0]
        probs = np.abs(evecs[:, j_top]) ** 2
        top_positions = np.argsort(probs)[::-1][:5]
        print(f"  Most localized flat-band state (IPR×φ = {ipr_vals[j_top]:.4f}):")
        print(f"  Top 5 positions:")
        for rank, idx in enumerate(top_positions):
            r = residues[idx]
            p = probs[idx] * 100
            # Which block is this residue in?
            for blk_idx, (sig, indices) in enumerate(blocks.items()):
                if idx in indices:
                    break
            print(f"    r={r:5d} (block {blk_idx:2d}): {p:.2f}%")

    # ============================================================
    # COMBINED DIAGNOSTIC: Does |A|_ew WITHIN each block predict IPR?
    # ============================================================
    print(f"\n  ======= COMBINED: |A|_ew intra-block spectrum =======")
    # Compute eigenvalue spread of |A|_ew within each CRT block
    block_list = list(blocks.values())
    block_keys = list(blocks.keys())
    spreads = []
    for k in range(n_blocks):
        sub = A_abs_ew[np.ix_(block_list[k], block_list[k])]
        eigvals = np.linalg.eigvalsh(sub)
        spreads.append(eigvals.max() - eigvals.min())
    spreads = np.array(spreads)
    print(f"  |A|_ew intra-block spectral spread:")
    print(f"    Mean: {spreads.mean():.4f}")
    print(f"    Std:  {spreads.std():.4f}")
    print(f"    Min:  {spreads.min():.4f} (block {np.argmin(spreads)})")
    print(f"    Max:  {spreads.max():.4f} (block {np.argmax(spreads)})")
    # If all spreads are equal → |A|_ew is isospectral across blocks
    cv = spreads.std() / spreads.mean() if spreads.mean() > 0 else 0
    print(f"    CV:   {cv:.6f}")
    if cv < 0.01:
        print(f"    ✓ |A|_ew is ISOSPECTRAL across CRT blocks (CV < 1%)")
    else:
        print(f"    ~ |A|_ew spectrum varies across blocks (CV = {cv*100:.2f}%)")

    return {
        'phi': phi,
        'n_blocks': n_blocks,
        'block_dim': block_dim,
        'identity_err': identity_rel,
        'identity2_err': id2_rel,
        'A_abs_ew_off_diag': ratio,
        'D_sym_off_diag': ratio_d,
        'comm_off_diag': ratio_c,
        'kurtosis_all': kurt_all,
        'kurtosis_flat': kurt_flat,
        'avg_bpr_top20': avg_bpr,
        'n_flat': n_flat,
        'A_abs_ew_cv': cv,
    }


# ==============================================================
if __name__ == '__main__':
    results = {}
    results['even'] = run_experiment(2310, "EVEN PRIMORIAL")
    results['odd'] = run_experiment(1155, "ODD NON-PRIMORIAL")

    # ============================================================
    # SCORECARD
    # ============================================================
    print(f"\n{'='*70}")
    print(f"  SCORECARD: Corrected CRT-Chirality Bridge")
    print(f"{'='*70}")

    print(f"\n  ─── PATH A: Identity & Block Structure ───")
    for key, lab in [('even', 'm=2310'), ('odd', 'm=1155')]:
        r = results[key]
        print(f"\n  {lab} ({r['n_blocks']} blocks × {r['block_dim']}D):")
        id1 = "✓" if r['identity_err'] < 1e-13 else "✗"
        id2 = "✓" if r['identity2_err'] < 1e-13 else "✗"
        bd = "✓ EXACT" if r['A_abs_ew_off_diag'] < 1e-10 else \
             f"✗ {r['A_abs_ew_off_diag']*100:.2f}%"
        print(f"    Identity 1 (D_sym = m/2(J-I) - |A|_ew):    {id1} ({r['identity_err']:.1e})")
        print(f"    Identity 2 ([D_sym,P_tau] = -[|A|_ew,P_tau]): {id2} ({r['identity2_err']:.1e})")
        print(f"    |A|_ew block-diagonal in CRT:                {bd}")
        print(f"    D_sym block-diagonal in CRT:                 {r['D_sym_off_diag']*100:.2f}%")

    print(f"\n  ─── PATH B: Flat-Band Autopsy ───")
    for key, lab in [('even', 'm=2310'), ('odd', 'm=1155')]:
        r = results[key]
        print(f"\n  {lab}:")
        print(f"    Flat band: {r['n_flat']}/{r['phi']} states")
        print(f"    Kurtosis (all): {r['kurtosis_all']:.1f}")
        print(f"    Kurtosis (flat): {r['kurtosis_flat']:.1f}")
        print(f"    Avg BPR top-20: {r['avg_bpr_top20']:.2f} / {r['n_blocks']} blocks")
        print(f"    |A|_ew intra-block spectral CV: {r['A_abs_ew_cv']*100:.4f}%")

    print(f"\n  ─── KEY QUESTION ANSWERS ───")
    r_even = results['even']
    r_odd = results['odd']
    print(f"\n  Q1: Does element-wise |A| respect CRT blocks?")
    if r_even['A_abs_ew_off_diag'] < 1e-10 and r_odd['A_abs_ew_off_diag'] < 1e-10:
        print(f"    → YES! Additive topology (|A|_ew) preserves multiplicative blocks (CRT)")
        print(f"    → The chirality identity chain IS a CRT-native decomposition!")
    elif r_even['A_abs_ew_off_diag'] < 0.01 or r_odd['A_abs_ew_off_diag'] < 0.01:
        print(f"    → APPROXIMATELY — small leakage, not exact")
    else:
        print(f"    → NO — additive and multiplicative structures are genuinely orthogonal")
        print(f"    → CRT blocks and |A|_ew are independent decompositions")

    print(f"\n  Q2: Where do extreme-IPR modes live in CRT space?")
    if r_odd['avg_bpr_top20'] < 3:
        print(f"    → TRAPPED in single CRT blocks (Outcome 1)")
    elif r_odd['avg_bpr_top20'] / r_odd['n_blocks'] > 0.7:
        print(f"    → SPREAD across all blocks (Outcome 2)")
        print(f"    → Position-localized but CRT-delocalized")
    else:
        print(f"    → INTERMEDIATE — partial block trapping")

    print(f"\n  Q3: Is |A|_ew isospectral across CRT blocks?")
    if r_even['A_abs_ew_cv'] < 0.01 and r_odd['A_abs_ew_cv'] < 0.01:
        print(f"    → YES — every block sees the same |A|_ew spectrum")
    else:
        print(f"    → NO — spectral spread varies across blocks")

    print()
