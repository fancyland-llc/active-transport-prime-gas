"""
CRT–Chirality Bridge Test
=========================
Tests three predictions connecting the Chirality Identity Chain to CRT qudit architecture:

P1: |A| block-diagonalizes in the CRT basis (exact invariant subspaces)
P2: Flat band = N-1 passive CRT blocks; 1 active block contains Perron + anti-Perron
P3: At odd m, V_2 becomes a block label → 60 blocks × 8D (vs 30 × 16D at even m)

Additionally tests:
P4 (Gemini): [|A|, P_tau] is a sparse routing bus, not dense noise
P5: Kurtosis hotspots hoard into the active block at odd m

Runtime: < 2 min on consumer PC
"""

import numpy as np
from math import gcd, sqrt
from scipy.linalg import sqrtm


def get_coprime_residues(m):
    return sorted([r for r in range(1, m) if gcd(r, m) == 1])


def build_D_sym(residues, m):
    r = np.array(residues)
    diff = np.abs(r[:, None] - r[None, :])
    return np.minimum(diff, m - diff).astype(float)


def build_P_tau(residues, m):
    n = len(residues)
    res_to_idx = {r: i for i, r in enumerate(residues)}
    perm = np.array([res_to_idx[pow(r, -1, m)] for r in residues])
    P = np.zeros((n, n))
    P[np.arange(n), perm] = 1.0
    return P


def build_A_abs(D_sym):
    """Compute |A| = sqrt(A^T A) where A = (D_directed - D_directed^T)/2.
    For symmetric D_sym, use the relation: the antisymmetric part A of the
    directed distance matrix satisfies D_sym = m/2*(J-I) - |A|,
    so |A| = m/2*(J-I) - D_sym."""
    # Actually we need the directed D first, then A = (D - D^T)/2, then |A| = sqrtm(A.T @ A)
    # But from Identity 2.5: |A| = m/2*(J-I) - D_sym
    # Let's verify both ways
    n = D_sym.shape[0]
    # Method via identity (fast):
    J = np.ones((n, n))
    I = np.eye(n)
    # We need m - but it's implicit in D_sym. Extract from max entry.
    # Actually let's compute m from D_sym: max possible distance is m/2
    m_half = D_sym.max()
    A_abs_identity = m_half * (J - I) - D_sym
    return A_abs_identity, m_half * 2


def build_A_abs_exact(residues, m):
    """Build |A| via sqrtm(A^T A) for the actual antisymmetric matrix."""
    n = len(residues)
    r = np.array(residues)
    # Directed distance: d(r_i, r_j) = (r_j - r_i) mod m
    D_dir = ((r[None, :] - r[:, None]) % m).astype(float)
    A = (D_dir - D_dir.T) / 2.0
    # |A| = sqrtm(A^T @ A) = sqrtm(A.T @ A)  (A is real skew-symmetric, so A^T = -A)
    # A^T A = (-A)(A) = -A^2, which is positive semidefinite
    ATA = A.T @ A  # = -A^2
    A_abs = np.real(sqrtm(ATA))
    return A_abs, A


def build_V_diag(residues, p):
    return np.array([np.cos(2 * np.pi * r / p) for r in residues])


def get_crt_blocks(residues, block_primes):
    """Group residues by joint cosine signature."""
    blocks = {}
    for i, r in enumerate(residues):
        sig = tuple(round(np.cos(2 * np.pi * r / p), 8) for p in block_primes)
        blocks.setdefault(sig, []).append(i)
    return blocks


def project_into_block(M, indices):
    """Extract the block submatrix M[indices][:, indices]."""
    return M[np.ix_(indices, indices)]


def off_diagonal_block_norm(M, blocks):
    """Compute sum of squared Frobenius norms of off-diagonal blocks."""
    block_list = list(blocks.values())
    total_off = 0.0
    total_on = 0.0
    for i, bi in enumerate(block_list):
        for j, bj in enumerate(block_list):
            sub = M[np.ix_(bi, bj)]
            norm_sq = np.sum(sub ** 2)
            if i == j:
                total_on += norm_sq
            else:
                total_off += norm_sq
    total = total_on + total_off
    return total_off, total_on, total_off / total if total > 0 else 0.0


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
    P_tau = build_P_tau(residues, m)

    # Build |A| via identity (fast)
    A_abs_id, m_recovered = build_A_abs(D_sym)
    print(f"  m recovered from D_sym max: {m_recovered:.0f} (expected {m})")

    # Build |A| via sqrtm (exact, for verification at this size)
    print(f"  Computing |A| via sqrtm (n={phi})...")
    A_abs_exact, A_skew = build_A_abs_exact(residues, m)

    # Verify identity
    identity_err = np.linalg.norm(A_abs_id - A_abs_exact, 'fro') / np.linalg.norm(A_abs_exact, 'fro')
    print(f"  Identity 2.5 verification: ||A_abs_id - A_abs_exact||/||A_abs_exact|| = {identity_err:.2e}")

    # ---- V_p analysis ----
    print(f"\n  --- V_p scalar test ---")
    for p in [2, 3, 5, 7, 11, 13]:
        v = build_V_diag(residues, p)
        spread = v.max() - v.min()
        n_distinct = len(set(np.round(v, 8)))
        print(f"  V_{p:2d}: spread = {spread:.6f}, distinct values = {n_distinct}"
              + (f"  → SCALAR ({v[0]:.4f})" if spread < 1e-10 else ""))

    # ---- Determine block primes ----
    # At even m: block primes = in-factorization primes > 3
    # At odd m: V_2 is NOT scalar → add it to block labels
    is_even = (m % 2 == 0)

    if is_even:
        block_primes = [5, 7, 11]
        print(f"\n  Even m → block primes: {block_primes}")
    else:
        v2 = build_V_diag(residues, 2)
        v2_spread = v2.max() - v2.min()
        if v2_spread > 0.1:
            block_primes = [2, 5, 7, 11]
            print(f"\n  Odd m → V_2 is non-trivial (spread={v2_spread:.4f}) → block primes: {block_primes}")
        else:
            block_primes = [5, 7, 11]
            print(f"\n  Odd m but V_2 scalar → block primes: {block_primes}")

    blocks = get_crt_blocks(residues, block_primes)
    n_blocks = len(blocks)
    dims = [len(v) for v in blocks.values()]
    print(f"  CRT blocks: {n_blocks} blocks")
    print(f"  Block dimensions: min={min(dims)}, max={max(dims)}, "
          f"all equal={len(set(dims))==1}")
    if len(set(dims)) == 1:
        print(f"  → {n_blocks} blocks × {dims[0]}D = {n_blocks * dims[0]} total")

    # ==== PREDICTION 1: |A| block-diagonalizes in CRT basis ====
    print(f"\n  ======= P1: |A| block-diagonalization =======")
    off, on, ratio = off_diagonal_block_norm(A_abs_exact, blocks)
    print(f"  ||A_abs||² on-diagonal blocks:  {on:.4f}")
    print(f"  ||A_abs||² off-diagonal blocks: {off:.4f}")
    print(f"  Off-diagonal fraction: {ratio:.2e}")
    if ratio < 1e-10:
        print(f"  ✓ |A| EXACTLY block-diagonalizes in CRT basis!")
    elif ratio < 0.01:
        print(f"  ~ |A| APPROXIMATELY block-diagonal (< 1% leakage)")
    else:
        print(f"  ✗ |A| does NOT block-diagonalize ({ratio*100:.1f}% off-diagonal)")

    # Also test D_sym
    off_d, on_d, ratio_d = off_diagonal_block_norm(D_sym, blocks)
    print(f"  D_sym off-diagonal fraction: {ratio_d:.2e}")

    # ==== PREDICTION 4 (Gemini): [|A|, P_tau] sparsity in block basis ====
    print(f"\n  ======= P4: Commutator [|A|, P_tau] block sparsity =======")
    comm = A_abs_exact @ P_tau - P_tau @ A_abs_exact
    off_c, on_c, ratio_c = off_diagonal_block_norm(comm, blocks)
    print(f"  ||[|A|, P_tau]||² on-diagonal:  {on_c:.4f}")
    print(f"  ||[|A|, P_tau]||² off-diagonal: {off_c:.4f}")
    print(f"  Off-diagonal fraction: {ratio_c:.4f}")

    # Detailed block-to-block routing map
    block_list = list(blocks.values())
    n_b = len(block_list)
    routing = np.zeros((n_b, n_b))
    for i in range(n_b):
        for j in range(n_b):
            sub = comm[np.ix_(block_list[i], block_list[j])]
            routing[i, j] = np.linalg.norm(sub, 'fro')

    # Sparsity: fraction of block-pairs with negligible coupling
    threshold = routing.max() * 0.01
    n_pairs = n_b * (n_b - 1)
    n_zero = sum(1 for i in range(n_b) for j in range(n_b)
                 if i != j and routing[i, j] < threshold)
    sparsity = n_zero / n_pairs if n_pairs > 0 else 0
    print(f"  Block routing matrix: {n_b}×{n_b}")
    print(f"  Max coupling: {routing.max():.4f}")
    print(f"  Structural zeros (< 1% of max): {n_zero}/{n_pairs} = {sparsity*100:.1f}%")
    if sparsity > 0.5:
        print(f"  ✓ Sparse routing bus — deterministic topology, not dense noise")
    else:
        print(f"  ~ Dense coupling — commutator mixes most blocks")

    # ==== EIGENVECTOR ANALYSIS ====
    print(f"\n  ======= Eigendecomposition of H(alpha_c) =======")
    alpha_c = sqrt(135 / 88)
    lambda_P = np.max(np.linalg.eigvalsh(D_sym))
    C = D_sym @ P_tau - P_tau @ D_sym
    H = D_sym / lambda_P + 1j * alpha_c * C / lambda_P
    evals, evecs = np.linalg.eigh(H)

    # Identify Perron (max eigenvalue), anti-Perron (min), flat band
    idx_sorted = np.argsort(evals)
    idx_perron = idx_sorted[-1]
    idx_anti = idx_sorted[0]
    perron_eval = evals[idx_perron]
    anti_eval = evals[idx_anti]
    print(f"  Perron eigenvalue:      {perron_eval:.6f}")
    print(f"  Anti-Perron eigenvalue:  {anti_eval:.6f}")
    print(f"  Gap: {perron_eval - anti_eval:.6f}")

    # Flat band: all eigenvalues with |λ| < 0.01
    flat_mask = np.abs(evals) < 0.05
    n_flat = np.sum(flat_mask)
    print(f"  Flat band count (|λ| < 0.05): {n_flat} / {phi}")

    # ==== PREDICTION 2: Flat band = N-1 passive blocks ====
    print(f"\n  ======= P2: Flat band / CRT block alignment =======")
    flat_indices = np.where(flat_mask)[0]
    flat_evecs = evecs[:, flat_indices]

    # For each CRT block, compute overlap with flat-band subspace
    # Overlap = sum of |<block_basis_i | flat_evec_j>|^2 over all i in block, j in flat band
    # Equivalently: project each block's position-basis vectors into flat-band subspace
    block_flat_overlap = []
    block_perron_overlap = []
    block_anti_overlap = []

    v_perron = evecs[:, idx_perron]
    v_anti = evecs[:, idx_anti]

    for k, (sig, indices) in enumerate(blocks.items()):
        dim = len(indices)
        # Block's contribution to flat band
        # Sum of |v_flat_j[i]|^2 for i in block, j in flat band
        flat_mass = np.sum(np.abs(flat_evecs[indices, :]) ** 2)
        # Normalize by block dimension (each block "owns" dim basis vectors)
        flat_frac = flat_mass / dim

        perron_mass = np.sum(np.abs(v_perron[indices]) ** 2)
        anti_mass = np.sum(np.abs(v_anti[indices]) ** 2)

        block_flat_overlap.append(flat_frac)
        block_perron_overlap.append(perron_mass)
        block_anti_overlap.append(anti_mass)

    # Sort by flat-band overlap
    order = np.argsort(block_flat_overlap)
    print(f"\n  Block-by-block analysis (sorted by flat-band overlap):")
    print(f"  {'Block':>5} {'Dim':>4} {'Flat%':>8} {'Perron%':>9} {'Anti-P%':>9}")

    block_keys = list(blocks.keys())
    n_passive = 0
    active_idx = -1
    for rank, k in enumerate(order):
        ff = block_flat_overlap[k] * 100
        pf = block_perron_overlap[k] * 100
        af = block_anti_overlap[k] * 100
        dim = len(blocks[block_keys[k]])
        tag = ""
        if ff < 90:
            tag = " ← ACTIVE"
            active_idx = k
        else:
            n_passive += 1
        if rank < 3 or rank >= n_blocks - 3 or tag:
            print(f"  {k:5d} {dim:4d} {ff:7.2f}% {pf:8.4f}% {af:8.4f}%{tag}")
        elif rank == 3:
            print(f"  {'...':>5}")

    print(f"\n  Passive blocks (>90% flat band): {n_passive} / {n_blocks}")
    print(f"  Active blocks (<90% flat band):  {n_blocks - n_passive}")

    if n_passive == n_blocks - 1:
        print(f"  ✓ Exactly 1 active block — flat band IS the passive CRT blocks!")
    elif n_passive == n_blocks:
        print(f"  ~ ALL blocks are passive — Perron/anti-Perron spread across blocks")
        # Check if Perron is distributed equally
        po = np.array(block_perron_overlap)
        print(f"    Perron mass: max block = {po.max()*100:.2f}%, "
              f"min = {po.min()*100:.2f}%, ratio = {po.max()/po.min():.2f}×")
    else:
        print(f"  ~ Multiple active blocks ({n_blocks - n_passive})")

    # ==== PREDICTION 5: Kurtosis per block ====
    print(f"\n  ======= P5: Kurtosis per CRT block =======")
    # Compute IPR and kurtosis for anti-Perron eigenvector per block
    v = np.abs(v_anti) ** 2  # probability distribution
    # Per-block kurtosis of the anti-Perron weight distribution
    for k_idx in [order[0], order[-1]]:  # least and most flat-band overlap
        indices = blocks[block_keys[k_idx]]
        w = v[indices]
        w_norm = w / w.sum() if w.sum() > 0 else w
        mean_w = np.mean(w_norm)
        var_w = np.var(w_norm)
        if var_w > 0:
            kurt = np.mean((w_norm - mean_w)**4) / var_w**2 - 3
        else:
            kurt = 0
        ipr = np.sum(w_norm**2)
        print(f"  Block {k_idx:3d} (dim={len(indices):2d}): "
              f"anti-Perron mass={np.sum(w[:])*100:.4f}%, "
              f"within-block IPR={ipr:.4f}, kurtosis={kurt:.1f}")

    # Global eigenvector kurtosis (all eigenvectors)
    all_ipr = []
    for j in range(phi):
        v_j = np.abs(evecs[:, j]) ** 2
        all_ipr.append(np.sum(v_j ** 2) * phi)
    all_ipr = np.array(all_ipr)
    from scipy.stats import kurtosis as scipy_kurtosis
    try:
        kurt_all = scipy_kurtosis(all_ipr, fisher=True)
    except ImportError:
        kurt_all = np.mean((all_ipr - all_ipr.mean())**4) / np.var(all_ipr)**2 - 3
    print(f"  Global IPR×φ: mean={all_ipr.mean():.4f}, max={all_ipr.max():.4f}, kurtosis={kurt_all:.1f}")

    return {
        'phi': phi,
        'n_blocks': n_blocks,
        'block_dim': dims[0] if len(set(dims)) == 1 else dims,
        'A_abs_off_diag_ratio': ratio,
        'comm_off_diag_ratio': ratio_c,
        'comm_sparsity': sparsity,
        'n_passive': n_passive,
        'identity_err': identity_err,
        'kurtosis': kurt_all,
    }


if __name__ == '__main__':
    results = {}

    # Even primorial
    results['even'] = run_experiment(2310, "EVEN PRIMORIAL")

    # Odd matched pair
    results['odd'] = run_experiment(1155, "ODD NON-PRIMORIAL")

    # Summary
    print(f"\n{'='*70}")
    print(f"  SUMMARY: Prediction Scorecard")
    print(f"{'='*70}")

    print(f"\n  P1: |A| block-diagonalizes in CRT basis")
    print(f"    Even m=2310: off-diag fraction = {results['even']['A_abs_off_diag_ratio']:.2e}"
          + (" ✓" if results['even']['A_abs_off_diag_ratio'] < 1e-10 else " ✗"))
    print(f"    Odd  m=1155: off-diag fraction = {results['odd']['A_abs_off_diag_ratio']:.2e}"
          + (" ✓" if results['odd']['A_abs_off_diag_ratio'] < 1e-10 else " ✗"))

    print(f"\n  P2: Flat band = N-1 passive CRT blocks")
    print(f"    Even: {results['even']['n_passive']}/{results['even']['n_blocks']} passive"
          + (" ✓" if results['even']['n_passive'] == results['even']['n_blocks'] - 1 else " ?"))
    print(f"    Odd:  {results['odd']['n_passive']}/{results['odd']['n_blocks']} passive"
          + (" ✓" if results['odd']['n_passive'] == results['odd']['n_blocks'] - 1 else " ?"))

    print(f"\n  P3: Block structure at odd m")
    print(f"    Even m=2310: {results['even']['n_blocks']} blocks × {results['even']['block_dim']}D")
    print(f"    Odd  m=1155: {results['odd']['n_blocks']} blocks × {results['odd']['block_dim']}D")
    expected_odd = 60
    print(f"    Expected odd: {expected_odd} blocks × 8D"
          + (" ✓" if results['odd']['n_blocks'] == expected_odd else " ✗"))

    print(f"\n  P4: [|A|, P_tau] routing sparsity")
    print(f"    Even: {results['even']['comm_sparsity']*100:.1f}% structural zeros"
          + (" ✓ sparse" if results['even']['comm_sparsity'] > 0.5 else " ~ dense"))
    print(f"    Odd:  {results['odd']['comm_sparsity']*100:.1f}% structural zeros"
          + (" ✓ sparse" if results['odd']['comm_sparsity'] > 0.5 else " ~ dense"))

    print(f"\n  P5: Kurtosis")
    print(f"    Even: {results['even']['kurtosis']:.1f}")
    print(f"    Odd:  {results['odd']['kurtosis']:.1f}")

    print(f"\n  Identity 2.5 verification:")
    print(f"    Even: {results['even']['identity_err']:.2e}")
    print(f"    Odd:  {results['odd']['identity_err']:.2e}")
