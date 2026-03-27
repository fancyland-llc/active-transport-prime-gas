"""
Orbit Hamiltonian — Arithmetic Crystal Band Structure
=====================================================
Projects H(α_c) from the φ-dimensional residue space into the 
Klein-4 orbit space and cross-tabulates orbits with CRT blocks.

Tests:
  1. NORM CAPTURE: What fraction of ||H||²_F lives in orbit space?
     - Orbit projector P_orb averages within each orbit
     - Complement (I - P_orb) captures intra-orbit fluctuations
  2. ORBIT HAMILTONIAN SPECTRUM: Does H_tilde = P_orb H P_orb simplify?
     - Block-diagonal by CRT grouping? (Brillouin zone test)
     - Eigenvalue structure — does it reproduce the Perron gap?
  3. O(2) COUPLING: Is O(2) bound state or resonance?
     - Diagonal vs off-diagonal power in orbit Hamiltonian row
  4. CROSS-TABULATION: Where do Klein-4 orbits land in CRT blocks?
     - How many CRT blocks does each orbit span?
     - Structure factor: is there a pattern?
  5. IRREP DECOMPOSITION: Project onto the 4 Klein-4 irreps within each orbit
     - (++), (+-), (-+), (--) parity sectors
     - Which sector carries the physics?
  6. FRUSTRATED CRYSTAL TEST: Measure the CRT-orbit commutator mismatch

Runtime: < 60s
"""

import numpy as np
from math import gcd, sqrt


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


def build_klein4_orbits(residues, m):
    """Partition coprime residues into Klein-4 orbits under <sigma, tau>."""
    res_set = set(residues)
    res_to_idx = {r: i for i, r in enumerate(residues)}
    visited = set()
    orbits = []

    for r in residues:
        if r in visited:
            continue
        r_inv = pow(r, -1, m)
        mr = m - r
        mr_inv = m - r_inv
        orbit_residues = sorted({x for x in [r, mr, r_inv, mr_inv] if x in res_set})
        orbit_indices = [res_to_idx[x] for x in orbit_residues]
        orbits.append({
            'residues': orbit_residues,
            'indices': orbit_indices,
            'size': len(orbit_residues),
            'seed': min(orbit_residues),
        })
        visited.update(orbit_residues)

    return orbits


def get_crt_blocks(residues, block_primes):
    blocks = {}
    for i, r in enumerate(residues):
        sig = tuple(round(np.cos(2 * np.pi * r / p), 8) for p in block_primes)
        blocks.setdefault(sig, []).append(i)
    return blocks


def build_orbit_projector(orbits, phi):
    """Build the orbit-averaging projector P_orb.
    P_orb = sum_k |O_k><O_k| where |O_k> = (1/sqrt(|O_k|)) sum_{r in O_k} |r>.
    P_orb projects onto the symmetric subspace of each orbit."""
    P = np.zeros((phi, phi))
    for orb in orbits:
        idx = orb['indices']
        sz = orb['size']
        for i in idx:
            for j in idx:
                P[i, j] += 1.0 / sz
    return P


def build_orbit_basis(orbits, phi):
    """Build the N_orb × phi orbit basis matrix U.
    Row k = |O_k> = (1/sqrt(|O_k|)) * indicator on orbit k."""
    n_orb = len(orbits)
    U = np.zeros((n_orb, phi))
    for k, orb in enumerate(orbits):
        for i in orb['indices']:
            U[k, i] = 1.0 / np.sqrt(orb['size'])
    return U


def build_irrep_projectors(orbits, residues, m, phi):
    """Build the 4 Klein-4 irrep projectors for each orbit.
    For a size-4 orbit {r, m-r, r^{-1}, (m-r)^{-1}}:
      |++> = (|r> + |mr> + |ri> + |mri>) / 2
      |+-> = (|r> + |mr> - |ri> - |mri>) / 2
      |-+> = (|r> - |mr> + |ri> - |mri>) / 2
      |--> = (|r> - |mr> - |ri> + |mri>) / 2
    For size-2 orbits, only 2 irreps are available.
    Returns a dict of 4 projectors (phi × phi matrices).
    """
    res_to_idx = {r: i for i, r in enumerate(residues)}
    # Signs: (sigma_sign, tau_sign) for each irrep
    irrep_labels = ['++', '+-', '-+', '--']
    projectors = {label: np.zeros((phi, phi)) for label in irrep_labels}
    # Also build the full irrep basis vectors for the orbit Hamiltonian in irrep space
    irrep_basis = {label: [] for label in irrep_labels}

    for orb in orbits:
        res = orb['residues']
        seed = orb['seed']
        sz = orb['size']

        if sz == 4:
            # Identify the four elements
            r = seed
            r_inv = pow(r, -1, m)
            mr = m - r
            mr_inv = m - r_inv
            idx_r = res_to_idx[r]
            idx_mr = res_to_idx[mr]
            idx_ri = res_to_idx[r_inv]
            idx_mri = res_to_idx[mr_inv]

            # sigma(r) = m-r, tau(r) = r^{-1}
            # sigma sign: +1 for even, -1 for odd
            # tau sign: +1 for even, -1 for odd
            vectors = {
                '++': np.array([1, 1, 1, 1]) / 2.0,
                '+-': np.array([1, 1, -1, -1]) / 2.0,
                '-+': np.array([1, -1, 1, -1]) / 2.0,
                '--': np.array([1, -1, -1, 1]) / 2.0,
            }
            indices = [idx_r, idx_mr, idx_ri, idx_mri]

            for label in irrep_labels:
                v = np.zeros(phi)
                for k, idx in enumerate(indices):
                    v[idx] = vectors[label][k]
                projectors[label] += np.outer(v, v)
                irrep_basis[label].append(v)

        elif sz == 2:
            # Only sigma or tau collapses the orbit to size 2
            idx0 = orb['indices'][0]
            idx1 = orb['indices'][1]
            r0 = res[0]
            r1 = res[1]

            # Determine which symmetry fixes: r0 + r1 = m means sigma-fixed
            if r0 + r1 == m:
                # sigma(r) = r (i.e., r = m-r mod m, impossible for coprime)
                # Actually r0 = r, r1 = m-r, and r = r^{-1} (so tau also fixes)
                pass
            # Just do symmetric/antisymmetric combination
            v_sym = np.zeros(phi)
            v_sym[idx0] = 1.0 / np.sqrt(2)
            v_sym[idx1] = 1.0 / np.sqrt(2)
            v_asym = np.zeros(phi)
            v_asym[idx0] = 1.0 / np.sqrt(2)
            v_asym[idx1] = -1.0 / np.sqrt(2)

            projectors['++'] += np.outer(v_sym, v_sym)
            projectors['--'] += np.outer(v_asym, v_asym)
            irrep_basis['++'].append(v_sym)
            irrep_basis['--'].append(v_asym)

    return projectors, irrep_basis


def run_experiment(m, label):
    print(f"\n{'='*70}")
    print(f"  {label}: m = {m}")
    print(f"{'='*70}")

    residues = get_coprime_residues(m)
    phi = len(residues)
    res_to_idx = {r: i for i, r in enumerate(residues)}
    print(f"  phi = {phi}")

    # Build operators
    D_sym = build_D_sym(residues, m)
    P_tau = build_P_tau(residues, m)

    # Build H(alpha_c)
    alpha_c = sqrt(135 / 88)
    lambda_P = np.max(np.linalg.eigvalsh(D_sym))
    C = D_sym @ P_tau - P_tau @ D_sym
    H = D_sym / lambda_P + 1j * alpha_c * C / lambda_P

    H_norm_sq = np.sum(np.abs(H) ** 2)  # ||H||²_F
    print(f"  ||H||²_F = {H_norm_sq:.4f}")

    # Build Klein-4 orbits
    orbits = build_klein4_orbits(residues, m)
    n_orbits = len(orbits)
    size_counts = {}
    for o in orbits:
        size_counts[o['size']] = size_counts.get(o['size'], 0) + 1
    print(f"  Orbits: {n_orbits} ({dict(sorted(size_counts.items()))})")

    o2_exists = any(o['seed'] == 2 for o in orbits)
    print(f"  O(2) exists: {o2_exists}")

    # ============================================================
    # 1. NORM CAPTURE: orbit projection
    # ============================================================
    print(f"\n  ======= 1. NORM CAPTURE =======")

    U = build_orbit_basis(orbits, phi)  # n_orb × phi
    # Orbit Hamiltonian: H_tilde = U H U^T  (n_orb × n_orb)
    H_tilde = U @ H @ U.T

    H_tilde_norm_sq = np.sum(np.abs(H_tilde) ** 2)
    capture_frac = H_tilde_norm_sq / H_norm_sq
    print(f"  ||H_tilde||²_F = {H_tilde_norm_sq:.4f}")
    print(f"  Orbit-space capture: {capture_frac*100:.1f}%  ({H_tilde_norm_sq:.4f} / {H_norm_sq:.4f})")
    print(f"  Intra-orbit complement: {(1-capture_frac)*100:.1f}%")

    # Also check: what does the complement carry?
    P_orb = build_orbit_projector(orbits, phi)
    H_complement = H - P_orb @ H @ P_orb
    H_comp_norm_sq = np.sum(np.abs(H_complement) ** 2)
    # Note: ||H||² = ||P H P||² + ||H - PHP||² only if P is orthogonal projector
    # actually ||H||² ≠ ||PHP||² + ||H-PHP||² in general; let's check
    cross_term = H_norm_sq - H_tilde_norm_sq - H_comp_norm_sq
    print(f"  Cross-term (should be ~0 for orthogonal decomp): {cross_term:.4f}")

    # Better: use Parseval-like decomposition
    # Project H into orbit-symmetric and orbit-antisymmetric parts
    H_sym = P_orb @ H @ P_orb
    H_sym_norm = np.sum(np.abs(H_sym) ** 2)
    print(f"  P_orb H P_orb norm²: {H_sym_norm:.4f} ({H_sym_norm/H_norm_sq*100:.1f}%)")

    # ============================================================
    # 2. ORBIT HAMILTONIAN SPECTRUM
    # ============================================================
    print(f"\n  ======= 2. ORBIT HAMILTONIAN SPECTRUM =======")

    evals_tilde = np.linalg.eigvalsh(H_tilde)
    evals_full, evecs_full = np.linalg.eigh(H)

    print(f"  H_tilde eigenvalues ({n_orbits} total):")
    print(f"    Perron:      {evals_tilde[-1]:.6f}  (full H: {evals_full[-1]:.6f})")
    print(f"    Anti-Perron: {evals_tilde[0]:.6f}  (full H: {evals_full[0]:.6f})")
    print(f"    Gap:         {evals_tilde[-1] - evals_tilde[0]:.6f}"
          f"  (full H: {evals_full[-1] - evals_full[0]:.6f})")

    # Flat band in orbit space
    flat_tilde = np.sum(np.abs(evals_tilde) < 0.05)
    flat_full = np.sum(np.abs(evals_full) < 0.05)
    print(f"    Flat band (|λ|<0.05): {flat_tilde}/{n_orbits}"
          f"  (full H: {flat_full}/{phi})")

    # Eigenvalue histogram comparison
    print(f"\n  Eigenvalue distribution comparison:")
    for threshold in [0.01, 0.05, 0.1, 0.5]:
        n_t = np.sum(np.abs(evals_tilde) < threshold)
        n_f = np.sum(np.abs(evals_full) < threshold)
        print(f"    |λ| < {threshold}: orbit = {n_t}/{n_orbits},"
              f"  full = {n_f}/{phi}")

    # ============================================================
    # 3. O(2) COUPLING IN ORBIT SPACE
    # ============================================================
    print(f"\n  ======= 3. O(2) COUPLING =======")

    if o2_exists:
        o2_idx = next(k for k, o in enumerate(orbits) if o['seed'] == 2)
        o2_diag = np.abs(H_tilde[o2_idx, o2_idx]) ** 2
        o2_row = np.sum(np.abs(H_tilde[o2_idx, :]) ** 2)
        o2_off = o2_row - o2_diag
        print(f"  O(2) orbit index: {o2_idx}")
        print(f"  |H_tilde[O2,O2]|²: {o2_diag:.6f}")
        print(f"  |H_tilde[O2,:]|²:  {o2_row:.6f}")
        print(f"  Off-diagonal power: {o2_off:.6f} ({o2_off/o2_row*100:.1f}%)")
        if o2_off / o2_row > 0.9:
            print(f"  → RESONANCE: O(2) strongly coupled to other orbits")
        elif o2_off / o2_row < 0.1:
            print(f"  → BOUND STATE: O(2) nearly decoupled")
        else:
            print(f"  → INTERMEDIATE: partial coupling")

        # Top 5 orbits that O(2) couples to
        o2_couplings = np.abs(H_tilde[o2_idx, :]) ** 2
        o2_couplings[o2_idx] = 0  # exclude self
        top_coupled = np.argsort(o2_couplings)[::-1][:5]
        print(f"  Top 5 orbits coupled to O(2):")
        for k in top_coupled:
            if o2_couplings[k] > 1e-10:
                print(f"    O({orbits[k]['seed']:4d}): |H_tilde[O2,Ok]|² = {o2_couplings[k]:.6f}"
                      f" ({orbits[k]['residues']})")
    else:
        print(f"  O(2) does not exist at this modulus")
        # Check the orbit closest to the boundary instead
        min_seed = min(o['seed'] for o in orbits)
        boundary_idx = next(k for k, o in enumerate(orbits) if o['seed'] == min_seed)
        b_diag = np.abs(H_tilde[boundary_idx, boundary_idx]) ** 2
        b_row = np.sum(np.abs(H_tilde[boundary_idx, :]) ** 2)
        b_off = b_row - b_diag
        print(f"  Boundary orbit O({min_seed}): off-diag = {b_off/b_row*100:.1f}%")

    # ============================================================
    # 4. CROSS-TABULATION: Orbits vs CRT Blocks
    # ============================================================
    print(f"\n  ======= 4. ORBIT × CRT CROSS-TABULATION =======")

    is_even = (m % 2 == 0)
    if is_even:
        block_primes = [5, 7, 11]
    else:
        v2 = np.array([np.cos(2 * np.pi * r / 2) for r in residues])
        if v2.max() - v2.min() > 0.1:
            block_primes = [2, 5, 7, 11]
        else:
            block_primes = [5, 7, 11]
    print(f"  CRT block primes: {block_primes}")

    blocks = get_crt_blocks(residues, block_primes)
    n_blocks = len(blocks)
    block_list = list(blocks.values())
    print(f"  CRT blocks: {n_blocks}")

    # Map each residue index to its CRT block
    idx_to_block = {}
    for b_idx, indices in enumerate(block_list):
        for i in indices:
            idx_to_block[i] = b_idx

    # For each orbit, count how many CRT blocks its members span
    orbit_block_spans = []
    orbits_in_single_block = 0
    for orb in orbits:
        blocks_hit = set(idx_to_block[i] for i in orb['indices'])
        span = len(blocks_hit)
        orbit_block_spans.append(span)
        if span == 1:
            orbits_in_single_block += 1

    span_counts = {}
    for s in orbit_block_spans:
        span_counts[s] = span_counts.get(s, 0) + 1

    print(f"  Orbits contained in single CRT block: {orbits_in_single_block}/{n_orbits}")
    print(f"  Span distribution (# CRT blocks per orbit): {dict(sorted(span_counts.items()))}")

    if orbits_in_single_block == 0:
        print(f"  → EVERY orbit scatters across multiple CRT blocks")
        print(f"  → Klein-4 orbits are MAXIMALLY non-local in CRT space")
    elif orbits_in_single_block == n_orbits:
        print(f"  → Every orbit fits inside one CRT block — clean tensor product!")
    else:
        print(f"  → Mixed: {orbits_in_single_block} fit, {n_orbits - orbits_in_single_block} scatter")

    # Detailed span for O(2) and boundary orbits
    if o2_exists:
        o2_orb = next(o for o in orbits if o['seed'] == 2)
        o2_blocks = set(idx_to_block[i] for i in o2_orb['indices'])
        print(f"  O(2) = {o2_orb['residues']} spans {len(o2_blocks)} CRT blocks")

    # Show first 10 orbits with their CRT block distribution
    print(f"\n  First 10 orbits: residues → CRT blocks")
    for orb in orbits[:10]:
        blocks_hit = [idx_to_block[i] for i in orb['indices']]
        unique_blocks = sorted(set(blocks_hit))
        print(f"    O({orb['seed']:4d}) [{orb['size']}]: "
              f"residues {orb['residues']} → blocks {unique_blocks} "
              f"(span {len(unique_blocks)})")

    # ============================================================
    # 5. H_TILDE BLOCK STRUCTURE (does orbit H respect CRT grouping?)
    # ============================================================
    print(f"\n  ======= 5. ORBIT H vs CRT BLOCK STRUCTURE =======")

    # Group orbits by their CRT block span pattern
    # For each orbit, determine which CRT blocks it touches
    # Then group the orbit indices by block membership to see if H_tilde
    # has any block-diagonal structure when orbits are sorted by CRT affinity

    # Build a CRT-affinity grouping: group orbits by which set of CRT blocks they span
    orbit_crt_groups = {}
    for k, orb in enumerate(orbits):
        blocks_hit = frozenset(idx_to_block[i] for i in orb['indices'])
        orbit_crt_groups.setdefault(blocks_hit, []).append(k)

    n_crt_groups = len(orbit_crt_groups)
    print(f"  Distinct CRT block patterns: {n_crt_groups}")
    crt_group_sizes = sorted([len(v) for v in orbit_crt_groups.values()], reverse=True)
    print(f"  Group sizes (top 10): {crt_group_sizes[:10]}")

    # Compute block-diagonality of H_tilde when orbits are grouped by CRT pattern
    total_on_crt = 0.0
    total_off_crt = 0.0
    groups_list = list(orbit_crt_groups.values())
    for i, gi in enumerate(groups_list):
        for j, gj in enumerate(groups_list):
            sub = H_tilde[np.ix_(gi, gj)]
            norm_sq = np.sum(np.abs(sub) ** 2)
            if i == j:
                total_on_crt += norm_sq
            else:
                total_off_crt += norm_sq

    total_crt = total_on_crt + total_off_crt
    off_frac_crt = total_off_crt / total_crt if total_crt > 0 else 0
    print(f"  H_tilde grouped by CRT pattern:")
    print(f"    On-diagonal:  {total_on_crt:.4f} ({total_on_crt/total_crt*100:.1f}%)")
    print(f"    Off-diagonal: {total_off_crt:.4f} ({off_frac_crt*100:.1f}%)")

    if off_frac_crt < 0.01:
        print(f"  → H_tilde is BLOCK-DIAGONAL by CRT pattern (bands respect Brillouin zones)")
    elif off_frac_crt < 0.1:
        print(f"  → H_tilde has WEAK CRT-block structure ({off_frac_crt*100:.1f}% leakage)")
    else:
        print(f"  → H_tilde has NO CRT-block structure (frustrated)")

    # ============================================================
    # 6. IRREP DECOMPOSITION
    # ============================================================
    print(f"\n  ======= 6. IRREP DECOMPOSITION =======")

    projectors, irrep_basis = build_irrep_projectors(orbits, residues, m, phi)

    # How much of H lives in each irrep sector?
    print(f"  H norm² by Klein-4 irrep sector (P_chi H P_chi):")
    irrep_norms = {}
    for label in ['++', '+-', '-+', '--']:
        P_chi = projectors[label]
        H_chi = P_chi @ H @ P_chi
        norm_sq = np.sum(np.abs(H_chi) ** 2)
        frac = norm_sq / H_norm_sq
        irrep_norms[label] = norm_sq
        sigma_label = "σ-even" if label[0] == '+' else "σ-odd"
        tau_label = "τ-even" if label[1] == '+' else "τ-odd"
        print(f"    ({label}) [{sigma_label}, {tau_label}]: "
              f"{norm_sq:.4f} ({frac*100:.1f}%)")

    total_irrep = sum(irrep_norms.values())
    print(f"  Total from 4 diagonal irrep blocks: {total_irrep:.4f} "
          f"({total_irrep/H_norm_sq*100:.1f}%)")
    print(f"  Cross-irrep (off-diagonal between sectors): "
          f"{H_norm_sq - total_irrep:.4f} ({(1-total_irrep/H_norm_sq)*100:.1f}%)")

    # Which irreps carry the active transport?
    # The commutator [|A|, P_tau] should flip tau-parity
    A = (build_directed_D_helper(residues, m) - build_directed_D_helper(residues, m).T) / 2.0
    A_abs = np.abs(A)
    comm = A_abs @ P_tau - P_tau @ A_abs
    print(f"\n  Commutator [|A|_ew, P_tau] by irrep sector:")
    for label in ['++', '+-', '-+', '--']:
        P_chi = projectors[label]
        comm_chi = P_chi @ comm @ P_chi
        norm_sq = np.sum(np.abs(comm_chi) ** 2)
        comm_norm_sq = np.sum(np.abs(comm) ** 2)
        print(f"    ({label}): {norm_sq:.4f} ({norm_sq/comm_norm_sq*100:.1f}%)")

    # Cross-sector: (++) → (+-) and (-+) → (--)  [tau-flip]
    for src, dst in [('++', '+-'), ('-+', '--'), ('++', '-+'), ('+-', '--')]:
        P_src = projectors[src]
        P_dst = projectors[dst]
        cross = P_dst @ comm @ P_src
        norm_sq = np.sum(np.abs(cross) ** 2)
        comm_norm_sq = np.sum(np.abs(comm) ** 2)
        if norm_sq / comm_norm_sq > 0.01:
            print(f"    ({src})→({dst}): {norm_sq:.4f} ({norm_sq/comm_norm_sq*100:.1f}%) "
                  f"{'← τ-flip' if src[0]==dst[0] and src[1]!=dst[1] else ''}"
                  f"{'← σ-flip' if src[0]!=dst[0] and src[1]==dst[1] else ''}"
                  f"{'← both flip' if src[0]!=dst[0] and src[1]!=dst[1] else ''}")

    # ============================================================
    # 7. FULL EIGENDECOMPOSITION COMPARISON
    # ============================================================
    print(f"\n  ======= 7. EIGENSTATE COMPARISON =======")

    # Full H eigenstates projected into orbit space
    evals_full_sorted = np.sort(evals_full)
    evals_tilde_sorted = np.sort(evals_tilde)

    # For Perron and anti-Perron: how well does orbit projection work?
    idx_sorted = np.argsort(evals_full)
    for name, idx in [("Perron", idx_sorted[-1]), ("Anti-Perron", idx_sorted[0])]:
        v = evecs_full[:, idx]
        # Project into orbit space
        v_orb = U @ v  # n_orb-dimensional
        v_reconstructed = U.T @ v_orb  # back to phi-dimensional
        overlap = np.abs(np.vdot(v, v_reconstructed)) ** 2
        print(f"  {name} (λ={evals_full[idx]:.6f}):")
        print(f"    Orbit projection overlap: {overlap:.6f} ({overlap*100:.1f}%)")
        print(f"    Orbit complement: {(1-overlap)*100:.1f}%")

    # Top 5 most localized flat-band states
    flat_mask = np.abs(evals_full) < 0.05
    flat_indices = np.where(flat_mask)[0]
    ipr_vals = np.array([np.sum(np.abs(evecs_full[:, j]) ** 4) * phi
                         for j in range(phi)])
    flat_ipr = [(ipr_vals[j], j) for j in flat_indices]
    flat_ipr.sort(reverse=True)

    print(f"\n  Top 5 localized flat-band states → orbit projection overlap:")
    for rank in range(min(5, len(flat_ipr))):
        ipr_val, j = flat_ipr[rank]
        v = evecs_full[:, j]
        v_orb = U @ v
        v_reconstructed = U.T @ v_orb
        overlap = np.abs(np.vdot(v, v_reconstructed)) ** 2
        print(f"    Rank {rank+1} (IPR×φ={ipr_val:.2f}): "
              f"orbit overlap = {overlap:.4f} ({overlap*100:.1f}%)")

    return {
        'phi': phi,
        'n_orbits': n_orbits,
        'n_blocks': n_blocks,
        'capture_frac': capture_frac,
        'o2_exists': o2_exists,
        'orbits_in_single_block': orbits_in_single_block,
        'orbit_block_spans': orbit_block_spans,
        'off_frac_crt': off_frac_crt,
        'perron_gap_tilde': evals_tilde[-1] - evals_tilde[0],
        'perron_gap_full': evals_full[-1] - evals_full[0],
        'irrep_norms': irrep_norms,
    }


def build_directed_D_helper(residues, m):
    r = np.array(residues)
    return ((r[None, :] - r[:, None]) % m).astype(float)


# ==============================================================
if __name__ == '__main__':
    results = {}
    results['odd'] = run_experiment(1155, "ODD: m = 1155")
    results['even'] = run_experiment(2310, "EVEN: m = 2310")

    # ==== SCORECARD ====
    print(f"\n{'='*70}")
    print(f"  SCORECARD: Arithmetic Crystal Band Structure")
    print(f"{'='*70}")

    r_odd = results['odd']
    r_even = results['even']

    print(f"\n  Test 1: Orbit-space norm capture")
    print(f"    Odd:  {r_odd['capture_frac']*100:.1f}%")
    print(f"    Even: {r_even['capture_frac']*100:.1f}%")
    if r_odd['capture_frac'] > 0.9:
        print(f"    → Orbits capture >90% of physics (clean coarse-graining)")
    elif r_odd['capture_frac'] > 0.5:
        print(f"    → Orbits capture majority but significant intra-orbit physics")
    else:
        print(f"    → Orbits capture <50% — intra-orbit fluctuations dominate")

    print(f"\n  Test 2: Perron gap preserved?")
    print(f"    Odd:  orbit gap = {r_odd['perron_gap_tilde']:.6f},"
          f"  full gap = {r_odd['perron_gap_full']:.6f}"
          f"  (ratio = {r_odd['perron_gap_tilde']/r_odd['perron_gap_full']:.4f})")
    print(f"    Even: orbit gap = {r_even['perron_gap_tilde']:.6f},"
          f"  full gap = {r_even['perron_gap_full']:.6f}"
          f"  (ratio = {r_even['perron_gap_tilde']/r_even['perron_gap_full']:.4f})")

    print(f"\n  Test 3: Do orbits fit inside CRT blocks?")
    print(f"    Odd:  {r_odd['orbits_in_single_block']}/{r_odd['n_orbits']} orbits in single block")
    print(f"    Even: {r_even['orbits_in_single_block']}/{r_even['n_orbits']} orbits in single block")

    print(f"\n  Test 4: Does orbit H respect CRT grouping?")
    print(f"    Odd:  {r_odd['off_frac_crt']*100:.1f}% off-diagonal by CRT pattern")
    print(f"    Even: {r_even['off_frac_crt']*100:.1f}% off-diagonal by CRT pattern")

    print(f"\n  Test 5: Irrep sector balance (odd m)")
    for label in ['++', '+-', '-+', '--']:
        print(f"    ({label}): {r_odd['irrep_norms'][label]:.4f}")

    print(f"\n  VERDICT:")
    if (r_odd['capture_frac'] > 0.9 and
        r_odd['orbits_in_single_block'] == r_odd['n_orbits']):
        print(f"    → CRYSTAL: clean tensor product CRT ⊗ orbit")
    elif r_odd['capture_frac'] > 0.9:
        print(f"    → ORBITS are the right basis (but not CRT-aligned)")
    elif r_odd['off_frac_crt'] < 0.1:
        print(f"    → CRT blocks are the right basis (orbits are too coarse)")
    else:
        print(f"    → FRUSTRATED: neither orbits nor CRT blocks cleanly diagonalize H")
        print(f"       Capture={r_odd['capture_frac']*100:.1f}%,"
              f" CRT-leak={r_odd['off_frac_crt']*100:.1f}%,"
              f" orbits-in-block={r_odd['orbits_in_single_block']}/{r_odd['n_orbits']}")
    print()
