"""
Chirality analysis Part 3: The eigenvalue spectrum of H(alpha) for matched pairs.

If |A|/m is isospectral and P_tau is conjugate, are the H(alpha) eigenvalues
also identical? If so, the torsion gap is PURELY an eigenvector phenomenon.
If not, the spectral difference tells us where the chirality lives.

Also: direct test of the "spatial chirality" hypothesis by measuring
eigenvector structure in the coprime position basis.
"""

import numpy as np
from math import gcd

def get_coprimes(m):
    return sorted([r for r in range(1, m) if gcd(r, m) == 1])

def forward_distance_matrix(coprimes, m):
    n = len(coprimes)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                D[i, j] = m
            else:
                D[i, j] = (coprimes[j] - coprimes[i]) % m
    return D

def palindromic_distance_matrix(coprimes, m):
    n = len(coprimes)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                diff = abs(coprimes[j] - coprimes[i])
                D[i, j] = min(diff, m - diff)
    return D

def multiplicative_inverse_perm(coprimes, m):
    n = len(coprimes)
    idx = {r: i for i, r in enumerate(coprimes)}
    P = np.zeros((n, n))
    for i, r in enumerate(coprimes):
        r_inv = pow(r, -1, m)
        P[i, idx[r_inv]] = 1.0
    return P


alpha_c = np.sqrt(135 / 88)

# ============================================================
# PART A: Are H(alpha) eigenvalues identical for matched pairs?
# ============================================================
print("=" * 70)
print("PART A: H(alpha) eigenvalue comparison for matched pairs")
print("=" * 70)

pairs = [
    (15, 30, 8),
    (105, 210, 48),
    (1155, 2310, 480),
]

for m_odd, m_even, phi in pairs:
    print(f"\n--- phi = {phi} ---")
    
    eigs_all = {}
    for m in [m_odd, m_even]:
        coprimes = get_coprimes(m)
        n = len(coprimes)
        
        D_sym = palindromic_distance_matrix(coprimes, m)
        P_tau = multiplicative_inverse_perm(coprimes, m)
        
        lambda_P = np.max(np.linalg.eigvalsh(D_sym))
        B = D_sym @ P_tau - P_tau @ D_sym
        H = D_sym / lambda_P + 1j * alpha_c * B / lambda_P
        
        eigs = np.sort(np.linalg.eigvalsh(H))
        eigs_all[m] = eigs
    
    eigs_odd = eigs_all[m_odd]
    eigs_even = eigs_all[m_even]
    
    # Compare eigenvalue spectra
    max_diff = np.max(np.abs(eigs_odd - eigs_even))
    mean_diff = np.mean(np.abs(eigs_odd - eigs_even))
    
    print(f"  max |eig_odd - eig_even|:  {max_diff:.8e}")
    print(f"  mean |eig_odd - eig_even|: {mean_diff:.8e}")
    
    # Show top and bottom eigenvalues
    print(f"  Top 5 eigenvalues:")
    print(f"    {'ODD':>14s} {'EVEN':>14s} {'Diff':>14s}")
    for k in range(-1, -6, -1):
        print(f"    {eigs_odd[k]:14.8f} {eigs_even[k]:14.8f} {eigs_odd[k]-eigs_even[k]:14.2e}")
    
    # If eigenvalues differ, where is the largest difference?
    if max_diff > 1e-10:
        diff_array = np.abs(eigs_odd - eigs_even)
        max_idx = np.argmax(diff_array)
        print(f"\n  EIGENVALUES DIFFER!")
        print(f"  Largest difference at index {max_idx}: odd={eigs_odd[max_idx]:.8f}, even={eigs_even[max_idx]:.8f}")
        print(f"  This is in the {'bulk' if max_idx < phi * 0.9 else 'top tail'} of the spectrum")
    else:
        print(f"\n  EIGENVALUES IDENTICAL (machine precision)")


# ============================================================
# PART B: The isospectrality of |A|/m — is it exact or approximate?
# ============================================================
print("\n" + "=" * 70)
print("PART B: Testing exactness of |A|/m isospectrality")
print("=" * 70)

for m_odd, m_even, phi in pairs:
    print(f"\n--- phi = {phi} ---")
    
    for m in [m_odd, m_even]:
        coprimes = get_coprimes(m)
        D_fwd = forward_distance_matrix(coprimes, m)
        A = (D_fwd - D_fwd.T) / 2
        A_abs = np.abs(A) / m
        
        eigs = np.sort(np.linalg.eigvalsh(A_abs))
        
        if m == m_odd:
            eigs_ref = eigs
        else:
            max_diff = np.max(np.abs(eigs - eigs_ref))
            mean_diff = np.mean(np.abs(eigs - eigs_ref))
            print(f"  max |eig(|A_odd|/m_odd) - eig(|A_even|/m_even)|: {max_diff:.4e}")
            print(f"  mean:                                              {mean_diff:.4e}")


# ============================================================  
# PART C: The actual spatial structure — gap sequences
# ============================================================
print("\n" + "=" * 70)
print("PART C: Coprime gap structure (the spatial chirality)")
print("=" * 70)

for m_odd, m_even, phi in pairs[:2]:  # just show small ones
    print(f"\n--- phi = {phi} ---")
    
    for m in [m_odd, m_even]:
        coprimes = get_coprimes(m)
        n = len(coprimes)
        
        gaps = []
        for i in range(n):
            j = (i + 1) % n
            if j == 0:
                gap = m - coprimes[-1] + coprimes[0]
            else:
                gap = coprimes[j] - coprimes[i]
            gaps.append(gap)
        
        is_even = (m % 2 == 0)
        label = "EVEN" if is_even else "ODD"
        
        print(f"  m = {m} ({label}): coprimes = {coprimes[:20]}{'...' if n > 20 else ''}")
        print(f"    Gaps: {gaps[:20]}{'...' if n > 20 else ''}")
        print(f"    Gap values: {sorted(set(gaps))}")
        print(f"    Gap statistics: mean={np.mean(gaps):.3f}, std={np.std(gaps):.3f}, "
              f"CV={np.std(gaps)/np.mean(gaps):.4f}")
        
        # Parity of gaps
        if not is_even:
            even_gaps = sum(1 for g in gaps if g % 2 == 0)
            odd_gaps = sum(1 for g in gaps if g % 2 == 1)
            print(f"    Even gaps: {even_gaps}, Odd gaps: {odd_gaps}")


# ============================================================
# PART D: The |A| matrix structure — are odd/even actually conjugate?
# ============================================================
print("\n" + "=" * 70)
print("PART D: Testing matrix conjugacy (not just isospectrality)")
print("=" * 70)

for m_odd, m_even, phi in pairs[:2]:
    print(f"\n--- phi = {phi} ---")
    
    coprimes_odd = get_coprimes(m_odd)
    coprimes_even = get_coprimes(m_even)
    n = phi
    
    # Build |A|/m for both
    D_fwd_odd = forward_distance_matrix(coprimes_odd, m_odd)
    A_odd = (D_fwd_odd - D_fwd_odd.T) / 2
    Aabs_odd = np.abs(A_odd) / m_odd
    
    D_fwd_even = forward_distance_matrix(coprimes_even, m_even)
    A_even = (D_fwd_even - D_fwd_even.T) / 2
    Aabs_even = np.abs(A_even) / m_even
    
    # Build the bijection permutation matrix
    idx_even = {r: i for i, r in enumerate(coprimes_even)}
    sigma = np.zeros((n, n))
    for i, r in enumerate(coprimes_odd):
        if r % 2 == 1:
            mapped = r
        else:
            mapped = r + m_odd
        sigma[i, idx_even[mapped]] = 1.0
    
    # Test: does sigma^T * Aabs_even * sigma = Aabs_odd?
    transported = sigma @ Aabs_even @ sigma.T
    conjugacy_error = np.max(np.abs(transported - Aabs_odd))
    frob_error = np.linalg.norm(transported - Aabs_odd, 'fro')
    
    print(f"  ||sigma * |A_even|/m_even * sigma^T  - |A_odd|/m_odd||_max = {conjugacy_error:.8e}")
    print(f"  ||sigma * |A_even|/m_even * sigma^T  - |A_odd|/m_odd||_F   = {frob_error:.8e}")
    
    if conjugacy_error > 1e-10:
        print(f"  → NOT conjugate! Isospectral but NOT unitarily equivalent via sigma.")
        
        # Show where they differ
        diff = np.abs(transported - Aabs_odd)
        max_ij = np.unravel_index(np.argmax(diff), diff.shape)
        r_i = coprimes_odd[max_ij[0]]
        r_j = coprimes_odd[max_ij[1]]
        print(f"  Largest discrepancy at ({r_i}, {r_j}): "
              f"transported={transported[max_ij]:.6f}, actual={Aabs_odd[max_ij]:.6f}")
    else:
        print(f"  → CONJUGATE via sigma! Same matrix up to permutation.")


# ============================================================
# PART E: Eigenvector localization structure
# ============================================================
print("\n" + "=" * 70)
print("PART E: Eigenvector localization — where on the coprime ring?")
print("=" * 70)

for m_odd, m_even, phi in pairs:
    print(f"\n--- phi = {phi} ---")
    
    for m in [m_odd, m_even]:
        coprimes = get_coprimes(m)
        n = len(coprimes)
        
        D_sym = palindromic_distance_matrix(coprimes, m)
        P_tau = multiplicative_inverse_perm(coprimes, m)
        
        lambda_P = np.max(np.linalg.eigvalsh(D_sym))
        B = D_sym @ P_tau - P_tau @ D_sym
        H = D_sym / lambda_P + 1j * alpha_c * B / lambda_P
        
        eigs, vecs = np.linalg.eigh(H)
        
        # IPR per eigenvector
        ipr_per_vec = np.sum(np.abs(vecs)**4, axis=0)
        
        # Top-3 most localized modes (highest IPR)
        top_idx = np.argsort(-ipr_per_vec)[:3]
        
        # Most delocalized modes (lowest IPR)
        bot_idx = np.argsort(ipr_per_vec)[:3]
        
        is_even = (m % 2 == 0)
        label = "EVEN" if is_even else "ODD"
        
        # Overall statistics
        ipr = np.mean(ipr_per_vec)
        ipr_ratio = ipr / (1.0 / n)
        ipr_median = np.median(ipr_per_vec)
        
        print(f"  m = {m} ({label}): IPR ratio = {ipr_ratio:.4f}x")
        print(f"    IPR distribution: mean={ipr*n:.4f}/n, median={ipr_median*n:.4f}/n, "
              f"max={np.max(ipr_per_vec)*n:.4f}/n, min={np.min(ipr_per_vec)*n:.4f}/n")
        
        # Kurtosis of the IPR distribution
        ipr_normalized = ipr_per_vec * n  # should be ~1 if fully delocalized
        kurtosis = np.mean((ipr_normalized - np.mean(ipr_normalized))**4) / np.var(ipr_normalized)**2
        skew = np.mean((ipr_normalized - np.mean(ipr_normalized))**3) / np.var(ipr_normalized)**1.5
        print(f"    IPR*n statistics: mean={np.mean(ipr_normalized):.4f}, "
              f"std={np.std(ipr_normalized):.4f}, skew={skew:.4f}, kurtosis={kurtosis:.4f}")


# ============================================================
# PART F: The key experiment — what if we use P_tau_even at odd m?
# ============================================================
print("\n" + "=" * 70)
print("PART F: Counterfactual — swap P_tau between matched pairs")
print("=" * 70)

for m_odd, m_even, phi in pairs:
    print(f"\n--- phi = {phi} ---")
    
    coprimes_odd = get_coprimes(m_odd)
    coprimes_even = get_coprimes(m_even)
    n = phi
    
    # Build the bijection sigma
    idx_even = {r: i for i, r in enumerate(coprimes_even)}
    sigma = np.zeros((n, n))
    for i, r in enumerate(coprimes_odd):
        mapped = r if r % 2 == 1 else r + m_odd
        sigma[i, idx_even[mapped]] = 1.0
    
    # Get P_tau for both moduli
    P_tau_odd = multiplicative_inverse_perm(coprimes_odd, m_odd)
    P_tau_even = multiplicative_inverse_perm(coprimes_even, m_even)
    
    # Transport P_tau_even to the odd basis
    P_tau_even_transported = sigma @ P_tau_even @ sigma.T
    
    # Check: is P_tau_odd = P_tau_even_transported?
    ptau_diff = np.max(np.abs(P_tau_odd - P_tau_even_transported))
    print(f"  ||P_tau_odd - sigma*P_tau_even*sigma^T||_max = {ptau_diff:.8e}")
    
    if ptau_diff < 1e-10:
        print(f"  P_tau is exactly conjugate via sigma")
    else:
        n_diff = int(np.sum(np.abs(P_tau_odd - P_tau_even_transported) > 0.5))
        print(f"  P_tau differs at {n_diff} entries")
    
    # Now: build H(alpha) at odd m but with P_tau_even (transported)
    D_sym_odd = palindromic_distance_matrix(coprimes_odd, m_odd)
    lambda_P = np.max(np.linalg.eigvalsh(D_sym_odd))
    
    # Normal H
    B_normal = D_sym_odd @ P_tau_odd - P_tau_odd @ D_sym_odd
    H_normal = D_sym_odd / lambda_P + 1j * alpha_c * B_normal / lambda_P
    eigs_n, vecs_n = np.linalg.eigh(H_normal)
    ipr_normal = np.mean(np.sum(np.abs(vecs_n)**4, axis=0))
    
    # With swapped P_tau
    B_swapped = D_sym_odd @ P_tau_even_transported - P_tau_even_transported @ D_sym_odd
    H_swapped = D_sym_odd / lambda_P + 1j * alpha_c * B_swapped / lambda_P
    eigs_s, vecs_s = np.linalg.eigh(H_swapped)
    ipr_swapped = np.mean(np.sum(np.abs(vecs_s)**4, axis=0))
    
    # Also: build H at odd m but with D_sym_even (transported)
    D_sym_even = palindromic_distance_matrix(coprimes_even, m_even)
    D_sym_even_transported = sigma @ D_sym_even @ sigma.T
    
    lambda_P_even = np.max(np.linalg.eigvalsh(D_sym_even_transported))
    B_dsym_swapped = D_sym_even_transported @ P_tau_odd - P_tau_odd @ D_sym_even_transported
    H_dsym_swapped = D_sym_even_transported / lambda_P_even + 1j * alpha_c * B_dsym_swapped / lambda_P_even
    eigs_ds, vecs_ds = np.linalg.eigh(H_dsym_swapped)
    ipr_dsym_swapped = np.mean(np.sum(np.abs(vecs_ds)**4, axis=0))
    
    # Reference: H at even m
    D_sym_even_direct = palindromic_distance_matrix(coprimes_even, m_even)
    P_tau_even_direct = multiplicative_inverse_perm(coprimes_even, m_even)
    lambda_P_even_d = np.max(np.linalg.eigvalsh(D_sym_even_direct))
    B_even = D_sym_even_direct @ P_tau_even_direct - P_tau_even_direct @ D_sym_even_direct
    H_even = D_sym_even_direct / lambda_P_even_d + 1j * alpha_c * B_even / lambda_P_even_d
    eigs_e, vecs_e = np.linalg.eigh(H_even)
    ipr_even = np.mean(np.sum(np.abs(vecs_e)**4, axis=0))
    
    print(f"\n  IPR ratios:")
    print(f"    Original odd m:                    {ipr_normal * n:.4f}x")
    print(f"    Odd D_sym + Even P_tau (swapped):  {ipr_swapped * n:.4f}x")
    print(f"    Even D_sym + Odd P_tau (swapped):  {ipr_dsym_swapped * n:.4f}x")
    print(f"    Original even m:                   {ipr_even * n:.4f}x")
    
    print(f"\n  → D_sym determines {'' if abs(ipr_swapped - ipr_normal) > abs(ipr_dsym_swapped - ipr_even) else 'NOT '}the IPR")
    
    # Eigenvalue comparison of normal and swapped
    eig_diff_ptau = np.max(np.abs(np.sort(eigs_n) - np.sort(eigs_s)))
    eig_diff_dsym = np.max(np.abs(np.sort(eigs_n) - np.sort(eigs_ds)))
    print(f"  Eigenvalue change from swapping P_tau: {eig_diff_ptau:.6e}")
    print(f"  Eigenvalue change from swapping D_sym: {eig_diff_dsym:.6e}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
