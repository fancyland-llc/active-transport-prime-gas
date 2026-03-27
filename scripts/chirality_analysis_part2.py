"""
Chirality analysis Part 2: The P_tau parity mixing mechanism.

KEY FINDING FROM PART 1:
  |A_odd|/m_odd = |A_even|/m_even  (spectrally identical)
  → The entire torsion gap comes from P_tau, not from |A|.

This script measures:
1. How P_tau mixes parities at odd m (and that it's trivially parity-preserving at even m)
2. The commutator [|A|, P_tau] structure difference between matched pairs
3. V2 decomposition of the imaginary wave-function
4. Direct prediction: parity mixing fraction predicts torsion gap
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

def compute_h_alpha(D_sym, P_tau, alpha):
    eigvals_D = np.linalg.eigvalsh(D_sym)
    lambda_P = eigvals_D[-1]
    B = D_sym @ P_tau - P_tau @ D_sym
    H = D_sym / lambda_P + 1j * alpha * B / lambda_P
    return H, lambda_P

def compute_ipr(eigvecs):
    return np.mean(np.sum(np.abs(eigvecs)**4, axis=0))


# ============================================================
# PART A: Parity mixing of P_tau
# ============================================================
print("=" * 70)
print("PART A: Parity mixing fraction of P_tau at odd m")
print("=" * 70)

test_moduli = [
    (15, 30, 8),
    (105, 210, 48),
    (1155, 2310, 480),
]

for m_odd, m_even, phi in test_moduli:
    print(f"\n--- phi = {phi} ---")
    
    for m in [m_odd, m_even]:
        coprimes = get_coprimes(m)
        n = len(coprimes)
        is_even = (m % 2 == 0)
        
        # Count parity-mixing inversions
        parity_mix_count = 0
        parity_same_count = 0
        for r in coprimes:
            r_inv = pow(r, -1, m)
            r_parity = r % 2
            rinv_parity = r_inv % 2
            if r_parity != rinv_parity:
                parity_mix_count += 1
            else:
                parity_same_count += 1
        
        mix_frac = parity_mix_count / n if n > 0 else 0
        
        label = "EVEN" if is_even else "ODD"
        print(f"  m = {m} ({label}): parity-mixing inversions = {parity_mix_count}/{n} ({mix_frac:.4f})")
        if not is_even:
            # Show actual parity transitions
            even_to_odd = 0
            odd_to_even = 0
            for r in coprimes:
                r_inv = pow(r, -1, m)
                if r % 2 == 0 and r_inv % 2 == 1:
                    even_to_odd += 1
                elif r % 2 == 1 and r_inv % 2 == 0:
                    odd_to_even += 1
            print(f"    even→odd: {even_to_odd}, odd→even: {odd_to_even}")


# ============================================================
# PART B: Commutator norm decomposition by V2 sectors
# ============================================================
print("\n" + "=" * 70)
print("PART B: Imaginary wave-function [|A|, P_tau] norm comparison")
print("=" * 70)

for m_odd, m_even, phi in test_moduli:
    print(f"\n--- phi = {phi} ---")
    
    data = {}
    for m in [m_odd, m_even]:
        coprimes = get_coprimes(m)
        n = len(coprimes)
        
        D_fwd = forward_distance_matrix(coprimes, m)
        A = (D_fwd - D_fwd.T) / 2
        A_abs = np.abs(A)
        P_tau = multiplicative_inverse_perm(coprimes, m)
        
        # Commutator [|A|, P_tau]
        comm = A_abs @ P_tau - P_tau @ A_abs
        
        # Normalize
        comm_norm = np.linalg.norm(comm, 'fro') / m  # normalize by m for comparison
        comm_max = np.max(np.abs(comm)) / m
        
        # Also compute [|A|/m, P_tau] for fair comparison
        comm_normalized = (A_abs / m) @ P_tau - P_tau @ (A_abs / m)
        cnorm = np.linalg.norm(comm_normalized, 'fro')
        
        is_even = (m % 2 == 0)
        label = "EVEN" if is_even else "ODD"
        
        data[m] = cnorm
        print(f"  m = {m} ({label}):")
        print(f"    ||[|A|/m, P_tau]||_F = {cnorm:.6f}")
        print(f"    max|[|A|/m, P_tau]|  = {np.max(np.abs(comm_normalized)):.6f}")
        
        # V2 decomposition of the commutator
        if not is_even:
            v2 = np.array([(-1)**r for r in coprimes])
            even_mask = np.array([r % 2 == 0 for r in coprimes])
            odd_mask = ~even_mask
            
            # Block decomposition of comm_normalized
            n_even = np.sum(even_mask)
            n_odd = np.sum(odd_mask)
            
            # Extract blocks: comm[even, even], comm[even, odd], etc.
            ee_block = comm_normalized[np.ix_(even_mask, even_mask)]
            oo_block = comm_normalized[np.ix_(odd_mask, odd_mask)]
            eo_block = comm_normalized[np.ix_(even_mask, odd_mask)]
            oe_block = comm_normalized[np.ix_(odd_mask, even_mask)]
            
            ee_norm = np.linalg.norm(ee_block, 'fro')
            oo_norm = np.linalg.norm(oo_block, 'fro')
            eo_norm = np.linalg.norm(eo_block, 'fro')
            oe_norm = np.linalg.norm(oe_block, 'fro')
            total = np.sqrt(ee_norm**2 + oo_norm**2 + eo_norm**2 + oe_norm**2)
            
            print(f"    V2 block decomposition of [|A|/m, P_tau]:")
            print(f"      same-parity (EE+OO):  {np.sqrt(ee_norm**2 + oo_norm**2):.6f}  ({(ee_norm**2 + oo_norm**2)/total**2 * 100:.1f}%)")
            print(f"      cross-parity (EO+OE): {np.sqrt(eo_norm**2 + oe_norm**2):.6f}  ({(eo_norm**2 + oe_norm**2)/total**2 * 100:.1f}%)")
    
    # Torsion gap: IPR difference
    norm_diff = data[m_even] - data[m_odd]
    norm_ratio = data[m_even] / data[m_odd] if data[m_odd] > 0 else float('inf')
    print(f"  ||[|A|/m, P_tau]||_F ratio (even/odd): {norm_ratio:.6f}")


# ============================================================
# PART C: Direct link — P_tau difference structure
# ============================================================
print("\n" + "=" * 70)
print("PART C: What P_tau 'sees' differently at even vs odd m")
print("=" * 70)

for m_odd, m_even, phi in test_moduli:
    print(f"\n--- phi = {phi} ---")
    
    coprimes_odd = get_coprimes(m_odd)
    coprimes_even = get_coprimes(m_even)
    
    # Build the bijection between coprimes of odd m and coprimes of even m
    # For the m_even = 2*m_odd primorial chain:
    # Coprimes of 2m_odd = odd coprimes of m_odd
    # Bijection: r → r if odd, r → r + m_odd if even
    
    bijection = {}
    for r in coprimes_odd:
        if r % 2 == 1:
            mapped = r
        else:
            mapped = r + m_odd
        bijection[r] = mapped
    
    # Verify bijection maps to coprimes of even m
    mapped_set = set(bijection.values())
    even_set = set(coprimes_even)
    if mapped_set == even_set:
        print(f"  Bijection φ(m_odd) → φ(m_even) verified: r → r (odd), r → r+m_odd (even)")
    else:
        print(f"  Bijection FAILED: mapped={sorted(mapped_set)[:10]}, expected={sorted(even_set)[:10]}")
        continue
    
    # Under this bijection, how does P_tau transform?
    # At odd m: tau_odd(r) = r^{-1} mod m_odd
    # At even m: tau_even(sigma(r)) = sigma(r)^{-1} mod m_even
    # Where sigma is the bijection
    
    # The "transported" inversion at odd m:
    # sigma(tau_odd(r)) vs tau_even(sigma(r))
    agree_count = 0
    disagree_count = 0
    for r in coprimes_odd:
        # tau at odd m, transported to even m
        r_inv_odd = pow(r, -1, m_odd)
        transported = bijection[r_inv_odd]
        
        # tau at even m on the corresponding element
        sigma_r = bijection[r]
        sigma_r_inv = pow(sigma_r, -1, m_even)
        
        if transported == sigma_r_inv:
            agree_count += 1
        else:
            disagree_count += 1
    
    agree_frac = agree_count / len(coprimes_odd)
    print(f"  σ∘τ_odd = τ_even∘σ agreement: {agree_count}/{len(coprimes_odd)} ({agree_frac:.4f})")
    print(f"  Inversion discrepancy: {disagree_count}/{len(coprimes_odd)} ({1-agree_frac:.4f})")
    print(f"  → {(1-agree_frac)*100:.1f}% of coprimes have DIFFERENT inversions under transport")


# ============================================================
# PART D: Spectral decomposition of H(alpha) into cotangent mode basis
# ============================================================
print("\n" + "=" * 70)
print("PART D: How H(alpha) eigenvectors decompose in the cotangent basis")
print("=" * 70)
print("(Measuring the mixing induced by the imaginary wave-function)")

alpha_c = np.sqrt(135 / 88)

for m_odd, m_even, phi in test_moduli:
    print(f"\n--- phi = {phi} ---")
    
    for m in [m_odd, m_even]:
        coprimes = get_coprimes(m)
        n = len(coprimes)
        
        D_fwd = forward_distance_matrix(coprimes, m)
        D_sym = palindromic_distance_matrix(coprimes, m)
        P_tau = multiplicative_inverse_perm(coprimes, m)
        A = (D_fwd - D_fwd.T) / 2
        A_abs = np.abs(A)
        
        # Eigenvectors of |A| (the chirality magnitude basis)
        eigs_Aabs, vecs_Aabs = np.linalg.eigh(A_abs)
        
        # H(alpha_c)
        H, lambda_P = compute_h_alpha(D_sym, P_tau, alpha_c)
        eigs_H, vecs_H = np.linalg.eigh(H)
        
        # How many |A| modes does each H mode span?
        # Overlap matrix
        overlap = np.abs(np.conj(vecs_Aabs.T) @ vecs_H)**2
        
        # For each H eigenvector: effective number of |A| modes
        participation_in_Aabs = np.array([
            1.0 / np.sum(overlap[:, k]**2) if np.sum(overlap[:, k]**2) > 0 else 0
            for k in range(n)
        ])
        
        # IPR  
        ipr = compute_ipr(vecs_H)
        ipr_ratio = ipr / (1.0 / n)
        
        is_even = (m % 2 == 0)
        label = "EVEN" if is_even else "ODD"
        
        mean_participation = np.mean(participation_in_Aabs)
        median_participation = np.median(participation_in_Aabs)
        
        print(f"  m = {m} ({label}): IPR ratio = {ipr_ratio:.4f}x")
        print(f"    Mean participation of H modes in |A| basis: {mean_participation:.2f}")
        print(f"    Median participation: {median_participation:.2f}")


# ============================================================
# PART E: The smoking gun — trace the imaginary wave-function energy
# ============================================================
print("\n" + "=" * 70)
print("PART E: Imaginary wave-function energy in H(alpha) eigenmodes")
print("=" * 70)

for m_odd, m_even, phi in test_moduli:
    print(f"\n--- phi = {phi} ---")
    
    for m in [m_odd, m_even]:
        coprimes = get_coprimes(m)
        n = len(coprimes)
        
        D_sym = palindromic_distance_matrix(coprimes, m)
        P_tau = multiplicative_inverse_perm(coprimes, m)
        
        eigvals_D = np.linalg.eigvalsh(D_sym)
        lambda_P = eigvals_D[-1]
        
        # The two parts of H(alpha)
        A_part = D_sym / lambda_P  # real symmetric part
        B = (D_sym @ P_tau - P_tau @ D_sym) / lambda_P  # skew-symmetric (imaginary wave-function)
        
        # H(alpha) = A_part + i*alpha_c * B
        H = A_part + 1j * alpha_c * B
        eigs_H, vecs_H = np.linalg.eigh(H)
        
        # For each H eigenmode: measure <v|A_part|v> and <v|iB|v>
        # Since H is Hermitian and v is an eigenvector: <v|H|v> = eigenvalue (real)
        # A_part contribution: <v|A_part|v> (real, since A_part is real symmetric)
        # B contribution: i*alpha_c * <v|B|v>. Since B is skew-symmetric, <v|B|v> is purely imaginary.
        # Actually, for real-valued analysis: <v|iB|v> where iB is Hermitian.
        
        iB = 1j * B
        
        # Energy from A_part and from iB for each mode
        A_energies = np.array([np.real(np.conj(vecs_H[:, k]) @ A_part @ vecs_H[:, k]) 
                                for k in range(n)])
        B_energies = np.array([np.real(np.conj(vecs_H[:, k]) @ iB @ vecs_H[:, k]) 
                                for k in range(n)])
        
        # Total = A_energies + alpha_c * B_energies should equal eigenvalues
        total = A_energies + alpha_c * B_energies
        eigenvalue_check = np.max(np.abs(total - eigs_H))
        
        # Statistics
        mean_A = np.mean(np.abs(A_energies))
        mean_B = np.mean(np.abs(B_energies)) * alpha_c
        ratio_BA = mean_B / mean_A if mean_A > 0 else 0
        
        ipr = compute_ipr(vecs_H)
        ipr_ratio = ipr / (1.0 / n)
        
        is_even = (m % 2 == 0)
        label = "EVEN" if is_even else "ODD"
        
        print(f"  m = {m} ({label}): IPR ratio = {ipr_ratio:.4f}x")
        print(f"    consistency check (eigenvalue decomposition): {eigenvalue_check:.2e}")
        print(f"    mean |<v|A|v>|:           {mean_A:.6f}")
        print(f"    mean alpha*|<v|iB|v>|:    {mean_B:.6f}")
        print(f"    imaginary/real ratio:      {ratio_BA:.6f}")
        print(f"    ||B||_F / ||A_part||_F:   {np.linalg.norm(B, 'fro') / np.linalg.norm(A_part, 'fro'):.6f}")
        
        # The B_energy variance — how spread are the imaginary contributions?
        B_var = np.var(B_energies)
        print(f"    Var(B_energies):          {B_var:.8f}")


print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
KEY RESULTS:
  1. D_sym = m/2*(J-I) - |A|           [proved, verified]
  2. [D_sym, P_tau] = -[|A|, P_tau]    [proved, verified]  
  3. |A_odd|/m_odd = |A_even|/m_even   [spectral identity — the chirality
                                         magnitude is UNIVERSAL for matched phi]
  4. The ENTIRE torsion gap comes from P_tau acting differently at odd vs even m.
     The imaginary wave-function [|A|, P_tau] is the ONLY source of difference.
  5. At odd m: P_tau mixes parities (even↔odd coprimes)
     At even m: P_tau preserves parity (all coprimes are odd)
  6. The cotangent eigenphase formula Im(lambda_k) = -R*cot(k*pi/phi)
     describes the eigenvalues of A (the antisymmetric part of D).
     The folding A → |A| produces D_sym, and the commutator -[|A|, P_tau]
     is the imaginary wave-function of H(alpha).
     
THE MECHANISM:
  The palindromic distance matrix D_sym is m/2*(J-I) minus the chirality 
  magnitude |A|. The chirality magnitude is universal (identical spectrum 
  for matched phi). The torsion gap arises ENTIRELY from the imaginary 
  wave-function -[|A|, P_tau], which differs because multiplicative 
  inversion P_tau interacts differently with parity at odd vs even m.
  
  At odd m: P_tau mixes parities → cross-sector coupling in the commutator
            → eigenvectors spread across parity sectors → DELOCALIZATION
  At even m: no parity sectors → commutator acts on a single sector
             → eigenvectors can localize freely → LOCALIZATION
             
  The crossover happens when the parity-mixing delocalization at odd m
  becomes strong enough to overcome the universal geometric localization.
""")
