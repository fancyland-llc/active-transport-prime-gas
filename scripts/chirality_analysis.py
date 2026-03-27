"""
Chirality analysis: connecting the Spectral Isotropy eigenphase formula
to the Active Transport torsion crossover.

Key identity: [D_sym, P_tau] = -[|A|, P_tau]
where A = (D - D^T)/2 is the antisymmetric part of the forward distance matrix
and |A| is the element-wise absolute value.

The active transport operator H(alpha) is built entirely from the chirality
magnitude |A| and its commutator with multiplicative inversion.
"""

import numpy as np
from math import gcd
from itertools import combinations
import sys

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

def v2_diagonal(coprimes):
    return np.diag([np.cos(np.pi * r) for r in coprimes])

def compute_ipr(eigvecs):
    """Mean IPR over all eigenvectors."""
    return np.mean(np.sum(np.abs(eigvecs)**4, axis=0))


# ============================================================
# PART 1: Verify the key identity [D_sym, P_tau] = -[|A|, P_tau]
# ============================================================
print("=" * 70)
print("PART 1: Verify [D_sym, P_tau] = -[|A|, P_tau]")
print("=" * 70)

for m in [30, 105, 210, 1155, 2310]:
    coprimes = get_coprimes(m)
    phi = len(coprimes)
    
    D_fwd = forward_distance_matrix(coprimes, m)
    D_sym = palindromic_distance_matrix(coprimes, m)
    P_tau = multiplicative_inverse_perm(coprimes, m)
    
    # Antisymmetric part
    A = (D_fwd - D_fwd.T) / 2
    A_abs = np.abs(A)
    
    # Verify D_sym = m/2 * (J - I) - |A|
    n = phi
    J = np.ones((n, n))
    I = np.eye(n)
    D_sym_predicted = (m / 2) * (J - I) - A_abs
    identity_error = np.max(np.abs(D_sym - D_sym_predicted))
    
    # Verify commutator identity
    comm_Dsym = D_sym @ P_tau - P_tau @ D_sym
    comm_Aabs = A_abs @ P_tau - P_tau @ A_abs
    comm_error = np.max(np.abs(comm_Dsym + comm_Aabs))  # should be 0
    
    # Check A entry type
    A_offdiag = A[~np.eye(n, dtype=bool)]
    is_integer = np.allclose(A_offdiag, np.round(A_offdiag))
    is_half_int = np.allclose(A_offdiag * 2, np.round(A_offdiag * 2)) and not is_integer
    
    print(f"\nm = {m:5d} (phi = {phi:4d}, {'even' if m % 2 == 0 else 'odd'})")
    print(f"  D_sym = m/2*(J-I) - |A| error:    {identity_error:.2e}")
    print(f"  [D_sym, P_tau] + [|A|, P_tau] error: {comm_error:.2e}")
    print(f"  A entries are: {'integers' if is_integer else 'half-integers' if is_half_int else 'mixed'}")
    print(f"  min |A_ij| (off-diag): {np.min(np.abs(A_offdiag)):.2f}")
    print(f"  Has zero entries in |A|: {np.any(np.abs(A_offdiag) < 0.01)}")


# ============================================================
# PART 2: Eigenvalues of A and the cotangent formula
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Eigenvalues of A vs cotangent formula")
print("=" * 70)

for m in [30, 105]:
    coprimes = get_coprimes(m)
    phi = len(coprimes)
    n = phi
    
    D_fwd = forward_distance_matrix(coprimes, m)
    A = (D_fwd - D_fwd.T) / 2
    
    # Eigenvalues of A (should be purely imaginary)
    eigs_A = np.linalg.eigvals(A)
    imag_parts = np.sort(np.imag(eigs_A))[::-1]
    real_parts = np.real(eigs_A)
    
    # Cotangent prediction: eigenvalues of nA should have imag = nm/2 * cot(k*pi/phi)
    # The eigenvalues of A itself: ±(m/2)*cot(k*pi/phi)
    predicted = []
    for k in range(1, phi // 2):
        val = (m / 2) * 1.0 / np.tan(k * np.pi / phi)
        predicted.extend([val, -val])
    predicted.append(0)  # zero eigenvalue if phi is even
    if len(predicted) < phi:
        predicted.append(0)
    predicted = sorted(predicted, reverse=True)
    
    # Compare
    computed_sorted = sorted(imag_parts, reverse=True)
    
    print(f"\nm = {m} (phi = {phi}):")
    print(f"  Max |Re(eig)|: {np.max(np.abs(real_parts)):.2e} (should be ~0)")
    print(f"  Computed    | Predicted (cot) | Error")
    for c, p in zip(computed_sorted[:6], predicted[:6]):
        print(f"  {c:12.4f} | {p:12.4f}      | {abs(c-p):.2e}")


# ============================================================
# PART 3: V2 projection of cotangent modes at odd m
# ============================================================
print("\n" + "=" * 70)
print("PART 3: V2 projection of cotangent modes (even vs odd m)")
print("=" * 70)

# For matched pairs with same phi
pairs = [
    (15, 30, 8),      # phi = 8
    (105, 210, 48),    # phi = 48
    (1155, 2310, 480), # phi = 480
]

for m_odd, m_even, phi in pairs:
    print(f"\n--- Matched pair: m_odd={m_odd}, m_even={m_even}, phi={phi} ---")
    
    for m in [m_odd, m_even]:
        coprimes = get_coprimes(m)
        n = len(coprimes)
        assert n == phi, f"Expected phi={phi}, got {n} for m={m}"
        
        D_fwd = forward_distance_matrix(coprimes, m)
        A = (D_fwd - D_fwd.T) / 2
        A_abs = np.abs(A)
        
        # V2 diagonal
        v2 = np.array([np.cos(np.pi * r) for r in coprimes])
        V2 = np.diag(v2)
        
        # Eigendecompose A (it's antisymmetric, so eigenvalues are purely imaginary)
        # Use A directly: eigenvalues ib_k. Work with iA which is Hermitian.
        iA = 1j * A  # Hermitian: (iA)^dag = -i A^T = -i(-A) = iA
        eigs_iA, vecs_iA = np.linalg.eigh(iA.astype(complex))
        # eigs_iA are real, and equal to -b_k (imaginary parts of A's eigenvalues)
        
        # Sort by magnitude of eigenvalue
        idx = np.argsort(-np.abs(eigs_iA))
        eigs_sorted = eigs_iA[idx]
        vecs_sorted = vecs_iA[:, idx]
        
        # V2 chirality of each mode: <v_k| V2 |v_k>
        v2_chirality = np.array([
            np.real(np.conj(vecs_sorted[:, k]) @ V2 @ vecs_sorted[:, k])
            for k in range(n)
        ])
        
        # Group by chirality magnitude
        even_count = np.sum(np.array(v2) > 0)
        odd_count = np.sum(np.array(v2) < 0)
        
        # Eigenvalues of |A|
        eigs_Aabs = np.sort(np.linalg.eigvalsh(A_abs))[::-1]
        
        is_even = (m % 2 == 0)
        label = "EVEN" if is_even else "ODD"
        
        print(f"\n  m = {m} ({label}):  V2 eigenvalues: +1 count={even_count}, -1 count={odd_count}")
        
        if not is_even:
            # Show V2 chirality for strongest cotangent modes
            print(f"  Top cotangent modes and their V2 chirality <v|V2|v>:")
            print(f"  {'Mode k':>8} {'|eigenval|':>12} {'V2 chirality':>14}")
            for k in range(min(10, n)):
                print(f"  {k+1:>8} {abs(eigs_sorted[k]):>12.4f} {v2_chirality[k]:>14.6f}")
            
            # Overall V2 mixing metric
            # For strong modes (top quartile by |eigenvalue|):
            nq = max(n // 4, 1)
            strong_v2 = np.mean(np.abs(v2_chirality[:nq]))
            weak_v2 = np.mean(np.abs(v2_chirality[-nq:])) if n > nq else 0
            print(f"  Mean |V2 chirality| of strongest {nq} modes: {strong_v2:.6f}")
            print(f"  Mean |V2 chirality| of weakest {nq} modes:   {weak_v2:.6f}")
        else:
            print(f"  V2 = -I (trivial): all V2 chirality = {v2_chirality[0]:.4f}")
        
        print(f"  |A| spectrum: top 5 = {eigs_Aabs[:5].round(3)}")
        print(f"  |A| spectrum: bot 5 = {eigs_Aabs[-5:].round(3)}")


# ============================================================
# PART 4: The IPR decomposition by V2 sectors
# ============================================================
print("\n" + "=" * 70)
print("PART 4: IPR decomposition by V2 parity sectors")
print("=" * 70)

for m_odd, m_even, phi in pairs:
    print(f"\n--- phi = {phi}: m_odd = {m_odd}, m_even = {m_even} ---")
    
    for m in [m_odd, m_even]:
        coprimes = get_coprimes(m)
        n = len(coprimes)
        
        D_sym = palindromic_distance_matrix(coprimes, m)
        P_tau = multiplicative_inverse_perm(coprimes, m)
        
        # Perron eigenvalue
        eigvals_D = np.linalg.eigvalsh(D_sym)
        lambda_P = eigvals_D[-1]
        
        # H(alpha) at alpha_c = sqrt(135/88)
        alpha = np.sqrt(135 / 88)
        B = D_sym @ P_tau - P_tau @ D_sym  # the commutator [D_sym, P_tau]
        H = D_sym / lambda_P + 1j * alpha * B / lambda_P
        
        eigvals_H, eigvecs_H = np.linalg.eigh(H)
        
        # Mean IPR
        ipr = np.mean(np.sum(np.abs(eigvecs_H)**4, axis=0))
        ipr_ratio = ipr / (1.0 / n)
        
        # V2 analysis
        v2 = np.array([np.cos(np.pi * r) for r in coprimes])
        is_even_m = (m % 2 == 0)
        
        if not is_even_m:
            even_mask = v2 > 0  # even coprimes
            odd_mask = v2 < 0   # odd coprimes
            
            # For each H(alpha) eigenvector: how much weight in even vs odd sectors?
            even_weight = np.array([np.sum(np.abs(eigvecs_H[even_mask, k])**2) 
                                     for k in range(n)])
            odd_weight = np.array([np.sum(np.abs(eigvecs_H[odd_mask, k])**2) 
                                    for k in range(n)])
            
            # Parity imbalance per mode
            parity_imbalance = np.abs(even_weight - odd_weight)
            
            print(f"  m = {m} (ODD):  IPR ratio = {ipr_ratio:.4f}x")
            print(f"    Even coprimes: {np.sum(even_mask)}, Odd coprimes: {np.sum(odd_mask)}")
            print(f"    Mean parity imbalance of H eigvecs: {np.mean(parity_imbalance):.6f}")
            print(f"    Max parity imbalance: {np.max(parity_imbalance):.6f}")
            
            # Cross-sector IPR: how much does parity mixing reduce IPR?
            # If all eigenvectors were confined to one parity sector,
            # IPR would be ~2/n (since each sector has ~n/2 components)
            # Parity mixing brings it down toward 1/n
            n_even = np.sum(even_mask)
            n_odd = np.sum(odd_mask)
            ipr_if_confined = 1.0 / min(n_even, n_odd) if n_even > 0 and n_odd > 0 else 1.0/n 
            print(f"    IPR if confined to one sector: {ipr_if_confined:.6f}")
            print(f"    Actual IPR: {ipr:.6f}")
            print(f"    Delocalization ratio: {ipr / ipr_if_confined:.4f}")
        else:
            print(f"  m = {m} (EVEN): IPR ratio = {ipr_ratio:.4f}x")
            print(f"    V2 = -I (no parity sectors)")


# ============================================================
# PART 5: The chirality folding spectrum comparison
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Spectral comparison: |A|_even vs |A|_odd")
print("=" * 70)

for m_odd, m_even, phi in pairs:
    print(f"\n--- phi = {phi}: m_odd = {m_odd}, m_even = {m_even} ---")
    
    for m in [m_odd, m_even]:
        coprimes = get_coprimes(m)
        n = len(coprimes)
        
        D_fwd = forward_distance_matrix(coprimes, m)
        A = (D_fwd - D_fwd.T) / 2
        A_abs = np.abs(A)
        
        # Normalize by modulus for comparison
        A_abs_norm = A_abs / m
        
        eigs = np.sort(np.linalg.eigvalsh(A_abs_norm))[::-1]
        
        is_even = (m % 2 == 0)
        label = "EVEN" if is_even else "ODD"
        
        # Spectral statistics
        spectral_gap = eigs[0] - eigs[1] if n > 1 else 0
        trace = np.sum(eigs)  # = trace of |A|/m
        
        print(f"  m = {m} ({label}):")
        print(f"    |A|/m top 3 eigs:  {eigs[:3].round(6)}")
        print(f"    |A|/m spectral gap: {spectral_gap:.6f}")
        print(f"    trace(|A|/m):       {trace:.4f}")
        print(f"    Frobenius ||A||/m:  {np.sqrt(np.sum(A_abs_norm**2)):.4f}")
        
        # Entropy of |A| spectrum (normalized)
        eigs_pos = eigs[eigs > 1e-10]
        eigs_norm = eigs_pos / np.sum(eigs_pos)
        entropy = -np.sum(eigs_norm * np.log(eigs_norm))
        print(f"    Spectral entropy:   {entropy:.4f}")


# ============================================================
# PART 6: Key prediction — R*cot(k*pi/phi) vs |A| eigval overlap
# ============================================================
print("\n" + "=" * 70)
print("PART 6: How cotangent modes project onto |A| eigenvectors")
print("=" * 70)

for m_odd, m_even, phi in [(15, 30, 8), (105, 210, 48)]:
    print(f"\n--- phi = {phi} ---")
    
    for m in [m_odd, m_even]:
        coprimes = get_coprimes(m)
        n = len(coprimes)
        
        D_fwd = forward_distance_matrix(coprimes, m)
        A = (D_fwd - D_fwd.T) / 2
        A_abs = np.abs(A)
        
        # Eigenvectors of iA (Hermitian form of A)
        iA = 1j * A
        eigs_iA, vecs_A = np.linalg.eigh(iA.astype(complex))
        
        # Eigenvectors of |A|
        eigs_Aabs, vecs_Aabs = np.linalg.eigh(A_abs)
        idx_aabs = np.argsort(-eigs_Aabs)
        vecs_Aabs = vecs_Aabs[:, idx_aabs]
        eigs_Aabs = eigs_Aabs[idx_aabs]
        
        # Overlap matrix: how do A modes project onto |A| modes?
        overlap = np.abs(np.conj(vecs_A.T) @ vecs_Aabs)**2  # |<A_mode|Aabs_mode>|^2
        
        # For each |A| mode: how many A modes contribute significantly?
        participation = []
        for j in range(min(8, n)):
            col = overlap[:, j]
            # Effective number of contributing A modes (inverse participation in A basis)
            p = 1.0 / np.sum(col**2) if np.sum(col**2) > 0 else 0
            participation.append(p)
        
        is_even = (m % 2 == 0)
        label = "EVEN" if is_even else "ODD"
        print(f"\n  m = {m} ({label}):")
        print(f"    |A| mode | Effective # of A modes contributing")
        for j, p in enumerate(participation):
            print(f"    {j+1:>7}  | {p:.2f}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
