"""
Error Protection Test — Quantifying the Four-Layer Arithmetic Error Protection
===============================================================================
Tests how the σ-parity conservation, flat-band structure, and τ-mixing
degrade under random Hermitian noise at increasing amplitudes.

Key question: Is "four layers of error protection" a quantified claim or
aspirational language?

Tests:
  1. σ-parity breaking under random noise: ||[H+εE, P_σ]||_F vs ε
  2. τ-mixing stability: does the cross-irrep share at α_c survive noise?
  3. Flat-band survival: fraction of eigenvalues in |λ|<0.01 vs ε
  4. Perron eigenvalue shift under noise
  5. Collective encoding: single-site perturbation vs random bulk noise

Runtime: < 5 min (m=2310)
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


def build_P_sigma(residues, m):
    n = len(residues)
    res_to_idx = {r: i for i, r in enumerate(residues)}
    P = np.zeros((n, n))
    for r in residues:
        mr = m - r
        if mr in res_to_idx:
            P[res_to_idx[r], res_to_idx[mr]] = 1.0
    return P


def build_H(D_sym, P_tau, lambda_P, alpha):
    """Build H(α) = D_sym/λ_P + iα[|A|_ew, P_τ]/λ_P"""
    r_arr_placeholder = None  # We pass C directly
    return None  # Use build_H_from_components instead


def build_H_from_components(D_sym, C, lambda_P, alpha):
    """H(α) = D_sym/λ_P + iα·C/λ_P where C = [|A|_ew, P_τ]"""
    return D_sym / lambda_P + 1j * alpha * C / lambda_P


def random_hermitian(n, rng):
    """Generate a random Hermitian matrix with unit Frobenius norm."""
    X = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    H = (X + X.conj().T) / 2.0
    H /= np.linalg.norm(H, 'fro')
    return H


def single_site_perturbation(n, site_idx):
    """Perturbation that affects only one diagonal site."""
    E = np.zeros((n, n), dtype=complex)
    E[site_idx, site_idx] = 1.0
    return E


def symmetric_noise(n, rng, P_sigma):
    """Random Hermitian noise that preserves σ-symmetry: P_σ E P_σ = E."""
    E = random_hermitian(n, rng)
    # Symmetrize under σ: E_sym = (E + P_σ E P_σ) / 2
    E_sym = (E + P_sigma @ E @ P_sigma) / 2.0
    E_sym /= np.linalg.norm(E_sym, 'fro')
    return E_sym


def measure_sigma_breaking(H_perturbed, P_sigma):
    """||[H', P_σ]||_F — measures how much σ-parity is broken."""
    comm = H_perturbed @ P_sigma - P_sigma @ H_perturbed
    return np.linalg.norm(comm, 'fro')


def measure_tau_mixing(H_perturbed, P_sigma, P_tau, phi):
    """Measure cross-irrep share (fraction of H in τ-off-diagonal blocks within σ-even)."""
    S_plus = (np.eye(phi) + P_sigma) / 2.0
    T_plus = (np.eye(phi) + P_tau) / 2.0
    T_minus = (np.eye(phi) - P_tau) / 2.0

    # Project into σ-even sector
    H_sigma_even = S_plus @ H_perturbed @ S_plus

    # τ-off-diagonal block within σ-even
    H_cross = T_plus @ H_sigma_even @ T_minus + T_minus @ H_sigma_even @ T_plus
    H_diag = T_plus @ H_sigma_even @ T_plus + T_minus @ H_sigma_even @ T_minus

    cross_norm = np.linalg.norm(H_cross, 'fro')
    total_norm = np.linalg.norm(H_sigma_even, 'fro')

    if total_norm < 1e-15:
        return 0.0
    return (cross_norm / total_norm * 100.0)


def measure_flatband_survival(H_perturbed, threshold=0.01):
    """Fraction of eigenvalues with |λ| < threshold."""
    evals = np.linalg.eigvalsh(H_perturbed)
    n_flat = np.sum(np.abs(evals) < threshold)
    return n_flat, len(evals)


def measure_perron_shift(H_perturbed):
    """Max eigenvalue of the perturbed Hamiltonian."""
    evals = np.linalg.eigvalsh(H_perturbed)
    return np.max(evals)


def measure_perron_mode_eigenvalue(H_perturbed, v_perron_clean):
    """
    Find the eigenvalue of the state with maximum overlap to clean Perron.
    This tracks the actual collective mode, not spurious defect states.
    """
    evals, evecs = np.linalg.eigh(H_perturbed)
    overlaps = np.abs(v_perron_clean.conj() @ evecs)**2
    best_idx = np.argmax(overlaps)
    return evals[best_idx], overlaps[best_idx]


def run_test():
    m = 2310
    print(f"{'='*72}")
    print(f"  ERROR PROTECTION TEST — m = {m}")
    print(f"{'='*72}")

    residues = get_coprime_residues(m)
    phi = len(residues)
    print(f"  φ = {phi}")

    D_sym = build_D_sym(residues, m)
    P_tau = build_P_tau(residues, m)
    P_sigma = build_P_sigma(residues, m)
    lambda_P = np.max(np.linalg.eigvalsh(D_sym))

    # Build commutator C = [|A|_ew, P_τ]
    r_arr = np.array(residues)
    D_dir = ((r_arr[None, :] - r_arr[:, None]) % m).astype(float)
    A = (D_dir - D_dir.T) / 2.0
    A_abs = np.abs(A)
    C = A_abs @ P_tau - P_tau @ A_abs

    alpha_c = sqrt(135 / 88)
    H_clean = build_H_from_components(D_sym, C, lambda_P, alpha_c)

    # Baseline measurements
    print(f"\n  --- BASELINE (ε = 0) ---")
    sigma_break_0 = measure_sigma_breaking(H_clean, P_sigma)
    tau_mix_0 = measure_tau_mixing(H_clean, P_sigma, P_tau, phi)
    n_flat_0, _ = measure_flatband_survival(H_clean)
    perron_0 = measure_perron_shift(H_clean)
    print(f"  ||[H, P_σ]||_F     = {sigma_break_0:.2e}")
    print(f"  Cross-irrep share  = {tau_mix_0:.1f}%")
    print(f"  Flat-band states   = {n_flat_0}/{phi}")
    print(f"  Perron eigenvalue  = {perron_0:.6f}")

    H_norm = np.linalg.norm(H_clean, 'fro')
    print(f"  ||H||_F            = {H_norm:.2f}")

    rng = np.random.default_rng(42)

    # ================================================================
    # TEST 1: Random Hermitian noise (breaks σ-symmetry)
    # ================================================================
    print(f"\n  {'='*60}")
    print(f"  TEST 1: RANDOM HERMITIAN NOISE (σ-breaking)")
    print(f"  {'='*60}")
    print(f"  {'ε':>10} {'ε/||H||':>10} {'||[H′,Pσ]||':>12} {'τ-mix%':>8} {'Flat':>6} {'λ_max':>10}")
    print(f"  {'-'*58}")

    epsilons = [0.0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]
    E_random = random_hermitian(phi, rng)

    for eps in epsilons:
        H_pert = H_clean + eps * E_random
        sigma_b = measure_sigma_breaking(H_pert, P_sigma)
        tau_m = measure_tau_mixing(H_pert, P_sigma, P_tau, phi)
        n_flat, _ = measure_flatband_survival(H_pert)
        lam_max = measure_perron_shift(H_pert)
        rel_eps = eps / H_norm
        print(f"  {eps:10.4f} {rel_eps:10.6f} {sigma_b:12.4e} {tau_m:8.1f} {n_flat:6d} {lam_max:10.6f}")

    # ================================================================
    # TEST 2: σ-symmetric noise (preserves σ, tests τ + flat band)
    # ================================================================
    print(f"\n  {'='*60}")
    print(f"  TEST 2: σ-SYMMETRIC NOISE (preserves σ-parity)")
    print(f"  {'='*60}")
    print(f"  {'ε':>10} {'||[H′,Pσ]||':>12} {'τ-mix%':>8} {'Flat':>6} {'λ_max':>10}")
    print(f"  {'-'*50}")

    E_sym = symmetric_noise(phi, rng, P_sigma)

    for eps in epsilons:
        H_pert = H_clean + eps * E_sym
        sigma_b = measure_sigma_breaking(H_pert, P_sigma)
        tau_m = measure_tau_mixing(H_pert, P_sigma, P_tau, phi)
        n_flat, _ = measure_flatband_survival(H_pert)
        lam_max = measure_perron_shift(H_pert)
        print(f"  {eps:10.4f} {sigma_b:12.4e} {tau_m:8.1f} {n_flat:6d} {lam_max:10.6f}")

    # ================================================================
    # TEST 3: Single-site perturbation (tests collective encoding)
    # ================================================================
    print(f"\n  {'='*60}")
    print(f"  TEST 3: SINGLE-SITE PERTURBATION (collective encoding)")
    print(f"  {'='*60}")
    print(f"  Perturbing site 0 (residue {residues[0]})")
    print(f"  {'ε':>10} {'||[H′,Pσ]||':>12} {'τ-mix%':>8} {'Flat':>6} {'λ_max':>10} {'Δλ_max':>10}")
    print(f"  {'-'*58}")

    E_site = single_site_perturbation(phi, 0)

    for eps in epsilons:
        H_pert = H_clean + eps * E_site
        sigma_b = measure_sigma_breaking(H_pert, P_sigma)
        tau_m = measure_tau_mixing(H_pert, P_sigma, P_tau, phi)
        n_flat, _ = measure_flatband_survival(H_pert)
        lam_max = measure_perron_shift(H_pert)
        delta_lam = lam_max - perron_0
        print(f"  {eps:10.4f} {sigma_b:12.4e} {tau_m:8.1f} {n_flat:6d} {lam_max:10.6f} {delta_lam:+10.6f}")

    # ================================================================
    # TEST 4: Eigenvector σ-parity degradation under noise
    # ================================================================
    print(f"\n  {'='*60}")
    print(f"  TEST 4: EIGENVECTOR σ-PARITY DEGRADATION")
    print(f"  {'='*60}")
    print(f"  How many eigenstates maintain |<v|P_σ|v>| > 0.99?")
    print(f"  {'ε':>10} {'n(|σ|>0.99)':>12} {'n(|σ|>0.9)':>12} {'min|σ|':>10} {'max dev':>10}")
    print(f"  {'-'*58}")

    for eps in [0.0, 0.001, 0.01, 0.1, 0.3, 1.0, 3.0, 10.0]:
        H_pert = H_clean + eps * E_random
        evals, evecs = np.linalg.eigh(H_pert)
        sigma_exp = np.array([np.real(evecs[:, j].conj() @ P_sigma @ evecs[:, j])
                               for j in range(phi)])
        n_99 = np.sum(np.abs(sigma_exp) > 0.99)
        n_90 = np.sum(np.abs(sigma_exp) > 0.90)
        min_abs = np.min(np.abs(sigma_exp))
        max_dev = np.max(1.0 - np.abs(sigma_exp))
        print(f"  {eps:10.4f} {n_99:12d} {n_90:12d} {min_abs:10.6f} {max_dev:10.6f}")

    # ================================================================
    # TEST 5: PERRON STATE — amplitude compensation under noise
    # ================================================================
    print(f"\n  {'='*60}")
    print(f"  TEST 5: PERRON STATE STABILITY UNDER NOISE")
    print(f"  {'='*60}")

    # Get clean Perron state
    evals_clean, evecs_clean = np.linalg.eigh(H_clean)
    perron_idx_clean = np.argmax(evals_clean)
    v_perron_clean = evecs_clean[:, perron_idx_clean]

    # Clean Perron stats
    probs_clean = np.abs(v_perron_clean)**2
    print(f"  Clean Perron: max prob = {probs_clean.max():.6e}, min = {probs_clean.min():.6e}")
    print(f"  Uniformity: max/min = {probs_clean.max()/probs_clean.min():.6f}")
    print(f"  σ-parity: {np.real(v_perron_clean.conj() @ P_sigma @ v_perron_clean):.6f}")

    print(f"\n  {'ε':>10} {'max_prob':>12} {'min_prob':>12} {'max/min':>10} {'σ-parity':>10} {'|<v0|v>|²':>10}")
    print(f"  {'-'*66}")

    for eps in [0.0, 0.001, 0.01, 0.1, 0.3, 1.0, 3.0]:
        H_pert = H_clean + eps * E_random
        evals_p, evecs_p = np.linalg.eigh(H_pert)
        perron_idx = np.argmax(evals_p)
        v_perron = evecs_p[:, perron_idx]
        probs = np.abs(v_perron)**2
        overlap = np.abs(v_perron_clean.conj() @ v_perron)**2
        sigma_p = np.real(v_perron.conj() @ P_sigma @ v_perron)
        print(f"  {eps:10.4f} {probs.max():12.6e} {probs.min():12.6e} {probs.max()/probs.min():10.4f} {sigma_p:10.6f} {overlap:10.6f}")

    # ================================================================
    # SUMMARY
    # ================================================================
    print(f"\n  {'='*60}")
    print(f"  SUMMARY: ERROR PROTECTION THRESHOLDS")
    print(f"  {'='*60}")
    print(f"""
  Layer 1 (σ-conservation):
    - Algebraic identity → survives ANY perturbation that preserves r↔m-r
    - Random noise: σ-breaking ∝ ε (linear, as expected)
    - Key threshold: see Test 1 and Test 4 for degradation curve

  Layer 2 (coprime sieve):
    - Structural deletion of Fourier modes — not tested by noise
    - Survives by construction (Hilbert space definition)

  Layer 3 (flat-band averaging):
    - See flat-band count column in Tests 1-3
    - Threshold: ε at which flat-band fraction drops below 90%

  Layer 4 (collective encoding):
    - Single-site perturbation (Test 3) vs bulk noise (Test 1)
    - If single-site Δλ_max << bulk Δλ_max at same ε, encoding works

  ||H||_F = {H_norm:.2f}
  Relative noise level ε/||H|| is the physically meaningful comparison.
""")


if __name__ == "__main__":
    run_test()
    
    # ================================================================
    # TEST 6: PERRON MODE TRACKING (fixes the λ_max vs λ_Perron confusion)
    # ================================================================
    print(f"\n  {'='*72}")
    print(f"  TEST 6: PERRON MODE EIGENVALUE UNDER SINGLE-SITE PERTURBATION")
    print(f"  {'='*72}")
    print(f"  This tracks the COLLECTIVE Perron mode, not spurious defect states.")
    print(f"  First-order theory predicts: Δλ_Perron ≈ ε/φ ≈ ε/480 ≈ 0.002 at ε=1.0")
    print()
    
    m = 2310
    residues = sorted([r for r in range(1, m) if gcd(r, m) == 1])
    phi = len(residues)
    
    D_sym = build_D_sym(residues, m)
    P_tau = build_P_tau(residues, m)
    P_sigma = build_P_sigma(residues, m)
    lambda_P = np.max(np.linalg.eigvalsh(D_sym))
    
    r_arr = np.array(residues)
    D_dir = ((r_arr[None, :] - r_arr[:, None]) % m).astype(float)
    A = (D_dir - D_dir.T) / 2.0
    A_abs = np.abs(A)
    C = A_abs @ P_tau - P_tau @ A_abs
    
    alpha_c = sqrt(135 / 88)
    H_clean = build_H_from_components(D_sym, C, lambda_P, alpha_c)
    
    # Get clean Perron eigenvector
    evals_clean, evecs_clean = np.linalg.eigh(H_clean)
    perron_idx = np.argmax(evals_clean)
    v_perron_clean = evecs_clean[:, perron_idx]
    lambda_perron_clean = evals_clean[perron_idx]
    
    print(f"  Clean Perron eigenvalue: λ₀ = {lambda_perron_clean:.6f}")
    print(f"  Expected shift per unit ε: 1/φ = {1/phi:.6f}")
    print()
    
    E_site = single_site_perturbation(phi, 0)
    
    print(f"  {'ε':>10} {'λ_max':>12} {'λ_Perron':>12} {'Δλ_Perron':>12} {'overlap':>10} {'theory':>10}")
    print(f"  {'-'*70}")
    
    for eps in [0.0, 0.001, 0.01, 0.1, 0.3, 1.0, 3.0, 10.0]:
        H_pert = H_clean + eps * E_site
        lam_max = np.max(np.linalg.eigvalsh(H_pert))
        lam_perron, overlap = measure_perron_mode_eigenvalue(H_pert, v_perron_clean)
        delta_lam = lam_perron - lambda_perron_clean
        theory = eps / phi
        print(f"  {eps:10.4f} {lam_max:12.6f} {lam_perron:12.6f} {delta_lam:+12.6f} {overlap:10.6f} {theory:+10.6f}")
    
    print()
    print(f"  CONCLUSION: Single-site Δλ_Perron ≈ ε/φ, NOT the λ_max values!")
    print(f"  The 0.046 in the paper was measuring defect state emergence, not Perron shift.")