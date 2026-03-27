#!/usr/bin/env python3
"""
perturbation_robustness.py — Is the flat-band caging topological or fine-tuned?

THE GATE TEST: Add random Hermitian noise to H(α) and measure whether the
Loschmidt echo survival persists.

  - If caging survives noise → structurally protected (topological)
  - If caging dies with small noise → fine-tuned, fragile, not useful

Design:
  1. Build H(α) at m=2310 (φ=480)
  2. Generate random Hermitian noise: N = (R + R†)/2, normalized to ||N||_F = 1
  3. H_noisy = H(α) + ε * N
  4. Sweep ε from 0 to significant fraction of ||H||
  5. Average over multiple noise realizations for statistics
  6. Compare α=0 (strongest caging) vs α=α_c vs α=2.5

Usage:
  python perturbation_robustness.py              # Standard: m=2310, 10 realizations
  python perturbation_robustness.py --quick      # Quick: m=210, 5 realizations
  python perturbation_robustness.py --deep       # Deep: m=2310, 30 realizations
"""

import argparse
import time
import sys
from math import gcd, sqrt

import numpy as np

# ── Matrix construction (shared with quantum_prime_dynamics.py) ──

def coprime_residues(m):
    return sorted(r for r in range(1, m) if gcd(r, m) == 1)

def build_D_sym(residues, m):
    r = np.array(residues, dtype=np.float64)
    diff = np.abs(r[:, None] - r[None, :])
    return np.minimum(diff, m - diff)

def inversion_perm(residues, m):
    idx = {r: i for i, r in enumerate(residues)}
    return np.array([idx[pow(r, -1, m)] for r in residues])

def commutator_fast(D, perm):
    return D[:, perm] - D[perm, :]

def perron_eigenvalue(D, tol=1e-14, max_iter=500):
    n = D.shape[0]
    v = np.ones(n, dtype=np.float64) / np.sqrt(n)
    lam = 0.0
    for _ in range(max_iter):
        w = D @ v
        lam_new = float(np.dot(v, w))
        nrm = np.linalg.norm(w)
        if nrm < 1e-15:
            break
        v = w / nrm
        if abs(lam_new - lam) < tol * max(abs(lam_new), 1.0):
            break
        lam = lam_new
    return lam_new

def build_operators(residues, m):
    """Build normalized A = D/λ_P and B = [D,P_τ]/λ_P."""
    D = build_D_sym(residues, m)
    lam_P = perron_eigenvalue(D)
    perm = inversion_perm(residues, m)
    C = commutator_fast(D, perm)
    return D / lam_P, C / lam_P, lam_P

def build_H(A, B, alpha):
    """H(α) = A + iα B."""
    return A + (1j * alpha) * B

# ── Noise generation ──

def random_hermitian(n, rng):
    """Generate a random Hermitian matrix with Frobenius norm 1."""
    R = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    N = (R + R.conj().T) / 2.0
    N /= np.linalg.norm(N, 'fro')
    return N

# ── Loschmidt echo (time-averaged) ──

def time_averaged_loschmidt(evals, evecs, site_idx, t_max=2000.0, n_t=2000):
    """Compute ⟨L⟩_∞ = time-averaged Loschmidt echo."""
    c_sq = np.abs(evecs[site_idx, :]) ** 2
    t_array = np.linspace(0, t_max, n_t)
    phases = np.exp(-1j * np.outer(t_array, evals))
    amplitude = phases @ c_sq
    L = np.abs(amplitude) ** 2
    return np.mean(L), L

def flat_band_fraction(evals, threshold=0.01):
    return np.sum(np.abs(evals) < threshold) / len(evals)

# ── Main experiment ──

def run_robustness_test(m, alpha_values, epsilon_values, n_realizations=10,
                         t_max=2000.0, n_t=2000, seed=42):
    """Sweep noise amplitude ε and measure echo survival."""
    residues = coprime_residues(m)
    phi = len(residues)
    noise_floor = 1.0 / phi
    site_idx = residues.index(1)

    A, B, lam_P = build_operators(residues, m)
    H_norm = np.linalg.norm(A, 'fro')  # characteristic scale

    rng = np.random.default_rng(seed)

    print(f"\n{'='*78}")
    print(f"  PERTURBATION ROBUSTNESS TEST — m = {m}, φ = {phi}")
    print(f"  ||A||_F = {H_norm:.4f} (energy scale of unperturbed kinetic term)")
    print(f"  Noise model: H_noisy = H(α) + ε × N, ||N||_F = 1")
    print(f"  ε range: [{epsilon_values[0]:.4f}, {epsilon_values[-1]:.4f}]")
    print(f"  Realizations per (α, ε): {n_realizations}")
    print(f"  Time evolution: t ∈ [0, {t_max}], {n_t} points")
    print(f"  Noise floor: 1/φ = {noise_floor:.6f}")
    print(f"{'='*78}")

    # Pre-generate all noise matrices for reproducibility
    n = phi
    noise_matrices = [random_hermitian(n, rng) for _ in range(n_realizations)]

    all_results = {}

    for alpha_name, alpha_val in alpha_values:
        print(f"\n  ── {alpha_name} (α = {alpha_val:.6f}) ──")

        H_clean = build_H(A, B, alpha_val)

        # Clean baseline
        evals_clean, evecs_clean = np.linalg.eigh(H_clean)
        L_inf_clean, _ = time_averaged_loschmidt(evals_clean, evecs_clean, site_idx, t_max, n_t)
        fbf_clean = flat_band_fraction(evals_clean)

        print(f"    Clean:  ⟨L⟩∞ = {L_inf_clean:.6f} ({L_inf_clean/noise_floor:.1f}× ergodic), "
              f"flat-band = {fbf_clean*100:.1f}%")

        results_by_eps = []

        for eps in epsilon_values:
            t0 = time.perf_counter()
            L_values = []
            fbf_values = []

            for noise_mat in noise_matrices:
                H_noisy = H_clean + eps * noise_mat
                evals_n, evecs_n = np.linalg.eigh(H_noisy)
                L_inf_n, _ = time_averaged_loschmidt(evals_n, evecs_n, site_idx, t_max, n_t)
                fbf_n = flat_band_fraction(evals_n)
                L_values.append(L_inf_n)
                fbf_values.append(fbf_n)

            L_mean = np.mean(L_values)
            L_std = np.std(L_values)
            fbf_mean = np.mean(fbf_values)
            dt = time.perf_counter() - t0

            # Relative to clean
            survival_frac = L_mean / L_inf_clean if L_inf_clean > 0 else 0
            enhancement = L_mean / noise_floor

            print(f"    ε={eps:8.4f}: ⟨L⟩∞ = {L_mean:.6f} ± {L_std:.6f} "
                  f"({enhancement:6.1f}× ergodic, {survival_frac*100:5.1f}% of clean) "
                  f"flat={fbf_mean*100:4.1f}%  [{dt:.1f}s]")

            results_by_eps.append({
                'eps': eps,
                'L_mean': L_mean,
                'L_std': L_std,
                'enhancement': enhancement,
                'survival_frac': survival_frac,
                'fbf_mean': fbf_mean,
            })

        all_results[alpha_name] = {
            'alpha': alpha_val,
            'clean_L': L_inf_clean,
            'clean_fbf': fbf_clean,
            'by_eps': results_by_eps,
        }

    # ── Summary ──
    print(f"\n{'='*78}")
    print(f"  ROBUSTNESS SUMMARY — m = {m}, φ = {phi}")
    print(f"{'='*78}")

    # For each α, find the critical ε where caging drops below 50% of clean
    for alpha_name, alpha_val in alpha_values:
        r = all_results[alpha_name]
        print(f"\n  {alpha_name} (α = {alpha_val:.6f}):")
        print(f"    Clean ⟨L⟩∞ = {r['clean_L']:.6f} ({r['clean_L']/noise_floor:.1f}× ergodic)")

        eps_50 = None
        eps_10x = None  # where enhancement drops below 10×
        for entry in r['by_eps']:
            if eps_50 is None and entry['survival_frac'] < 0.5:
                eps_50 = entry['eps']
            if eps_10x is None and entry['enhancement'] < 10:
                eps_10x = entry['eps']

        if eps_50 is not None:
            print(f"    50% survival at ε ≈ {eps_50:.4f} (ε/||A|| = {eps_50/H_norm:.4f})")
        else:
            print(f"    50% survival: NOT REACHED (caging is robust across all tested ε)")

        if eps_10x is not None:
            print(f"    10× ergodic at ε ≈ {eps_10x:.4f} (ε/||A|| = {eps_10x/H_norm:.4f})")
        else:
            print(f"    10× ergodic:  NOT REACHED (caging persists at 10× even at max ε)")

    # ── Verdict ──
    print(f"\n{'='*78}")
    print(f"  VERDICT")
    print(f"{'='*78}")

    # Compare the most-caged (α=0) with α_c — is one more robust than the other?
    r0 = all_results[alpha_values[0][0]]   # α=0
    rc = all_results[alpha_values[2][0]]   # α_c

    # Find survival fraction at the median ε
    mid_idx = len(epsilon_values) // 2
    surv_0_mid = r0['by_eps'][mid_idx]['survival_frac']
    surv_c_mid = rc['by_eps'][mid_idx]['survival_frac']
    eps_mid = epsilon_values[mid_idx]

    # At max ε
    surv_0_max = r0['by_eps'][-1]['survival_frac']
    surv_c_max = rc['by_eps'][-1]['survival_frac']
    enh_0_max = r0['by_eps'][-1]['enhancement']
    enh_c_max = rc['by_eps'][-1]['enhancement']

    print(f"\n  At ε = {eps_mid:.4f} (mid-range):")
    print(f"    α = 0:   {surv_0_mid*100:.1f}% of clean surviving")
    print(f"    α = α_c: {surv_c_mid*100:.1f}% of clean surviving")

    print(f"\n  At ε = {epsilon_values[-1]:.4f} (max tested):")
    print(f"    α = 0:   {surv_0_max*100:.1f}% of clean, {enh_0_max:.1f}× ergodic")
    print(f"    α = α_c: {surv_c_max*100:.1f}% of clean, {enh_c_max:.1f}× ergodic")

    # Is it robust?
    if enh_0_max > 10 and enh_c_max > 10:
        print(f"\n  ★ TOPOLOGICALLY ROBUST: Caging survives {epsilon_values[-1]:.2f} noise")
        print(f"    at BOTH α=0 and α=α_c. The flat-band structure is not fine-tuned.")
        print(f"    The geometry of coprime residues provides structural protection.")
    elif enh_0_max > 10 or enh_c_max > 10:
        print(f"\n  ◆ PARTIALLY ROBUST: Some α values survive, others don't.")
        print(f"    Caging is geometry-dependent, not universal.")
    elif enh_0_max > 3 or enh_c_max > 3:
        print(f"\n  ○ WEAKLY ROBUST: Significant caging remains but eroded by noise.")
        print(f"    The flat band is partially protected.")
    else:
        print(f"\n  ✗ FRAGILE: Caging destroyed by noise. The flat band is fine-tuned.")
        print(f"    No topological protection. The structure is accidental.")

    print()
    return all_results

def main():
    parser = argparse.ArgumentParser(description="Perturbation robustness of flat-band caging")
    parser.add_argument("--quick", action="store_true", help="Quick: m=210, 5 realizations")
    parser.add_argument("--deep", action="store_true", help="Deep: m=2310, 30 realizations")
    parser.add_argument("--m30030", action="store_true", help="Run at m=30030 (slow)")
    args = parser.parse_args()

    alpha_c = sqrt(135 / 88)

    alpha_values = [
        ("A: Thermal (α=0)",       0.0),
        ("B: Chaotic (α=2.5)",     2.5),
        ("C: Magic Angle (α_c)",   alpha_c),
    ]

    if args.quick:
        m = 210
        n_real = 5
        # ε values: 0.001 to 1.0 in log-spaced steps
        epsilon_values = np.array([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0])
    elif args.deep:
        m = 2310
        n_real = 30
        epsilon_values = np.array([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0])
    elif args.m30030:
        m = 30030
        n_real = 5
        epsilon_values = np.array([0.001, 0.01, 0.1, 1.0])
    else:
        m = 2310
        n_real = 10
        epsilon_values = np.array([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0])

    print("╔══════════════════════════════════════════════════════════════════════════════╗")
    print("║  PERTURBATION ROBUSTNESS TEST — Is the flat-band caging topological?        ║")
    print("║  H_noisy = H(α) + ε × N,  N = random Hermitian, ||N||_F = 1                ║")
    print("║  Observable: ⟨L⟩∞ = time-averaged Loschmidt echo                            ║")
    print("║  Question: Does caging survive noise?                                        ║")
    print("╚══════════════════════════════════════════════════════════════════════════════╝")

    run_robustness_test(m, alpha_values, epsilon_values, n_realizations=n_real)

if __name__ == "__main__":
    main()
