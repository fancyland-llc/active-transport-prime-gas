#!/usr/bin/env python3
"""
quantum_prime_dynamics.py — Quantum Time Evolution on the Prime Lattice

Deep Think Challenge (2026-03-25): If H(α) is a tight-binding Hamiltonian on the
coprime residue lattice, then a localized wavepacket |ψ(0)⟩ = |1⟩ should exhibit
dramatically different dynamics at three values of α:

  Test A (α = 0):           Thermal diffusion     → L(t) decays to 1/φ
  Test B (α = 2.5):         Chaotic scrambling     → L(t) crashes instantly
  Test C (α = α_c ≈ 1.239): Flat-band caging      → L(t) remains macroscopically high

L(t) = |⟨1|ψ(t)⟩|² is the Loschmidt echo (survival probability).
Time evolution: |ψ(t)⟩ = exp(-i H(α) t) |ψ(0)⟩ via exact diagonalization.

Usage:
  python quantum_prime_dynamics.py              # Full: m=2310 (φ=480), all three tests
  python quantum_prime_dynamics.py --quick      # Quick: m=210 (φ=48)
  python quantum_prime_dynamics.py --full       # Full: m=2310 + m=30030 comparison
"""

import argparse
import time
import sys
from math import gcd, sqrt

import numpy as np

# ── Matrix construction (from chiral_quench_test.py) ──

def coprime_residues(m):
    """Sorted coprime residues r with 1 <= r < m, gcd(r,m)=1."""
    return sorted(r for r in range(1, m) if gcd(r, m) == 1)

def build_D_sym(residues, m):
    """Palindromic distance matrix: D[i,j] = min(|ri-rj|, m-|ri-rj|)."""
    r = np.array(residues, dtype=np.float64)
    diff = np.abs(r[:, None] - r[None, :])
    return np.minimum(diff, m - diff)

def inversion_perm(residues, m):
    """Permutation sigma where residues[sigma[i]] = residues[i]^{-1} mod m."""
    idx = {r: i for i, r in enumerate(residues)}
    return np.array([idx[pow(r, -1, m)] for r in residues])

def commutator_fast(D, perm):
    """[D, P_tau] = D*P_tau - P_tau*D via column/row permutation."""
    return D[:, perm] - D[perm, :]

def perron_eigenvalue(D, tol=1e-14, max_iter=500):
    """Perron eigenvalue of non-negative symmetric D via power iteration."""
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

def build_H(residues, m, alpha):
    """Build the Active Transport Operator H(α) and return (H, eigenvalues, eigenvectors)."""
    D = build_D_sym(residues, m)
    lam_P = perron_eigenvalue(D)
    perm = inversion_perm(residues, m)
    C = commutator_fast(D, perm)

    A = D / lam_P
    B = C / lam_P
    H = A + (1j * alpha) * B

    evals, evecs = np.linalg.eigh(H)
    return H, evals, evecs, lam_P

# ── Quantum dynamics ──

def loschmidt_echo(evals, evecs, site_idx, t_array):
    """
    Compute L(t) = |⟨site|ψ(t)⟩|² for |ψ(0)⟩ = |site⟩.

    Using spectral decomposition:
      ⟨site|ψ(t)⟩ = Σ_k |⟨site|k⟩|² exp(-i λ_k t)
                   = Σ_k c_k exp(-i λ_k t)
    where c_k = evecs[site_idx, k] and the sum uses the expansion coefficients.

    More precisely:
      |ψ(t)⟩ = Σ_k ⟨k|site⟩ exp(-i λ_k t) |k⟩
      ⟨site|ψ(t)⟩ = Σ_k |⟨site|k⟩|² exp(-i λ_k t)

    Wait — let me be careful. If evecs columns are eigenvectors:
      |ψ(0)⟩ = |site⟩ = Σ_k ⟨k|site⟩ |k⟩ = Σ_k evecs[site,k]* |k⟩
      |ψ(t)⟩ = Σ_k evecs[site,k]* exp(-i λ_k t) |k⟩
      ⟨site|ψ(t)⟩ = Σ_k evecs[site,k] · evecs[site,k]* exp(-i λ_k t)
                   = Σ_k |evecs[site,k]|² exp(-i λ_k t)
    """
    # coefficients c_k = |⟨site|k⟩|² = |evecs[site_idx, k]|²
    c = np.abs(evecs[site_idx, :]) ** 2  # shape (n,)

    # phases exp(-i λ_k t) for each (t, k)
    phases = np.exp(-1j * np.outer(t_array, evals))  # shape (T, n)

    # amplitude ⟨site|ψ(t)⟩ = Σ_k c_k exp(-i λ_k t)
    amplitude = phases @ c  # shape (T,)

    return np.abs(amplitude) ** 2

def participation_ratio(evals, evecs, site_idx):
    """
    Effective number of eigenstates participating in the site state.
    PR = 1 / Σ_k |⟨site|k⟩|⁴
    High PR → state spreads over many eigenstates → fast thermalization.
    Low PR → state overlaps few eigenstates → caging.
    """
    c = np.abs(evecs[site_idx, :]) ** 2
    return 1.0 / np.sum(c ** 2)

def flat_band_fraction(evals, threshold=0.01):
    """Fraction of eigenvalues with |λ| < threshold."""
    return np.sum(np.abs(evals) < threshold) / len(evals)

def time_averaged_echo(L):
    """Running time average of L(t)."""
    return np.cumsum(L) / np.arange(1, len(L) + 1)

# ── Main experiment ──

def run_experiment(m, alpha_values, t_max=200.0, n_t=2000, label=""):
    """Run the quantum dynamics experiment at multiple α values."""
    residues = coprime_residues(m)
    phi = len(residues)
    noise_floor = 1.0 / phi

    # Find site index for r=1 (identity residue)
    site_idx = residues.index(1)

    t_array = np.linspace(0, t_max, n_t)

    print(f"\n{'='*72}")
    print(f"  QUANTUM PRIME DYNAMICS — m = {m}, φ = {phi}{label}")
    print(f"  Initial state: |ψ(0)⟩ = |r=1⟩ (identity residue, site {site_idx})")
    print(f"  Time range: [0, {t_max}], {n_t} points")
    print(f"  Noise floor: 1/φ = {noise_floor:.6f}")
    print(f"{'='*72}")

    results = {}

    for alpha_name, alpha_val in alpha_values:
        t0 = time.perf_counter()
        print(f"\n  [{alpha_name}] α = {alpha_val:.6f} ...", end=" ", flush=True)

        H, evals, evecs, lam_P = build_H(residues, m, alpha_val)

        # Compute Loschmidt echo
        L = loschmidt_echo(evals, evecs, site_idx, t_array)

        # Diagnostics
        pr = participation_ratio(evals, evecs, site_idx)
        fbf = flat_band_fraction(evals)
        L_avg = time_averaged_echo(L)
        L_inf = L_avg[-1]  # time-averaged echo at t_max (proxy for t→∞)
        L_min = np.min(L)
        L_max_late = np.max(L[n_t // 2:])  # max in second half (recurrences)

        dt = time.perf_counter() - t0
        print(f"done ({dt:.1f}s)")

        print(f"    λ_P = {lam_P:.4f}")
        print(f"    Flat-band fraction (|λ| < 0.01): {fbf*100:.1f}%")
        print(f"    Participation ratio: {pr:.1f} / {phi} ({pr/phi*100:.1f}% of Hilbert space)")
        print(f"    L(0)       = {L[0]:.6f}")
        print(f"    L_min      = {L_min:.6f}  (minimum echo)")
        print(f"    L_max(t>T/2) = {L_max_late:.6f}  (max recurrence in late half)")
        print(f"    ⟨L⟩_∞      = {L_inf:.6f}  (time-averaged → infinite-time survival)")
        print(f"    1/φ        = {noise_floor:.6f}  (ergodic floor)")
        print(f"    ⟨L⟩_∞ / (1/φ) = {L_inf * phi:.2f}×  (enhancement over ergodic)")

        # Verdict
        enhancement = L_inf * phi
        if enhancement > 10:
            verdict = "CAGED — wavepacket trapped by flat-band topology"
        elif enhancement > 3:
            verdict = "PARTIALLY LOCALIZED — significant survival enhancement"
        elif enhancement > 1.5:
            verdict = "WEAKLY ENHANCED — slight deviation from ergodic"
        else:
            verdict = "THERMALIZED — ergodic mixing, no caging"

        print(f"    VERDICT: {verdict}")

        results[alpha_name] = {
            "alpha": alpha_val,
            "L": L,
            "L_avg": L_avg,
            "L_inf": L_inf,
            "pr": pr,
            "fbf": fbf,
            "enhancement": enhancement,
            "evals": evals,
        }

    # ── Spectral decomposition insight ──
    print(f"\n{'─'*72}")
    print(f"  SPECTRAL DECOMPOSITION OF |r=1⟩")
    print(f"{'─'*72}")

    # Show how much of |r=1⟩ overlaps with the Perron mode, anti-Perron mode, and flat band
    for alpha_name, alpha_val in alpha_values:
        H, evals, evecs, lam_P = build_H(residues, m, alpha_val)
        c_sq = np.abs(evecs[site_idx, :]) ** 2

        # Identify Perron (largest eval), anti-Perron (most negative), flat band
        i_perron = np.argmax(evals)
        i_anti = np.argmin(evals)
        flat_mask = np.abs(evals) < 0.01
        flat_weight = np.sum(c_sq[flat_mask])
        perron_weight = c_sq[i_perron]
        anti_weight = c_sq[i_anti]

        print(f"\n  [{alpha_name}] α = {alpha_val:.6f}:")
        print(f"    Perron mode (λ={evals[i_perron]:+.4f}):      weight = {perron_weight:.6f}")
        print(f"    Anti-Perron (λ={evals[i_anti]:+.4f}):      weight = {anti_weight:.6f}")
        print(f"    Flat band ({np.sum(flat_mask)} modes, |λ|<0.01): weight = {flat_weight:.6f}")
        print(f"    Perron + anti-Perron:                   weight = {perron_weight + anti_weight:.6f}")
        print(f"    Sum check:                              total  = {np.sum(c_sq):.10f}")

    # ── Summary table ──
    print(f"\n{'='*72}")
    print(f"  SUMMARY — m = {m}, φ = {phi}")
    print(f"{'='*72}")
    print(f"  {'Test':<25} {'α':>8} {'⟨L⟩∞':>10} {'× ergodic':>10} {'PR':>8} {'Flat%':>6} {'Verdict'}")
    print(f"  {'─'*25} {'─'*8} {'─'*10} {'─'*10} {'─'*8} {'─'*6} {'─'*40}")
    for alpha_name, _ in alpha_values:
        r = results[alpha_name]
        enhancement = r["enhancement"]
        if enhancement > 10:
            v = "CAGED"
        elif enhancement > 3:
            v = "PARTIAL"
        elif enhancement > 1.5:
            v = "WEAK"
        else:
            v = "THERMALIZED"
        print(f"  {alpha_name:<25} {r['alpha']:8.4f} {r['L_inf']:10.6f} {r['enhancement']:9.1f}× {r['pr']:7.1f} {r['fbf']*100:5.1f}% {v}")

    print(f"\n  Ergodic floor: 1/φ = {noise_floor:.6f}")
    print()

    return results, t_array

def main():
    parser = argparse.ArgumentParser(description="Quantum time evolution on the prime lattice")
    parser.add_argument("--quick", action="store_true", help="Quick mode: m=210 only")
    parser.add_argument("--full", action="store_true", help="Full mode: m=2310 + m=30030")
    parser.add_argument("--tmax", type=float, default=200.0, help="Maximum time (default: 200)")
    parser.add_argument("--points", type=int, default=2000, help="Number of time points (default: 2000)")
    args = parser.parse_args()

    alpha_c = sqrt(135 / 88)

    alpha_values = [
        ("A: Thermal (α=0)",       0.0),
        ("B: Chaotic (α=2.5)",     2.5),
        ("C: Magic Angle (α_c)",   alpha_c),
    ]

    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  QUANTUM PRIME DYNAMICS — Aharonov-Bohm Caging on Prime Lattice    ║")
    print("║  Deep Think Challenge (2026-03-25)                                  ║")
    print("║  H(α) = D_sym/λ_P + iα[D_sym, P_τ]/λ_P                            ║")
    print("║  Observable: Loschmidt echo L(t) = |⟨1|exp(-iHt)|1⟩|²              ║")
    print(f"║  α_c = √(135/88) = {alpha_c:.10f}                              ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")

    if args.quick:
        moduli = [210]
    elif args.full:
        moduli = [2310, 30030]
    else:
        moduli = [2310]

    all_results = {}
    for m in moduli:
        label = ""
        if m == 30030:
            label = " [CAUTION: 5760×5760 eigenvector decomp, ~2 min]"
        results, t_array = run_experiment(m, alpha_values, t_max=args.tmax, n_t=args.points, label=label)
        all_results[m] = results

    # ── Final interpretation ──
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  INTERPRETATION                                                     ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")

    for m in moduli:
        r = all_results[m]
        e_thermal = r["A: Thermal (α=0)"]["enhancement"]
        e_chaotic = r["B: Chaotic (α=2.5)"]["enhancement"]
        e_magic = r["C: Magic Angle (α_c)"]["enhancement"]

        print(f"\n  m = {m}:")
        print(f"    Thermal (α=0):   {e_thermal:.1f}× ergodic")
        print(f"    Chaotic (α=2.5): {e_chaotic:.1f}× ergodic")
        print(f"    Magic (α_c):     {e_magic:.1f}× ergodic")

        if e_magic > 10 * max(e_thermal, e_chaotic):
            print(f"    → DARK MIRACLE CONFIRMED: α_c produces {e_magic/max(e_thermal, e_chaotic, 1):.0f}× more caging than thermal/chaotic")
            print(f"    → The Prime Gas freezes at the Magic Angle.")
        elif e_magic > 3 * max(e_thermal, e_chaotic):
            print(f"    → SIGNIFICANT CAGING: α_c enhances survival by {e_magic/max(e_thermal, e_chaotic, 1):.1f}× over thermal/chaotic")
        elif e_magic > max(e_thermal, e_chaotic):
            print(f"    → MILD CAGING: α_c shows {e_magic/max(e_thermal, e_chaotic, 1):.1f}× enhancement — suggestive but not dramatic")
        else:
            print(f"    → NO SPECIAL CAGING at α_c. The flat band does not trap the wavepacket.")
            print(f"    → Deep Think's prediction is FALSIFIED at this modulus.")

    print()

if __name__ == "__main__":
    main()
