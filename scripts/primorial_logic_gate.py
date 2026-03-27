#!/usr/bin/env python3
"""
primorial_logic_gate.py — Does V_13 act as a coherent gate or destructive noise?

Deep Think's proposal:
  1. Initialize |ψ_0⟩ = projection of |r=1⟩ into the flat-band subspace of H(α)
  2. Background noise: H_noise = ε·(V_5 + V_7)  [IN-primes, should be screened]
  3. Gate pulse: H_gate = A·V_13              [MISSING prime, should couple]
  4. Track: F(t) = |⟨ψ_0|e^{-iH_total·t}|ψ_0⟩|²

If F(t) shows clean Rabi oscillations under V_13 while ignoring V_5+V_7 noise,
  → V_13 is a coherent gate in the flat-band subspace.
If F(t) decays monotonically,
  → V_13 is just slower decoherence.

This script uses the FULL Active Transport Operator H(α), not just D_sym.
"""

import numpy as np
from math import gcd, sqrt, pi
import sys

# ── Matrix construction (same as all previous scripts) ──

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

def build_H(residues, m, alpha):
    D = build_D_sym(residues, m)
    lam_P = perron_eigenvalue(D)
    perm = inversion_perm(residues, m)
    C = commutator_fast(D, perm)
    A = D / lam_P
    B = C / lam_P
    H = A + (1j * alpha) * B
    return H, lam_P

def build_prime_potential(residues, p):
    """V_p = diag(cos(2π r / p)), normalized to ||V||_F = 1."""
    r = np.array(residues, dtype=np.float64)
    v = np.cos(2 * pi * r / p)
    V = np.diag(v)
    nrm = np.linalg.norm(V, 'fro')
    if nrm > 0:
        V /= nrm
    return V

def build_prime_potential_unnorm(residues, p):
    """V_p = diag(cos(2π r / p)), unnormalized (for checking constancy)."""
    r = np.array(residues, dtype=np.float64)
    v = np.cos(2 * pi * r / p)
    return np.diag(v)

# ── Flat-band projector ──

def flat_band_projector(evals, evecs, threshold=0.01):
    """Return projector onto flat-band subspace (eigenvalues with |λ| < threshold)."""
    flat_mask = np.abs(evals) < threshold
    P_flat = evecs[:, flat_mask]
    return P_flat, flat_mask

def project_into_flat_band(state, P_flat):
    """Project |state⟩ into flat-band subspace and normalize."""
    projected = P_flat @ (P_flat.conj().T @ state)
    nrm = np.linalg.norm(projected)
    if nrm < 1e-15:
        print("WARNING: State has zero overlap with flat band!")
        return projected
    return projected / nrm

# ── Time evolution via spectral decomposition ──

def fidelity_trace(H_total, psi_0, t_max, n_t):
    """
    Compute F(t) = |⟨ψ_0|e^{-iHt}|ψ_0⟩|² via exact diagonalization.
    Returns (times, fidelities).
    """
    evals, evecs = np.linalg.eigh(H_total)
    c = evecs.conj().T @ psi_0  # expansion coefficients
    c_sq = np.abs(c) ** 2

    times = np.linspace(0, t_max, n_t)
    fidelities = np.zeros(n_t)

    for it, t in enumerate(times):
        phases = np.exp(-1j * evals * t)
        amplitude = np.sum(c_sq * phases)  # this gives ⟨ψ_0|e^{-iHt}|ψ_0⟩ when c are real overlaps
        # Actually: ⟨ψ_0|e^{-iHt}|ψ_0⟩ = Σ_k |c_k|² e^{-iλ_k t}
        fidelities[it] = np.abs(amplitude) ** 2

    return times, fidelities

# ── Main experiment ──

def main():
    m = 2310
    alpha_c = sqrt(135 / 88)

    residues = coprime_residues(m)
    phi = len(residues)
    site_idx = residues.index(1)

    print("=" * 80)
    print(f"  PRIMORIAL LOGIC GATE DYNAMICS")
    print(f"  m = {m}, phi = {phi}, alpha = alpha_c = {alpha_c:.6f}")
    print("=" * 80)

    # ── Build H(α_c) ──
    print(f"\n  Building H(alpha_c)...")
    H0, lam_P = build_H(residues, m, alpha_c)
    evals0, evecs0 = np.linalg.eigh(H0)

    # ── Flat-band analysis ──
    P_flat, flat_mask = flat_band_projector(evals0, evecs0)
    n_flat = np.sum(flat_mask)
    print(f"  Flat-band dimension: {n_flat} / {phi} ({n_flat/phi*100:.1f}%)")
    print(f"  Perron eigenvalue (normalized): {evals0[-1]:.6f}")
    print(f"  Anti-Perron eigenvalue: {evals0[0]:.6f}")

    # ── Initialize: project |r=1⟩ into flat band ──
    psi_site = np.zeros(phi, dtype=complex)
    psi_site[site_idx] = 1.0

    psi_0 = project_into_flat_band(psi_site, P_flat)
    overlap_before = np.abs(np.vdot(psi_site, psi_0)) ** 2
    print(f"\n  |r=1> overlap with flat band: {overlap_before:.6f}")
    print(f"  |psi_0> = projection of |r=1> into flat band (normalized)")

    # ── Build perturbations ──
    # Deep Think says: normalize V_5+V_7 together as the "thermal environment"
    V5_raw = build_prime_potential_unnorm(residues, 5)
    V7_raw = build_prime_potential_unnorm(residues, 7)
    V_noise_raw = V5_raw + V7_raw
    V_noise = V_noise_raw / np.linalg.norm(V_noise_raw, 'fro')

    V_13 = build_prime_potential(residues, 13)  # already normalized

    # Quick sanity: check V_2 and V_3 are constant
    V2_raw = build_prime_potential_unnorm(residues, 2)
    V3_raw = build_prime_potential_unnorm(residues, 3)
    v2_diag = np.diag(V2_raw)
    v3_diag = np.diag(V3_raw)
    print(f"\n  Sanity check:")
    print(f"    V_2 diagonal: min={v2_diag.min():.6f}, max={v2_diag.max():.6f} (constant = {np.std(v2_diag) < 1e-10})")
    print(f"    V_3 diagonal: min={v3_diag.min():.6f}, max={v3_diag.max():.6f} (constant = {np.std(v3_diag) < 1e-10})")

    # ── Parameters ──
    eps_noise = 0.1   # noise amplitude (would kill 90% under random noise)
    A_gate = 0.5      # gate amplitude
    t_max = 200.0     # long enough to see recurrence
    n_t = 4000        # time points

    print(f"\n  Noise amplitude (eps):    {eps_noise}")
    print(f"  Gate amplitude (A):       {A_gate}")
    print(f"  Time range:               [0, {t_max}]")
    print(f"  Time points:              {n_t}")

    # ══════════════════════════════════════════════════════════════
    #  TEST A: Flat-band evolved under H0 alone (control)
    # ══════════════════════════════════════════════════════════════
    print(f"\n{'─'*80}")
    print("  TEST A: Pure H(alpha_c), no perturbation (flat-band control)")
    print(f"{'─'*80}")
    t_A, F_A = fidelity_trace(H0, psi_0, t_max, n_t)
    print(f"  F(0)   = {F_A[0]:.6f}")
    print(f"  F_min  = {min(F_A):.6f}")
    print(f"  F_avg  = {np.mean(F_A):.6f}")
    print(f"  F_max  = {max(F_A):.6f}")

    # ══════════════════════════════════════════════════════════════
    #  TEST B: H0 + IN-prime noise (thermal environment)
    # ══════════════════════════════════════════════════════════════
    print(f"\n{'─'*80}")
    print(f"  TEST B: H(alpha_c) + {eps_noise}*(V_5+V_7) [IN-prime thermal noise]")
    print(f"{'─'*80}")
    H_B = H0 + eps_noise * V_noise
    t_B, F_B = fidelity_trace(H_B, psi_0, t_max, n_t)
    print(f"  F(0)   = {F_B[0]:.6f}")
    print(f"  F_min  = {min(F_B):.6f}")
    print(f"  F_avg  = {np.mean(F_B):.6f}")
    print(f"  F_max  = {max(F_B):.6f}")

    # ══════════════════════════════════════════════════════════════
    #  TEST C: H0 + gate only (V_13, no noise)
    # ══════════════════════════════════════════════════════════════
    print(f"\n{'─'*80}")
    print(f"  TEST C: H(alpha_c) + {A_gate}*V_13 [missing-prime gate, no noise]")
    print(f"{'─'*80}")
    H_C = H0 + A_gate * V_13
    t_C, F_C = fidelity_trace(H_C, psi_0, t_max, n_t)
    print(f"  F(0)   = {F_C[0]:.6f}")
    print(f"  F_min  = {min(F_C):.6f}")
    print(f"  F_avg  = {np.mean(F_C):.6f}")
    print(f"  F_max  = {max(F_C):.6f}")

    # ══════════════════════════════════════════════════════════════
    #  TEST D: H0 + noise + gate (the full test)
    # ══════════════════════════════════════════════════════════════
    print(f"\n{'─'*80}")
    print(f"  TEST D: H(alpha_c) + {eps_noise}*(V_5+V_7) + {A_gate}*V_13 [FULL TEST]")
    print(f"{'─'*80}")
    H_D = H0 + eps_noise * V_noise + A_gate * V_13
    t_D, F_D = fidelity_trace(H_D, psi_0, t_max, n_t)
    print(f"  F(0)   = {F_D[0]:.6f}")
    print(f"  F_min  = {min(F_D):.6f}")
    print(f"  F_avg  = {np.mean(F_D):.6f}")
    print(f"  F_max  = {max(F_D):.6f}")

    # ══════════════════════════════════════════════════════════════
    #  TEST E: H0 + random noise at same amplitude (worst case)
    # ══════════════════════════════════════════════════════════════
    print(f"\n{'─'*80}")
    print(f"  TEST E: H(alpha_c) + {A_gate}*Random + {eps_noise}*Random [random chaos]")
    print(f"{'─'*80}")
    rng = np.random.default_rng(42)
    R1 = rng.standard_normal((phi, phi)) + 1j * rng.standard_normal((phi, phi))
    V_rand1 = (R1 + R1.conj().T) / 2.0
    V_rand1 /= np.linalg.norm(V_rand1, 'fro')
    R2 = rng.standard_normal((phi, phi)) + 1j * rng.standard_normal((phi, phi))
    V_rand2 = (R2 + R2.conj().T) / 2.0
    V_rand2 /= np.linalg.norm(V_rand2, 'fro')

    H_E = H0 + A_gate * V_rand1 + eps_noise * V_rand2
    t_E, F_E = fidelity_trace(H_E, psi_0, t_max, n_t)
    print(f"  F(0)   = {F_E[0]:.6f}")
    print(f"  F_min  = {min(F_E):.6f}")
    print(f"  F_avg  = {np.mean(F_E):.6f}")
    print(f"  F_max  = {max(F_E):.6f}")

    # ══════════════════════════════════════════════════════════════
    #  TEST F: Also try at α=0 for comparison
    # ══════════════════════════════════════════════════════════════
    print(f"\n{'─'*80}")
    print(f"  TEST F: H(alpha=0) + {eps_noise}*(V_5+V_7) + {A_gate}*V_13 [α=0 comparison]")
    print(f"{'─'*80}")
    H0_zero, _ = build_H(residues, m, 0.0)
    evals_z, evecs_z = np.linalg.eigh(H0_zero)
    P_flat_z, _ = flat_band_projector(evals_z, evecs_z)
    psi_0z = project_into_flat_band(psi_site, P_flat_z)

    H_F = H0_zero + eps_noise * V_noise + A_gate * V_13
    t_F, F_F = fidelity_trace(H_F, psi_0z, t_max, n_t)
    print(f"  F(0)   = {F_F[0]:.6f}")
    print(f"  F_min  = {min(F_F):.6f}")
    print(f"  F_avg  = {np.mean(F_F):.6f}")
    print(f"  F_max  = {max(F_F):.6f}")

    # ══════════════════════════════════════════════════════════════
    #  ANALYSIS
    # ══════════════════════════════════════════════════════════════
    print(f"\n{'='*80}")
    print("  ANALYSIS")
    print(f"{'='*80}")

    # Check for recurrence in Test D (the main test)
    # Look for peaks after the initial drop
    from_quarter = n_t // 4
    F_D_later = F_D[from_quarter:]
    max_recurrence = max(F_D_later) if len(F_D_later) > 0 else 0
    avg_later = np.mean(F_D_later) if len(F_D_later) > 0 else 0

    # Also check Test C (gate only, no noise)
    F_C_later = F_C[from_quarter:]
    max_recurrence_C = max(F_C_later) if len(F_C_later) > 0 else 0

    # Count oscillation peaks in Test D
    peaks_D = 0
    for i in range(1, len(F_D) - 1):
        if F_D[i] > F_D[i-1] and F_D[i] > F_D[i+1] and F_D[i] > 0.3:
            peaks_D += 1

    peaks_C = 0
    for i in range(1, len(F_C) - 1):
        if F_C[i] > F_C[i-1] and F_C[i] > F_C[i+1] and F_C[i] > 0.3:
            peaks_C += 1

    print(f"\n  {'Test':<45} {'F_avg':>8} {'F_max (late)':>14} {'Peaks>0.3':>10}")
    print(f"  {'─'*45} {'─'*8} {'─'*14} {'─'*10}")
    print(f"  {'A: Pure H(alpha_c) [control]':<45} {np.mean(F_A):8.4f} {max(F_A[from_quarter:]):14.4f} {'—':>10}")
    print(f"  {'B: + IN-prime noise':<45} {np.mean(F_B):8.4f} {max(F_B[from_quarter:]):14.4f} {'—':>10}")
    print(f"  {'C: + V_13 gate only':<45} {np.mean(F_C):8.4f} {max_recurrence_C:14.4f} {peaks_C:>10}")
    print(f"  {'D: + noise + V_13 gate [MAIN TEST]':<45} {np.mean(F_D):8.4f} {max_recurrence:14.4f} {peaks_D:>10}")
    print(f"  {'E: + random noise [chaos control]':<45} {np.mean(F_E):8.4f} {max(F_E[from_quarter:]):14.4f} {'—':>10}")
    print(f"  {'F: alpha=0 + noise + V_13':<45} {np.mean(F_F):8.4f} {max(F_F[from_quarter:]):14.4f} {'—':>10}")

    # ── Verdict ──
    print(f"\n{'='*80}")
    print("  VERDICT")
    print(f"{'='*80}")

    if max_recurrence > 0.8:
        print("  RABI OSCILLATION DETECTED.")
        print(f"  V_13 drives coherent rotation with recurrence to {max_recurrence:.2%}")
        if np.mean(F_B) > 0.8:
            print("  IN-prime noise is screened during gate operation.")
            print("  >>> COHERENT ARITHMETIC GATE in noise-protected flat band.")
        else:
            print("  But IN-prime noise also degrades the state.")
            print("  >>> Partial gate: rotation exists but armor is imperfect.")
    elif max_recurrence > 0.3:
        print("  PARTIAL RECURRENCE DETECTED.")
        print(f"  V_13 drives partial rotation (max return: {max_recurrence:.2%})")
        print("  >>> Intermediate regime: coherent but lossy.")
    else:
        print("  NO RECURRENCE.")
        print(f"  V_13 drives monotonic decay (max return: {max_recurrence:.2%})")
        print("  >>> V_13 is destructive noise, not a coherent gate.")

    # Does noise affect the gate quality?
    if max_recurrence_C > 0.3 and max_recurrence > 0.3:
        ratio = max_recurrence / max_recurrence_C
        print(f"\n  Gate quality with noise / without noise: {ratio:.3f}")
        if ratio > 0.9:
            print("  IN-prime noise has negligible effect on gate operation.")
        elif ratio > 0.5:
            print("  IN-prime noise partially degrades gate operation.")
        else:
            print("  IN-prime noise severely degrades gate operation.")

    # Ergodic comparison
    ergodic = 1.0 / phi
    print(f"\n  Ergodic baseline: 1/phi = {ergodic:.6f}")
    print(f"  Random noise F_avg = {np.mean(F_E):.6f} ({np.mean(F_E)/ergodic:.1f}x ergodic)")

    print(f"\n{'='*80}")
    print("  Done.")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
