#!/usr/bin/env python3
"""
REVISITING THE FLAT-BAND CAVITY — What Actually Worked
=======================================================

Earlier we found that psi_Cavity = (Perron + Anti-Perron)/√2 under V_13
showed perfect Rabi oscillation with 0.998 noise immunity.

But the multi-gate survey showed the Perron/Anti-Perron qubit is NOT screened.

The difference: the earlier test used the FLAT-BAND PROJECTED state.

Let's reconstruct exactly what worked and understand why.

Author: Claude + Tony (Fancyland LLC)
Date: 2026-03-25
"""

import numpy as np
from math import gcd, sqrt
from scipy.linalg import expm
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CORE CONSTRUCTION
# =============================================================================

def get_coprime_residues(m):
    return sorted([r for r in range(1, m) if gcd(r, m) == 1])

def build_D_sym(residues, m):
    n = len(residues)
    D = np.zeros((n, n))
    for i, r in enumerate(residues):
        for j, s in enumerate(residues):
            d_plus = (r - s) % m
            d_minus = (s - r) % m
            D[i, j] = min(d_plus, d_minus)
    return D

def build_P_tau(residues, m):
    n = len(residues)
    res_to_idx = {r: i for i, r in enumerate(residues)}
    P = np.zeros((n, n))
    for i, r in enumerate(residues):
        r_inv = pow(r, -1, m)
        j = res_to_idx[r_inv]
        P[i, j] = 1.0
    return P

def build_H_alpha(D_sym, P_tau, alpha, lambda_P):
    A = D_sym / lambda_P
    C = D_sym @ P_tau - P_tau @ D_sym
    return A + 1j * alpha * C / lambda_P

def build_V_diag(residues, p):
    return np.diag([np.cos(2 * np.pi * r / p) for r in residues])

# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 80)
    print("  REVISITING THE FLAT-BAND CAVITY — What Actually Worked")
    print("=" * 80)
    print()
    
    m = 2310
    alpha_c = sqrt(135/88)
    
    residues = get_coprime_residues(m)
    phi = len(residues)
    
    D_sym = build_D_sym(residues, m)
    P_tau = build_P_tau(residues, m)
    evals_D, _ = np.linalg.eigh(D_sym)
    lambda_P = np.max(evals_D)
    H = build_H_alpha(D_sym, P_tau, alpha_c, lambda_P)
    evals_H, evecs_H = np.linalg.eigh(H)
    
    V_5 = build_V_diag(residues, 5)
    V_7 = build_V_diag(residues, 7)
    V_11 = build_V_diag(residues, 11)
    V_13 = build_V_diag(residues, 13)
    
    # ==========================================================================
    # REPRODUCE THE EARLIER SUCCESSFUL TEST
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 1: Identify the flat band")
    print("=" * 80)
    print()
    
    # Find Perron and Anti-Perron
    idx_perron = np.argmax(evals_H)
    idx_antiperron = np.argmin(evals_H)
    
    lambda_perron = evals_H[idx_perron]
    lambda_antiperron = evals_H[idx_antiperron]
    
    print(f"  Perron eigenvalue: {lambda_perron:.6f}")
    print(f"  Anti-Perron eigenvalue: {lambda_antiperron:.6f}")
    print(f"  Gap: {lambda_perron - lambda_antiperron:.6f}")
    print()
    
    # Flat band = eigenvalues near 0
    flat_threshold = 0.1
    flat_indices = np.where(np.abs(evals_H) < flat_threshold)[0]
    flat_dim = len(flat_indices)
    
    print(f"  Flat band: {flat_dim} states with |λ| < {flat_threshold}")
    print(f"  Flat band range: [{evals_H[flat_indices].min():.4f}, {evals_H[flat_indices].max():.4f}]")
    print()
    
    # Project V_13 into flat band
    P_flat = evecs_H[:, flat_indices]
    V_13_flat = P_flat.conj().T @ V_13 @ P_flat
    evals_13_flat, evecs_13_flat = np.linalg.eigh(V_13_flat)
    
    print(f"  V_13 eigenvalues within flat band:")
    print(f"    Min: {evals_13_flat.min():.4f}")
    print(f"    Max: {evals_13_flat.max():.4f}")
    print(f"    Gap: {evals_13_flat.max() - evals_13_flat.min():.4f}")
    print()
    
    # ==========================================================================
    # STEP 2: Construct the flat-band cavity state
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 2: Flat-band cavity state")
    print("=" * 80)
    print()
    
    # The cavity state is (|λ_max⟩ + |λ_min⟩)/√2 where λ are V_13 eigenvalues IN THE FLAT BAND
    
    idx_max_13 = np.argmax(evals_13_flat)
    idx_min_13 = np.argmin(evals_13_flat)
    
    # These are in flat-band coordinates, transform to full space
    psi_max_flat = evecs_13_flat[:, idx_max_13]
    psi_min_flat = evecs_13_flat[:, idx_min_13]
    
    psi_max = P_flat @ psi_max_flat
    psi_min = P_flat @ psi_min_flat
    
    # Normalize
    psi_max = psi_max / np.linalg.norm(psi_max)
    psi_min = psi_min / np.linalg.norm(psi_min)
    
    # Cavity state
    psi_cavity = (psi_max + psi_min) / sqrt(2)
    psi_cavity = psi_cavity / np.linalg.norm(psi_cavity)
    
    # Check flat-band content
    fb_content = sum(np.abs(np.vdot(evecs_H[:, i], psi_cavity))**2 for i in flat_indices)
    
    print(f"  |0⟩ = max V_13 eigenstate in flat band")
    print(f"  |1⟩ = min V_13 eigenstate in flat band")
    print(f"  |ψ_cavity⟩ = (|0⟩ + |1⟩)/√2")
    print(f"  Flat-band content: {fb_content:.4f}")
    print()
    
    # ==========================================================================
    # STEP 3: Test V_13 rotation on this state
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 3: V_13 Rabi oscillation on flat-band cavity state")
    print("=" * 80)
    print()
    
    t_max = 100.0
    n_points = 2000
    times = np.linspace(0, t_max, n_points)
    gate_amplitude = 0.5
    
    H_gate = H + gate_amplitude * V_13
    evals_g, evecs_g = np.linalg.eigh(H_gate)
    psi_g = evecs_g.conj().T @ psi_cavity
    
    F_cav = np.zeros(n_points)
    F_0 = np.zeros(n_points)
    F_1 = np.zeros(n_points)
    
    for k, t in enumerate(times):
        phases = np.exp(-1j * evals_g * t)
        psi_t = evecs_g @ (phases * psi_g)
        F_cav[k] = np.abs(np.vdot(psi_cavity, psi_t))**2
        F_0[k] = np.abs(np.vdot(psi_max, psi_t))**2
        F_1[k] = np.abs(np.vdot(psi_min, psi_t))**2
    
    swing_cav = F_cav.max() - F_cav.min()
    
    print(f"  Fidelity with cavity state:")
    print(f"    Range: [{F_cav.min():.4f}, {F_cav.max():.4f}]")
    print(f"    Swing: {swing_cav:.4f}")
    print()
    
    # Find period
    peaks = []
    for i in range(1, len(F_cav)-1):
        if F_cav[i] > F_cav[i-1] and F_cav[i] > F_cav[i+1] and F_cav[i] > 0.5:
            peaks.append(times[i])
    
    if len(peaks) >= 2:
        period = np.mean(np.diff(peaks))
        print(f"  Period: {period:.2f}")
    
    print()
    
    # ==========================================================================
    # STEP 4: Noise immunity test
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 4: Noise immunity (V_5 + V_7 during V_13 gate)")
    print("=" * 80)
    print()
    
    noise_amp = 0.1
    H_noisy = H + gate_amplitude * V_13 + noise_amp * (V_5 + V_7)
    evals_n, evecs_n = np.linalg.eigh(H_noisy)
    psi_n = evecs_n.conj().T @ psi_cavity
    
    F_cav_noisy = np.zeros(n_points)
    
    for k, t in enumerate(times):
        phases = np.exp(-1j * evals_n * t)
        psi_t = evecs_n @ (phases * psi_n)
        F_cav_noisy[k] = np.abs(np.vdot(psi_cavity, psi_t))**2
    
    max_diff = np.max(np.abs(F_cav - F_cav_noisy))
    avg_diff = np.mean(np.abs(F_cav - F_cav_noisy))
    
    print(f"  Clean gate swing: {swing_cav:.4f}")
    print(f"  Noisy gate swing: {F_cav_noisy.max() - F_cav_noisy.min():.4f}")
    print(f"  Max |F_clean - F_noisy|: {max_diff:.4f}")
    print(f"  Avg |F_clean - F_noisy|: {avg_diff:.4f}")
    print()
    
    gate_quality = 1 - max_diff
    print(f"  Gate quality: {gate_quality:.4f}")
    print()
    
    # ==========================================================================
    # STEP 5: Test screening from V_5, V_7, V_11 alone (no V_13)
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 5: IN-prime screening (V_5, V_7, V_11 without V_13)")
    print("=" * 80)
    print()
    
    for p, V_p in [(5, V_5), (7, V_7), (11, V_11)]:
        H_p = H + gate_amplitude * V_p
        evals_p, evecs_p = np.linalg.eigh(H_p)
        psi_p = evecs_p.conj().T @ psi_cavity
        
        F_p = np.zeros(n_points)
        for k, t in enumerate(times):
            phases = np.exp(-1j * evals_p * t)
            psi_t = evecs_p @ (phases * psi_p)
            F_p[k] = np.abs(np.vdot(psi_cavity, psi_t))**2
        
        swing_p = F_p.max() - F_p.min()
        status = "✓ SCREENED" if swing_p < 0.1 else "✗ NOT SCREENED"
        print(f"  V_{p:2d}: swing = {swing_p:.4f}  {status}")
    
    print()
    
    # ==========================================================================
    # STEP 6: The critical comparison — Perron/Anti-Perron vs Flat-Band Cavity
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 6: Perron/Anti-Perron vs Flat-Band Cavity")
    print("=" * 80)
    print()
    
    # Perron/Anti-Perron cavity
    psi_perron = evecs_H[:, idx_perron]
    psi_antiperron = evecs_H[:, idx_antiperron]
    psi_edge_cav = (psi_perron + psi_antiperron) / sqrt(2)
    psi_edge_cav = psi_edge_cav / np.linalg.norm(psi_edge_cav)
    
    # Test noise on both
    # Edge cavity
    psi_edge_g = evecs_g.conj().T @ psi_edge_cav
    psi_edge_n = evecs_n.conj().T @ psi_edge_cav
    
    F_edge_clean = np.zeros(n_points)
    F_edge_noisy = np.zeros(n_points)
    
    for k, t in enumerate(times):
        psi_t_c = evecs_g @ (np.exp(-1j * evals_g * t) * psi_edge_g)
        psi_t_n = evecs_n @ (np.exp(-1j * evals_n * t) * psi_edge_n)
        F_edge_clean[k] = np.abs(np.vdot(psi_edge_cav, psi_t_c))**2
        F_edge_noisy[k] = np.abs(np.vdot(psi_edge_cav, psi_t_n))**2
    
    max_diff_edge = np.max(np.abs(F_edge_clean - F_edge_noisy))
    max_diff_flat = max_diff
    
    print(f"  Perron/Anti-Perron (edge modes):")
    print(f"    Flat-band content: {sum(np.abs(np.vdot(evecs_H[:, i], psi_edge_cav))**2 for i in flat_indices):.4f}")
    print(f"    Max |F_clean - F_noisy|: {max_diff_edge:.4f}")
    print()
    
    print(f"  Flat-Band V_13 Cavity:")
    print(f"    Flat-band content: {fb_content:.4f}")
    print(f"    Max |F_clean - F_noisy|: {max_diff_flat:.4f}")
    print()
    
    print(f"  Screening ratio (flat/edge): {max_diff_edge/max_diff_flat:.1f}×")
    print()
    
    # ==========================================================================
    # FINAL VERDICT
    # ==========================================================================
    
    print("=" * 80)
    print("  FINAL VERDICT")
    print("=" * 80)
    print()
    
    if max_diff_flat < 0.05 and swing_cav > 0.9:
        print("  ╔════════════════════════════════════════════════════════════════════╗")
        print("  ║  ✓ THE CONNECTION IS THE FLAT BAND                                 ║")
        print("  ║                                                                    ║")
        print("  ║  Edge modes (Perron/Anti-Perron): UNPROTECTED                      ║")
        print("  ║  Flat band states: PROTECTED                                       ║")
        print("  ║                                                                    ║")
        print("  ║  Qubit encoding: V_13 eigenstates WITHIN the flat band             ║")
        print("  ║  Gate: V_13 drives transitions between |λ_max⟩ and |λ_min⟩        ║")
        print("  ║  Protection: IN-primes act trivially on flat band                  ║")
        print("  ╚════════════════════════════════════════════════════════════════════╝")
    elif max_diff_flat < max_diff_edge / 2:
        print("  ╔════════════════════════════════════════════════════════════════════╗")
        print("  ║  ~ PARTIAL CONNECTION                                              ║")
        print("  ║  Flat band provides some protection, but not complete              ║")
        print("  ╚════════════════════════════════════════════════════════════════════╝")
    else:
        print("  ╔════════════════════════════════════════════════════════════════════╗")
        print("  ║  ✗ NO DIFFERENTIAL PROTECTION                                      ║")
        print("  ║  Flat band and edge modes equally affected by noise                ║")
        print("  ╚════════════════════════════════════════════════════════════════════╝")
    
    print()
    print("=" * 80)

if __name__ == "__main__":
    main()