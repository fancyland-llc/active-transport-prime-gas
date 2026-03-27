#!/usr/bin/env python3
"""
THE CORRECT FLAT-BAND TEST — Using threshold 0.01, not 0.1
==========================================================

The original test used |λ| < 0.01 for the flat band.
My tests used |λ| < 0.1 — 10× too loose!

Author: Claude + Tony (Fancyland LLC)
Date: 2026-03-25
"""

import numpy as np
from math import gcd, sqrt
import warnings
warnings.filterwarnings('ignore')

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

def main():
    print("=" * 80)
    print("  CORRECT FLAT-BAND TEST — Threshold |λ| < 0.01")
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
    V_13 = build_V_diag(residues, 13)
    
    # ==========================================================================
    # STEP 1: Compare thresholds
    # ==========================================================================
    
    print("  Comparing flat-band thresholds:")
    for thresh in [0.01, 0.05, 0.1]:
        n_flat = np.sum(np.abs(evals_H) < thresh)
        print(f"    |λ| < {thresh}: {n_flat} states ({100*n_flat/phi:.1f}%)")
    print()
    
    # Use CORRECT threshold
    flat_threshold = 0.01
    flat_mask = np.abs(evals_H) < flat_threshold
    n_flat = np.sum(flat_mask)
    
    print(f"  Using threshold {flat_threshold}: {n_flat} states")
    print()
    
    # Projector onto flat band
    U_flat = evecs_H[:, flat_mask]
    P_flat = U_flat @ U_flat.conj().T
    
    # ==========================================================================
    # STEP 2: Construct psi_Rabi in CORRECT flat band
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 2: V_13 eigenstates within TIGHT flat band")
    print("=" * 80)
    print()
    
    V_13_FB = U_flat.conj().T @ V_13 @ U_flat
    evals_13, evecs_13 = np.linalg.eigh(V_13_FB)
    
    print(f"  V_13^FB eigenvalue range: [{evals_13.min():.4f}, {evals_13.max():.4f}]")
    print(f"  V_13^FB gap: {evals_13.max() - evals_13.min():.4f}")
    print()
    
    # Max and min eigenstates
    idx_max = np.argmax(evals_13)
    idx_min = np.argmin(evals_13)
    
    psi_max = U_flat @ evecs_13[:, idx_max]
    psi_min = U_flat @ evecs_13[:, idx_min]
    
    psi_max = psi_max / np.linalg.norm(psi_max)
    psi_min = psi_min / np.linalg.norm(psi_min)
    
    # psi_Rabi
    psi_rabi = (psi_max + psi_min) / sqrt(2)
    psi_rabi = psi_rabi / np.linalg.norm(psi_rabi)
    
    # Verify flat-band content
    fb_content = np.abs(np.vdot(psi_rabi, P_flat @ psi_rabi))
    print(f"  psi_Rabi flat-band content: {fb_content:.6f}")
    print()
    
    # ==========================================================================
    # STEP 3: Test Rabi oscillation
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 3: Rabi oscillation under H + V_13")
    print("=" * 80)
    print()
    
    t_max = 100.0
    n_points = 2000
    times = np.linspace(0, t_max, n_points)
    gate_amplitude = 0.5
    
    # Start from psi_rabi, measure against itself
    H_gate = H + gate_amplitude * V_13
    evals_g, evecs_g = np.linalg.eigh(H_gate)
    psi_g = evecs_g.conj().T @ psi_rabi
    
    F_rabi = np.zeros(n_points)
    
    for k, t in enumerate(times):
        phases = np.exp(-1j * evals_g * t)
        psi_t = evecs_g @ (phases * psi_g)
        F_rabi[k] = np.abs(np.vdot(psi_rabi, psi_t))**2
    
    swing = F_rabi.max() - F_rabi.min()
    print(f"  F_rabi range: [{F_rabi.min():.4f}, {F_rabi.max():.4f}]")
    print(f"  Swing: {swing:.4f}")
    
    # Find period
    peaks = []
    for i in range(1, len(F_rabi)-1):
        if F_rabi[i] > F_rabi[i-1] and F_rabi[i] > F_rabi[i+1] and F_rabi[i] > 0.5:
            peaks.append(times[i])
    
    if len(peaks) >= 2:
        period = np.mean(np.diff(peaks))
        print(f"  Period: {period:.2f}")
    
    print()
    
    # ==========================================================================
    # STEP 4: Noise immunity with TIGHT flat band
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 4: Noise immunity (V_5 + V_7)")
    print("=" * 80)
    print()
    
    noise_amp = 0.1
    H_noisy = H + gate_amplitude * V_13 + noise_amp * (V_5 + V_7)
    evals_n, evecs_n = np.linalg.eigh(H_noisy)
    psi_n = evecs_n.conj().T @ psi_rabi
    
    F_rabi_noisy = np.zeros(n_points)
    
    for k, t in enumerate(times):
        phases = np.exp(-1j * evals_n * t)
        psi_t = evecs_n @ (phases * psi_n)
        F_rabi_noisy[k] = np.abs(np.vdot(psi_rabi, psi_t))**2
    
    max_diff = np.max(np.abs(F_rabi - F_rabi_noisy))
    avg_diff = np.mean(np.abs(F_rabi - F_rabi_noisy))
    
    print(f"  Clean gate swing: {swing:.4f}")
    print(f"  Noisy gate swing: {F_rabi_noisy.max() - F_rabi_noisy.min():.4f}")
    print(f"  Max |F_clean - F_noisy|: {max_diff:.4f}")
    print(f"  Avg |F_clean - F_noisy|: {avg_diff:.4f}")
    print()
    
    gate_quality = 1 - max_diff
    print(f"  Gate quality: {gate_quality:.4f}")
    print()
    
    # ==========================================================================
    # STEP 5: Test individual IN-primes
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 5: Individual IN-prime screening")
    print("=" * 80)
    print()
    
    for p, V_p in [(5, V_5), (7, V_7)]:
        H_p = H + gate_amplitude * V_p
        evals_p, evecs_p = np.linalg.eigh(H_p)
        psi_p = evecs_p.conj().T @ psi_rabi
        
        F_p = np.zeros(n_points)
        for k, t in enumerate(times):
            psi_t = evecs_p @ (np.exp(-1j * evals_p * t) * psi_p)
            F_p[k] = np.abs(np.vdot(psi_rabi, psi_t))**2
        
        swing_p = F_p.max() - F_p.min()
        status = "✓ SCREENED" if swing_p < 0.1 else "✗ NOT SCREENED"
        print(f"  V_{p}: swing = {swing_p:.4f}  {status}")
    
    print()
    
    # ==========================================================================
    # VERDICT
    # ==========================================================================
    
    print("=" * 80)
    print("  VERDICT")
    print("=" * 80)
    print()
    
    if swing > 0.9 and max_diff < 0.05:
        print("  ╔════════════════════════════════════════════════════════════════════╗")
        print("  ║  ✓ FLAT-BAND QUBIT CONFIRMED (with correct 0.01 threshold)         ║")
        print("  ║  THIS is what the original test found!                             ║")
        print("  ╚════════════════════════════════════════════════════════════════════╝")
    elif swing > 0.5 and max_diff < 0.1:
        print("  ╔════════════════════════════════════════════════════════════════════╗")
        print("  ║  ~ PARTIAL — Close but not matching original claim                 ║")
        print("  ╚════════════════════════════════════════════════════════════════════╝")
    else:
        print("  ╔════════════════════════════════════════════════════════════════════╗")
        print("  ║  ✗ STILL NOT WORKING — Original result may have been a bug         ║")
        print("  ╚════════════════════════════════════════════════════════════════════╝")
    
    print()
    print(f"  Comparison with transcript claim:")
    print(f"    Transcript: swing = 1.0, max_diff = 0.019")
    print(f"    This test:  swing = {swing:.4f}, max_diff = {max_diff:.4f}")
    print()
    print("=" * 80)

if __name__ == "__main__":
    main()