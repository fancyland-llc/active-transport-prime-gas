#!/usr/bin/env python3
"""
THE STATE THAT WORKED — Reconstructing psi_Rabi
================================================

From the transcript:
"psi_Rabi (flat-band optimal): 1.000 swing, max diff = 0.019"

psi_Rabi = (|λ_max⟩ + |λ_min⟩)/√2 where λ are V_13 eigenvalues
WITHIN THE FLAT BAND.

This is different from the Perron/Anti-Perron cavity (edge modes).

Author: Claude + Tony (Fancyland LLC)
Date: 2026-03-25
"""

import numpy as np
from math import gcd, sqrt
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
    print("  THE STATE THAT WORKED — Reconstructing psi_Rabi")
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
    # STEP 1: Identify flat band
    # ==========================================================================
    
    print("  STEP 1: Flat band identification")
    print()
    
    flat_threshold = 0.1
    flat_indices = np.where(np.abs(evals_H) < flat_threshold)[0]
    flat_dim = len(flat_indices)
    
    print(f"  Flat band: {flat_dim} states with |λ| < {flat_threshold}")
    
    # Projector onto flat band
    P_flat = evecs_H[:, flat_indices]
    
    # ==========================================================================
    # STEP 2: Construct psi_Rabi within flat band
    # ==========================================================================
    
    print()
    print("  STEP 2: Constructing psi_Rabi")
    print()
    
    # Project V_13 into flat band
    V_13_flat = P_flat.conj().T @ V_13 @ P_flat
    evals_13, evecs_13 = np.linalg.eigh(V_13_flat)
    
    # Max and min V_13 eigenstates within flat band
    idx_max = np.argmax(evals_13)
    idx_min = np.argmin(evals_13)
    
    psi_max_flat = evecs_13[:, idx_max]
    psi_min_flat = evecs_13[:, idx_min]
    
    # Transform to full space
    psi_max = P_flat @ psi_max_flat
    psi_min = P_flat @ psi_min_flat
    
    # Normalize
    psi_max = psi_max / np.linalg.norm(psi_max)
    psi_min = psi_min / np.linalg.norm(psi_min)
    
    # psi_Rabi = (|λ_max⟩ + |λ_min⟩)/√2
    psi_rabi = (psi_max + psi_min) / sqrt(2)
    psi_rabi = psi_rabi / np.linalg.norm(psi_rabi)
    
    print(f"  V_13 eigenvalues in flat band: [{evals_13.min():.4f}, {evals_13.max():.4f}]")
    print(f"  V_13 gap within flat band: {evals_13.max() - evals_13.min():.4f}")
    print()
    
    # Verify flat-band content
    fb_content = sum(np.abs(np.vdot(evecs_H[:, i], psi_rabi))**2 for i in flat_indices)
    print(f"  psi_Rabi flat-band content: {fb_content:.4f}")
    print()
    
    # ==========================================================================
    # STEP 3: Test V_13 Rabi oscillation starting from psi_max (not cavity!)
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 3: Rabi oscillation from |λ_max⟩ → |λ_min⟩")
    print("=" * 80)
    print()
    
    # Start from psi_max (one pole of the Rabi oscillation)
    psi_init = psi_max.copy()
    
    t_max = 100.0
    n_points = 2000
    times = np.linspace(0, t_max, n_points)
    gate_amplitude = 0.5
    
    H_gate = H + gate_amplitude * V_13
    evals_g, evecs_g = np.linalg.eigh(H_gate)
    psi_g = evecs_g.conj().T @ psi_init
    
    F_max = np.zeros(n_points)
    F_min = np.zeros(n_points)
    
    for k, t in enumerate(times):
        phases = np.exp(-1j * evals_g * t)
        psi_t = evecs_g @ (phases * psi_g)
        F_max[k] = np.abs(np.vdot(psi_max, psi_t))**2
        F_min[k] = np.abs(np.vdot(psi_min, psi_t))**2
    
    print(f"  Starting from |λ_max⟩:")
    print(f"  F_max range: [{F_max.min():.4f}, {F_max.max():.4f}], swing = {F_max.max()-F_max.min():.4f}")
    print(f"  F_min range: [{F_min.min():.4f}, {F_min.max():.4f}], swing = {F_min.max()-F_min.min():.4f}")
    print()
    
    # Find max transfer
    idx_flip = np.argmax(F_min)
    t_flip = times[idx_flip]
    transfer = F_min[idx_flip]
    
    print(f"  At t = {t_flip:.2f}: F_max = {F_max[idx_flip]:.4f}, F_min = {F_min[idx_flip]:.4f}")
    print(f"  Transfer efficiency: {transfer:.4f}")
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
    psi_n = evecs_n.conj().T @ psi_init
    
    F_max_noisy = np.zeros(n_points)
    
    for k, t in enumerate(times):
        phases = np.exp(-1j * evals_n * t)
        psi_t = evecs_n @ (phases * psi_n)
        F_max_noisy[k] = np.abs(np.vdot(psi_max, psi_t))**2
    
    max_diff = np.max(np.abs(F_max - F_max_noisy))
    avg_diff = np.mean(np.abs(F_max - F_max_noisy))
    
    print(f"  Max |F_clean - F_noisy|: {max_diff:.4f}")
    print(f"  Avg |F_clean - F_noisy|: {avg_diff:.4f}")
    print()
    
    # ==========================================================================
    # STEP 5: Compare with cavity state starting point
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 5: Compare — psi_Rabi (superposition) vs psi_max (eigenstate)")
    print("=" * 80)
    print()
    
    # Test starting from psi_Rabi (the superposition)
    psi_g_rabi = evecs_g.conj().T @ psi_rabi
    psi_n_rabi = evecs_n.conj().T @ psi_rabi
    
    F_rabi_clean = np.zeros(n_points)
    F_rabi_noisy = np.zeros(n_points)
    
    for k, t in enumerate(times):
        psi_t_c = evecs_g @ (np.exp(-1j * evals_g * t) * psi_g_rabi)
        psi_t_n = evecs_n @ (np.exp(-1j * evals_n * t) * psi_n_rabi)
        F_rabi_clean[k] = np.abs(np.vdot(psi_rabi, psi_t_c))**2
        F_rabi_noisy[k] = np.abs(np.vdot(psi_rabi, psi_t_n))**2
    
    swing_rabi = F_rabi_clean.max() - F_rabi_clean.min()
    max_diff_rabi = np.max(np.abs(F_rabi_clean - F_rabi_noisy))
    
    print(f"  psi_Rabi cavity fidelity:")
    print(f"    Swing: {swing_rabi:.4f}")
    print(f"    Max diff with noise: {max_diff_rabi:.4f}")
    print()
    
    # ==========================================================================
    # STEP 6: Test individual IN-prime screening on flat-band states
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 6: IN-prime screening on flat-band |λ_max⟩")
    print("=" * 80)
    print()
    
    for p, V_p in [(5, V_5), (7, V_7)]:
        H_p = H + gate_amplitude * V_p
        evals_p, evecs_p = np.linalg.eigh(H_p)
        psi_p = evecs_p.conj().T @ psi_max
        
        F_p = np.zeros(n_points)
        for k, t in enumerate(times):
            psi_t = evecs_p @ (np.exp(-1j * evals_p * t) * psi_p)
            F_p[k] = np.abs(np.vdot(psi_max, psi_t))**2
        
        swing_p = F_p.max() - F_p.min()
        print(f"  V_{p}: |λ_max⟩ fidelity swing = {swing_p:.4f}")
    
    print()
    
    # ==========================================================================
    # VERDICT
    # ==========================================================================
    
    print("=" * 80)
    print("  VERDICT")
    print("=" * 80)
    print()
    
    print(f"  Rabi from |λ_max⟩ → |λ_min⟩:")
    print(f"    Transfer: {transfer:.4f}")
    print(f"    Noise immunity (max diff): {max_diff:.4f}")
    print()
    
    print(f"  psi_Rabi cavity oscillation:")
    print(f"    Swing: {swing_rabi:.4f}")
    print(f"    Noise immunity (max diff): {max_diff_rabi:.4f}")
    print()
    
    if transfer > 0.9 and max_diff < 0.05:
        print("  ╔════════════════════════════════════════════════════════════════════╗")
        print("  ║  ✓ FLAT-BAND QUBIT CONFIRMED                                       ║")
        print("  ║  The V_13 eigenstates within flat band form a working qubit        ║")
        print("  ║  with noise immunity from IN-primes                                ║")
        print("  ╚════════════════════════════════════════════════════════════════════╝")
    elif transfer > 0.5 or max_diff < 0.1:
        print("  ╔════════════════════════════════════════════════════════════════════╗")
        print("  ║  ~ PARTIAL SUCCESS                                                 ║")
        print("  ║  Some rotation and some screening, but not complete                ║")
        print("  ╚════════════════════════════════════════════════════════════════════╝")
    else:
        print("  ╔════════════════════════════════════════════════════════════════════╗")
        print("  ║  ✗ FLAT-BAND QUBIT NOT WORKING                                     ║")
        print("  ║  Need different approach                                           ║")
        print("  ╚════════════════════════════════════════════════════════════════════╝")
    
    print()
    print("=" * 80)

if __name__ == "__main__":
    main()