#!/usr/bin/env python3
"""
FINDING THE CONNECTION — Qubit States WITHIN the Flat Band
===========================================================

The flat band (465D) is screened against IN-primes.
The edge modes (Perron/Anti-Perron) are NOT screened.

Question: Can we find a qubit WITHIN the flat band that:
1. Is screened against IN-primes (V_5, V_7, V_11)
2. Can be rotated by OUT-primes (V_13, V_17, etc.)

If yes → The connection exists. Arithmetic quantum computing is real.

Author: Claude + Tony (Fancyland LLC)
Date: 2026-03-25
"""

import numpy as np
from math import gcd, sqrt
from scipy.linalg import expm, null_space
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CORE CONSTRUCTION (same as before)
# =============================================================================

def get_coprime_residues(m):
    return sorted([r for r in range(1, m) if gcd(r, m) == 1])

def factorize(n):
    factors = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            if d not in factors:
                factors.append(d)
            n //= d
        d += 1
    if n > 1 and n not in factors:
        factors.append(n)
    return sorted(factors)

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
# MAIN ANALYSIS
# =============================================================================

def main():
    print("=" * 80)
    print("  FINDING THE CONNECTION — Qubit States Within the Flat Band")
    print("=" * 80)
    print()
    
    # Setup
    m = 2310
    alpha_c = sqrt(135/88)
    
    residues = get_coprime_residues(m)
    phi = len(residues)
    factors = factorize(m)
    
    print(f"  m = {m}, φ = {phi}")
    print(f"  IN-primes: {factors}")
    print()
    
    # Build H(α_c)
    D_sym = build_D_sym(residues, m)
    P_tau = build_P_tau(residues, m)
    evals_D, _ = np.linalg.eigh(D_sym)
    lambda_P = np.max(evals_D)
    H = build_H_alpha(D_sym, P_tau, alpha_c, lambda_P)
    evals_H, evecs_H = np.linalg.eigh(H)
    
    # Find flat band (eigenvalues clustered near 0)
    flat_threshold = 0.1
    flat_indices = np.where(np.abs(evals_H) < flat_threshold)[0]
    flat_dim = len(flat_indices)
    
    print(f"  Flat band dimension: {flat_dim} (threshold |λ| < {flat_threshold})")
    print(f"  Flat band eigenvalue range: [{evals_H[flat_indices].min():.6f}, {evals_H[flat_indices].max():.6f}]")
    print()
    
    # Extract flat band subspace
    P_flat = evecs_H[:, flat_indices]  # Projector columns
    
    # ==========================================================================
    # KEY INSIGHT: Project V_p into the flat band and find ITS eigenstates
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 1: V_p eigenstructure WITHIN the flat band")
    print("=" * 80)
    print()
    
    test_primes = [5, 7, 11, 13, 17, 19]
    
    for p in test_primes:
        V_p = build_V_diag(residues, p)
        
        # Project V_p into flat band: V_p_flat = P_flat^† V_p P_flat
        V_p_flat = P_flat.conj().T @ V_p @ P_flat
        
        # Eigendecompose V_p restricted to flat band
        evals_Vp, evecs_Vp = np.linalg.eigh(V_p_flat)
        
        eval_range = evals_Vp.max() - evals_Vp.min()
        eval_std = np.std(evals_Vp)
        
        p_type = "IN" if p in factors else "OUT"
        
        print(f"  V_{p:2d} ({p_type}): eigenvalue range = {eval_range:.4f}, std = {eval_std:.4f}")
        print(f"         min = {evals_Vp.min():.4f}, max = {evals_Vp.max():.4f}")
    
    print()
    
    # ==========================================================================
    # STEP 2: Find pairs of flat-band states with MAXIMAL V_p splitting
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 2: Finding qubit candidates in flat band")
    print("=" * 80)
    print()
    
    # For V_13 (OUT-prime), find the two most separated eigenstates
    V_13 = build_V_diag(residues, 13)
    V_13_flat = P_flat.conj().T @ V_13 @ P_flat
    evals_13, evecs_13 = np.linalg.eigh(V_13_flat)
    
    # Max and min eigenstates of V_13 within flat band
    idx_max = np.argmax(evals_13)
    idx_min = np.argmin(evals_13)
    
    # These are in flat-band coordinates, transform back to full space
    psi_plus_flat = evecs_13[:, idx_max]  # Flat-band coords
    psi_minus_flat = evecs_13[:, idx_min]
    
    psi_plus = P_flat @ psi_plus_flat  # Full space
    psi_minus = P_flat @ psi_minus_flat
    
    gap_13 = evals_13[idx_max] - evals_13[idx_min]
    
    print(f"  V_13 eigenstates in flat band:")
    print(f"    λ_max = {evals_13[idx_max]:.6f}")
    print(f"    λ_min = {evals_13[idx_min]:.6f}")
    print(f"    Gap = {gap_13:.6f}")
    print()
    
    # Define flat-band qubit: |0⟩ = |λ_max⟩, |1⟩ = |λ_min⟩
    psi_0 = psi_plus / np.linalg.norm(psi_plus)
    psi_1 = psi_minus / np.linalg.norm(psi_minus)
    
    # Verify orthogonality
    overlap = np.abs(np.vdot(psi_0, psi_1))**2
    print(f"  |⟨0|1⟩|² = {overlap:.6f} (should be ~0)")
    print()
    
    # ==========================================================================
    # STEP 3: Test screening of this flat-band qubit
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 3: Screening test on flat-band qubit")
    print("=" * 80)
    print()
    
    # Cavity state within flat band
    psi_cavity_fb = (psi_0 + psi_1) / sqrt(2)
    psi_cavity_fb = psi_cavity_fb / np.linalg.norm(psi_cavity_fb)
    
    # Check: is this actually in the flat band?
    flat_band_projection = np.sum([np.abs(np.vdot(evecs_H[:, i], psi_cavity_fb))**2 
                                   for i in flat_indices])
    print(f"  Cavity state flat-band content: {flat_band_projection:.4f} (should be 1.0)")
    print()
    
    # Test evolution under each V_p
    t_max = 50.0
    n_points = 500
    gate_amplitude = 0.5
    
    print(f"  {'Prime':<8} {'Type':<6} {'Swing':<10} {'Status'}")
    print(f"  {'-'*8} {'-'*6} {'-'*10} {'-'*20}")
    
    for p in [5, 7, 11, 13, 17, 19, 23]:
        V_p = build_V_diag(residues, p)
        H_gate = H + gate_amplitude * V_p
        
        # Diagonalize
        evals_gate, evecs_gate = np.linalg.eigh(H_gate)
        psi_0_eig = evecs_gate.conj().T @ psi_cavity_fb
        
        # Time evolution
        times = np.linspace(0, t_max, n_points)
        F = np.zeros(n_points)
        for i, t in enumerate(times):
            phases = np.exp(-1j * evals_gate * t)
            psi_t_eig = phases * psi_0_eig
            psi_t = evecs_gate @ psi_t_eig
            F[i] = np.abs(np.vdot(psi_cavity_fb, psi_t))**2
        
        swing = np.max(F) - np.min(F)
        p_type = "IN" if p in factors else "OUT"
        
        if p in factors:
            status = "✓ SCREENED" if swing < 0.10 else "✗ NOT SCREENED"
        else:
            status = "✓ ROTATES" if swing > 0.50 else "✗ NO ROTATION"
        
        print(f"  V_{p:<5}  {p_type:<6} {swing:<10.4f} {status}")
    
    print()
    
    # ==========================================================================
    # STEP 4: The definitive test — IN vs OUT discrimination
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 4: IN/OUT discrimination on flat-band qubit")
    print("=" * 80)
    print()
    
    # Measure V_p expectation values on psi_0 and psi_1
    print(f"  {'Prime':<8} {'Type':<6} {'⟨0|V_p|0⟩':<12} {'⟨1|V_p|1⟩':<12} {'Difference':<12}")
    print(f"  {'-'*8} {'-'*6} {'-'*12} {'-'*12} {'-'*12}")
    
    for p in [5, 7, 11, 13, 17, 19]:
        V_p = build_V_diag(residues, p)
        
        exp_0 = np.real(np.vdot(psi_0, V_p @ psi_0))
        exp_1 = np.real(np.vdot(psi_1, V_p @ psi_1))
        diff = abs(exp_0 - exp_1)
        
        p_type = "IN" if p in factors else "OUT"
        
        print(f"  V_{p:<5}  {p_type:<6} {exp_0:<12.4f} {exp_1:<12.4f} {diff:<12.4f}")
    
    print()
    
    # ==========================================================================
    # STEP 5: Rabi oscillation on flat-band qubit
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 5: Rabi oscillation |0⟩ ↔ |1⟩ via V_13")
    print("=" * 80)
    print()
    
    # Start in |0⟩ (V_13 max eigenstate in flat band)
    psi_init = psi_0.copy()
    
    V_13 = build_V_diag(residues, 13)
    H_gate = H + gate_amplitude * V_13
    
    # Diagonalize
    evals_gate, evecs_gate = np.linalg.eigh(H_gate)
    psi_init_eig = evecs_gate.conj().T @ psi_init
    
    # Time evolution
    t_max = 100.0
    n_points = 2000
    times = np.linspace(0, t_max, n_points)
    
    F_0 = np.zeros(n_points)  # Overlap with |0⟩
    F_1 = np.zeros(n_points)  # Overlap with |1⟩
    
    for i, t in enumerate(times):
        phases = np.exp(-1j * evals_gate * t)
        psi_t_eig = phases * psi_init_eig
        psi_t = evecs_gate @ psi_t_eig
        F_0[i] = np.abs(np.vdot(psi_0, psi_t))**2
        F_1[i] = np.abs(np.vdot(psi_1, psi_t))**2
    
    print(f"  Starting in |0⟩, evolving under H + 0.5·V_13")
    print(f"  F_0 range: [{F_0.min():.4f}, {F_0.max():.4f}]")
    print(f"  F_1 range: [{F_1.min():.4f}, {F_1.max():.4f}]")
    print(f"  F_0 swing: {F_0.max() - F_0.min():.4f}")
    print(f"  F_1 swing: {F_1.max() - F_1.min():.4f}")
    print()
    
    # Check for clean Rabi between |0⟩ and |1⟩
    # At some time, F_0 should go to 0 while F_1 goes to 1
    
    max_F1 = F_1.max()
    min_F0 = F_0.min()
    t_flip = times[np.argmax(F_1)]
    
    print(f"  At t = {t_flip:.2f}:")
    print(f"    F_0 = {F_0[np.argmax(F_1)]:.4f} (should be ~0)")
    print(f"    F_1 = {max_F1:.4f} (should be ~1)")
    print()
    
    if max_F1 > 0.90 and min_F0 < 0.10:
        print("  ✓ CLEAN RABI OSCILLATION WITHIN FLAT BAND")
    elif max_F1 > 0.50:
        print("  ~ PARTIAL RABI — need to investigate leakage")
    else:
        print("  ✗ NO CLEAN RABI")
    
    print()
    
    # ==========================================================================
    # STEP 6: The ultimate test — screening DURING gate operation
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 6: Screening during V_13 gate on flat-band qubit")
    print("=" * 80)
    print()
    
    # Compare V_13 gate with and without V_5+V_7 noise
    V_5 = build_V_diag(residues, 5)
    V_7 = build_V_diag(residues, 7)
    
    noise_amplitude = 0.1
    
    H_clean = H + gate_amplitude * V_13
    H_noisy = H + gate_amplitude * V_13 + noise_amplitude * (V_5 + V_7)
    
    # Clean evolution
    evals_c, evecs_c = np.linalg.eigh(H_clean)
    psi_c = evecs_c.conj().T @ psi_init
    
    # Noisy evolution
    evals_n, evecs_n = np.linalg.eigh(H_noisy)
    psi_n = evecs_n.conj().T @ psi_init
    
    F_clean = np.zeros(n_points)
    F_noisy = np.zeros(n_points)
    
    for i, t in enumerate(times):
        # Clean
        phases_c = np.exp(-1j * evals_c * t)
        psi_t_c = evecs_c @ (phases_c * psi_c)
        F_clean[i] = np.abs(np.vdot(psi_0, psi_t_c))**2
        
        # Noisy
        phases_n = np.exp(-1j * evals_n * t)
        psi_t_n = evecs_n @ (phases_n * psi_n)
        F_noisy[i] = np.abs(np.vdot(psi_0, psi_t_n))**2
    
    max_diff = np.max(np.abs(F_clean - F_noisy))
    avg_diff = np.mean(np.abs(F_clean - F_noisy))
    
    print(f"  Clean V_13 gate: swing = {F_clean.max() - F_clean.min():.4f}")
    print(f"  Noisy V_13 gate: swing = {F_noisy.max() - F_noisy.min():.4f}")
    print(f"  Max |F_clean - F_noisy|: {max_diff:.4f}")
    print(f"  Avg |F_clean - F_noisy|: {avg_diff:.4f}")
    print()
    
    if max_diff < 0.05:
        print("  ✓ IN-PRIME NOISE FULLY SCREENED DURING GATE")
        print("  ╔════════════════════════════════════════════════════════════════════╗")
        print("  ║  CONNECTION FOUND: FLAT-BAND QUBIT WITH ARITHMETIC IMMUNITY        ║")
        print("  ╚════════════════════════════════════════════════════════════════════╝")
    elif max_diff < 0.15:
        print("  ~ PARTIAL SCREENING — connection exists but imperfect")
    else:
        print("  ✗ NO SCREENING — connection not found in this encoding")
    
    print()
    print("=" * 80)

if __name__ == "__main__":
    main()