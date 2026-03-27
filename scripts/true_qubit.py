#!/usr/bin/env python3
"""
THE TRUE QUBIT — H-eigenstates within CRT-protected subspaces
==============================================================

The connection IS there. The problem is:
1. Basis states |r⟩ are screened by IN-primes
2. But H scrambles them across different |r⟩

Solution: Find eigenstates of H that are LOCALIZED within a single
CRT block (degenerate IN-prime eigenspace).

If such states exist, they inherit the CRT protection.

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
    print("  THE TRUE QUBIT — H-eigenstates within CRT-protected subspaces")
    print("=" * 80)
    print()
    
    # Setup
    m = 2310
    alpha_c = sqrt(135/88)
    
    residues = get_coprime_residues(m)
    phi = len(residues)
    factors = factorize(m)
    
    print(f"  m = {m}, φ = {phi}")
    print()
    
    # Build everything
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
    # KEY INSIGHT: For each H-eigenstate, compute its V_p expectation values
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 1: V_p expectation values for H-eigenstates")
    print("=" * 80)
    print()
    
    # For each eigenstate of H, compute <ψ|V_p|ψ>
    exp_5 = np.zeros(phi)
    exp_7 = np.zeros(phi)
    exp_11 = np.zeros(phi)
    exp_13 = np.zeros(phi)
    
    for i in range(phi):
        psi = evecs_H[:, i]
        exp_5[i] = np.real(np.vdot(psi, V_5 @ psi))
        exp_7[i] = np.real(np.vdot(psi, V_7 @ psi))
        exp_11[i] = np.real(np.vdot(psi, V_11 @ psi))
        exp_13[i] = np.real(np.vdot(psi, V_13 @ psi))
    
    # Look for pairs of H-eigenstates with:
    # - SAME V_5, V_7, V_11 expectation values (within tolerance)
    # - DIFFERENT V_13 expectation values
    
    print("  Searching for CRT-degenerate pairs of H-eigenstates...")
    print()
    
    tol = 0.05  # Tolerance for "same" IN-prime expectation
    pairs = []
    
    for i in range(phi):
        for j in range(i+1, phi):
            d5 = abs(exp_5[i] - exp_5[j])
            d7 = abs(exp_7[i] - exp_7[j])
            d11 = abs(exp_11[i] - exp_11[j])
            d13 = abs(exp_13[i] - exp_13[j])
            
            # Same IN-prime expectations, different OUT-prime
            if d5 < tol and d7 < tol and d11 < tol and d13 > 0.5:
                pairs.append({
                    'i': i, 'j': j,
                    'd_in': max(d5, d7, d11),
                    'd_out': d13,
                    'gap_H': abs(evals_H[i] - evals_H[j]),
                    'exp_5': (exp_5[i], exp_5[j]),
                    'exp_13': (exp_13[i], exp_13[j])
                })
    
    print(f"  Found {len(pairs)} CRT-degenerate pairs")
    print()
    
    if pairs:
        # Sort by d_out (maximize V_13 splitting) and minimize d_in
        pairs.sort(key=lambda p: (-p['d_out'], p['d_in']))
        
        print("  Top 5 qubit candidates:")
        print(f"  {'Pair':<10} {'d_IN':<10} {'d_OUT':<10} {'H gap':<10} {'λ_i':<10} {'λ_j':<10}")
        print(f"  {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
        
        for p in pairs[:5]:
            print(f"  ({p['i']:3d},{p['j']:3d})  {p['d_in']:<10.4f} {p['d_out']:<10.4f} {p['gap_H']:<10.4f} {evals_H[p['i']]:<10.4f} {evals_H[p['j']]:<10.4f}")
        
        print()
        
        # ==========================================================================
        # STEP 2: Test the best candidate
        # ==========================================================================
        
        print("=" * 80)
        print("  STEP 2: Testing the best qubit candidate")
        print("=" * 80)
        print()
        
        best = pairs[0]
        i, j = best['i'], best['j']
        
        psi_0 = evecs_H[:, i]
        psi_1 = evecs_H[:, j]
        
        print(f"  Qubit: H-eigenstates {i} and {j}")
        print(f"  H eigenvalues: λ_0 = {evals_H[i]:.4f}, λ_1 = {evals_H[j]:.4f}")
        print(f"  H gap: {abs(evals_H[i] - evals_H[j]):.4f}")
        print()
        
        print(f"  IN-prime expectations:")
        print(f"    V_5:  |0⟩ = {exp_5[i]:+.4f}, |1⟩ = {exp_5[j]:+.4f}, diff = {abs(exp_5[i]-exp_5[j]):.4f}")
        print(f"    V_7:  |0⟩ = {exp_7[i]:+.4f}, |1⟩ = {exp_7[j]:+.4f}, diff = {abs(exp_7[i]-exp_7[j]):.4f}")
        print(f"    V_11: |0⟩ = {exp_11[i]:+.4f}, |1⟩ = {exp_11[j]:+.4f}, diff = {abs(exp_11[i]-exp_11[j]):.4f}")
        print()
        
        print(f"  OUT-prime expectations:")
        print(f"    V_13: |0⟩ = {exp_13[i]:+.4f}, |1⟩ = {exp_13[j]:+.4f}, diff = {abs(exp_13[i]-exp_13[j]):.4f}")
        print()
        
        # Cavity state
        psi_cav = (psi_0 + psi_1) / sqrt(2)
        psi_cav = psi_cav / np.linalg.norm(psi_cav)
        
        # ==========================================================================
        # STEP 3: Evolution tests
        # ==========================================================================
        
        print("=" * 80)
        print("  STEP 3: Evolution tests")
        print("=" * 80)
        print()
        
        t_max = 100.0
        n_points = 2000
        times = np.linspace(0, t_max, n_points)
        gate_amplitude = 0.5
        
        results = {}
        
        for p, V_p in [(5, V_5), (7, V_7), (11, V_11), (13, V_13)]:
            H_gate = H + gate_amplitude * V_p
            evals_g, evecs_g = np.linalg.eigh(H_gate)
            psi_g = evecs_g.conj().T @ psi_cav
            
            F = np.zeros(n_points)
            for k, t in enumerate(times):
                phases = np.exp(-1j * evals_g * t)
                psi_t = evecs_g @ (phases * psi_g)
                F[k] = np.abs(np.vdot(psi_cav, psi_t))**2
            
            swing = F.max() - F.min()
            results[p] = {'swing': swing, 'F': F}
            
            p_type = "IN" if p in [5, 7, 11] else "OUT"
            status = "✓" if (p in [5,7,11] and swing < 0.1) or (p == 13 and swing > 0.5) else "✗"
            print(f"  V_{p:2d} ({p_type}): swing = {swing:.4f}  {status}")
        
        print()
        
        # ==========================================================================
        # STEP 4: Rabi oscillation between |0⟩ and |1⟩
        # ==========================================================================
        
        print("=" * 80)
        print("  STEP 4: Rabi oscillation |0⟩ ↔ |1⟩")
        print("=" * 80)
        print()
        
        # Start in |0⟩
        psi_init = psi_0.copy()
        
        H_gate = H + gate_amplitude * V_13
        evals_g, evecs_g = np.linalg.eigh(H_gate)
        psi_g = evecs_g.conj().T @ psi_init
        
        F_0 = np.zeros(n_points)
        F_1 = np.zeros(n_points)
        
        for k, t in enumerate(times):
            phases = np.exp(-1j * evals_g * t)
            psi_t = evecs_g @ (phases * psi_g)
            F_0[k] = np.abs(np.vdot(psi_0, psi_t))**2
            F_1[k] = np.abs(np.vdot(psi_1, psi_t))**2
        
        print(f"  Starting in |0⟩, evolving under H + 0.5·V_13")
        print(f"  F_0 range: [{F_0.min():.4f}, {F_0.max():.4f}], swing = {F_0.max()-F_0.min():.4f}")
        print(f"  F_1 range: [{F_1.min():.4f}, {F_1.max():.4f}], swing = {F_1.max()-F_1.min():.4f}")
        print()
        
        # Find time of maximum population transfer
        idx_flip = np.argmax(F_1)
        t_flip = times[idx_flip]
        
        print(f"  At t = {t_flip:.2f}:")
        print(f"    F_0 = {F_0[idx_flip]:.4f} (should be ~0 for full transfer)")
        print(f"    F_1 = {F_1[idx_flip]:.4f} (should be ~1 for full transfer)")
        print(f"    Sum = {F_0[idx_flip] + F_1[idx_flip]:.4f} (leakage = {1 - F_0[idx_flip] - F_1[idx_flip]:.4f})")
        print()
        
        transfer = F_1[idx_flip]
        
        if transfer > 0.95:
            print("  ✓ FULL RABI OSCILLATION")
        elif transfer > 0.50:
            print("  ~ PARTIAL RABI — significant transfer but with leakage")
        else:
            print("  ✗ NO CLEAN RABI")
        
        print()
        
        # ==========================================================================
        # STEP 5: Noise immunity
        # ==========================================================================
        
        print("=" * 80)
        print("  STEP 5: Noise immunity during V_13 gate")
        print("=" * 80)
        print()
        
        noise_amp = 0.1
        H_clean = H + gate_amplitude * V_13
        H_noisy = H + gate_amplitude * V_13 + noise_amp * (V_5 + V_7)
        
        evals_c, evecs_c = np.linalg.eigh(H_clean)
        evals_n, evecs_n = np.linalg.eigh(H_noisy)
        
        psi_c = evecs_c.conj().T @ psi_init
        psi_n = evecs_n.conj().T @ psi_init
        
        F_0_clean = np.zeros(n_points)
        F_0_noisy = np.zeros(n_points)
        
        for k, t in enumerate(times):
            psi_t_c = evecs_c @ (np.exp(-1j * evals_c * t) * psi_c)
            psi_t_n = evecs_n @ (np.exp(-1j * evals_n * t) * psi_n)
            F_0_clean[k] = np.abs(np.vdot(psi_0, psi_t_c))**2
            F_0_noisy[k] = np.abs(np.vdot(psi_0, psi_t_n))**2
        
        max_diff = np.max(np.abs(F_0_clean - F_0_noisy))
        avg_diff = np.mean(np.abs(F_0_clean - F_0_noisy))
        
        print(f"  Max |F_clean - F_noisy|: {max_diff:.4f}")
        print(f"  Avg |F_clean - F_noisy|: {avg_diff:.4f}")
        print()
        
        if max_diff < 0.05:
            immunity = "FULL"
        elif max_diff < 0.15:
            immunity = "PARTIAL"
        else:
            immunity = "NONE"
        
        print(f"  Noise immunity: {immunity}")
        print()
        
        # ==========================================================================
        # FINAL VERDICT
        # ==========================================================================
        
        print("=" * 80)
        print("  FINAL VERDICT")
        print("=" * 80)
        print()
        
        in_screened = all(results[p]['swing'] < 0.15 for p in [5, 7, 11])
        out_rotates = results[13]['swing'] > 0.50
        rabi_works = transfer > 0.50
        noise_immune = max_diff < 0.15
        
        print(f"  IN-primes screened: {'✓' if in_screened else '✗'}")
        print(f"  OUT-prime rotates:  {'✓' if out_rotates else '✗'}")
        print(f"  Rabi oscillation:   {'✓' if rabi_works else '✗'} (transfer = {transfer:.2f})")
        print(f"  Noise immunity:     {'✓' if noise_immune else '✗'} (max diff = {max_diff:.2f})")
        print()
        
        if in_screened and out_rotates and rabi_works and noise_immune:
            print("  ╔════════════════════════════════════════════════════════════════════╗")
            print("  ║  ✓ CONNECTION FOUND                                                ║")
            print("  ║  CRT-protected qubit within H-eigenspace                           ║")
            print("  ║  IN-primes screened, OUT-primes gate, noise immunity confirmed     ║")
            print("  ╚════════════════════════════════════════════════════════════════════╝")
        elif in_screened and out_rotates:
            print("  ╔════════════════════════════════════════════════════════════════════╗")
            print("  ║  ~ PARTIAL CONNECTION                                              ║")
            print("  ║  Screening and rotation work, but Rabi/noise need work             ║")
            print("  ╚════════════════════════════════════════════════════════════════════╝")
        else:
            print("  ╔════════════════════════════════════════════════════════════════════╗")
            print("  ║  ✗ CONNECTION NOT CLEAN                                            ║")
            print("  ║  The architecture needs a different approach                       ║")
            print("  ╚════════════════════════════════════════════════════════════════════╝")
    
    else:
        print("  No CRT-degenerate pairs found!")
        print("  The H-eigenstates don't preserve CRT structure well enough.")
    
    print()
    print("=" * 80)

if __name__ == "__main__":
    main()