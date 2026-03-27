#!/usr/bin/env python3
"""
THE INVARIANT SUBSPACE — Finding the True Connection
=====================================================

Key insight: V_p eigenstates are STATIONARY under V_p evolution.
We need states that are:
1. In a COMMON eigenspace of V_5, V_7, V_11 (screened)
2. NOT in a common eigenspace of V_13 (can rotate)

Strategy: Find the intersection structure.

Author: Claude + Tony (Fancyland LLC)
Date: 2026-03-25
"""

import numpy as np
from math import gcd, sqrt
from scipy.linalg import expm, block_diag
from itertools import product
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
    print("  THE INVARIANT SUBSPACE — Finding the True Connection")
    print("=" * 80)
    print()
    
    # Setup
    m = 2310
    alpha_c = sqrt(135/88)
    
    residues = get_coprime_residues(m)
    phi = len(residues)
    factors = factorize(m)
    in_primes = [p for p in factors if p > 2]  # Skip 2
    
    print(f"  m = {m}, φ = {phi}")
    print(f"  IN-primes (for screening): {in_primes}")
    print()
    
    # Build H(α_c)
    D_sym = build_D_sym(residues, m)
    P_tau = build_P_tau(residues, m)
    evals_D, _ = np.linalg.eigh(D_sym)
    lambda_P = np.max(evals_D)
    H = build_H_alpha(D_sym, P_tau, alpha_c, lambda_P)
    
    # ==========================================================================
    # KEY INSIGHT: The CRT structure creates BLOCK-DIAGONAL structure
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 1: CRT Block Structure Analysis")
    print("=" * 80)
    print()
    
    # For each residue r, compute (r mod 5, r mod 7, r mod 11)
    # This is a CRT signature
    
    crt_signatures = {}
    for i, r in enumerate(residues):
        sig = (r % 5, r % 7, r % 11)
        if sig not in crt_signatures:
            crt_signatures[sig] = []
        crt_signatures[sig].append(i)
    
    print(f"  Number of CRT blocks: {len(crt_signatures)}")
    print(f"  Block sizes: {sorted([len(v) for v in crt_signatures.values()])[:10]}...")
    print()
    
    # ==========================================================================
    # STEP 2: Find degenerate eigenspaces of (V_5, V_7, V_11) joint action
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 2: Joint eigenspaces of IN-primes")
    print("=" * 80)
    print()
    
    V_5 = build_V_diag(residues, 5)
    V_7 = build_V_diag(residues, 7)
    V_11 = build_V_diag(residues, 11)
    
    # For diagonal V_p, the eigenvalue of state |r⟩ is cos(2πr/p)
    # States with same (cos(2πr/5), cos(2πr/7), cos(2πr/11)) are in same joint eigenspace
    
    joint_eigenspaces = {}
    for i, r in enumerate(residues):
        sig = (round(np.cos(2*np.pi*r/5), 6), 
               round(np.cos(2*np.pi*r/7), 6), 
               round(np.cos(2*np.pi*r/11), 6))
        if sig not in joint_eigenspaces:
            joint_eigenspaces[sig] = []
        joint_eigenspaces[sig].append(i)
    
    print(f"  Number of joint eigenspaces: {len(joint_eigenspaces)}")
    
    # Find eigenspaces with dimension > 1 (these have internal degrees of freedom)
    degenerate_spaces = {k: v for k, v in joint_eigenspaces.items() if len(v) > 1}
    print(f"  Degenerate eigenspaces (dim > 1): {len(degenerate_spaces)}")
    
    if degenerate_spaces:
        dims = sorted([len(v) for v in degenerate_spaces.values()], reverse=True)
        print(f"  Largest degenerate subspace dimensions: {dims[:10]}")
    print()
    
    # ==========================================================================
    # STEP 3: Within degenerate subspaces, check V_13 action
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 3: V_13 action within degenerate subspaces")
    print("=" * 80)
    print()
    
    V_13 = build_V_diag(residues, 13)
    
    # Find the largest degenerate subspace
    largest_key = max(degenerate_spaces.keys(), key=lambda k: len(degenerate_spaces[k]))
    largest_indices = degenerate_spaces[largest_key]
    
    print(f"  Largest joint eigenspace: dim = {len(largest_indices)}")
    print(f"  IN-prime eigenvalues: {largest_key}")
    print()
    
    # Check V_13 eigenvalues within this subspace
    V_13_vals = [np.cos(2*np.pi*residues[i]/13) for i in largest_indices]
    print(f"  V_13 eigenvalues in this subspace: {sorted(set([round(v,4) for v in V_13_vals]))}")
    
    unique_13 = len(set([round(v, 6) for v in V_13_vals]))
    print(f"  Number of distinct V_13 eigenvalues: {unique_13}")
    print()
    
    if unique_13 > 1:
        print("  ✓ V_13 IS NOT CONSTANT on this subspace — potential for rotation!")
    else:
        print("  ✗ V_13 is constant — no rotation possible")
    
    print()
    
    # ==========================================================================
    # STEP 4: Construct qubit from degenerate subspace
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 4: Constructing qubit from degenerate subspace")
    print("=" * 80)
    print()
    
    # Find two states in largest degenerate subspace with DIFFERENT V_13 eigenvalues
    V_13_groups = {}
    for idx in largest_indices:
        val = round(np.cos(2*np.pi*residues[idx]/13), 6)
        if val not in V_13_groups:
            V_13_groups[val] = []
        V_13_groups[val].append(idx)
    
    print(f"  V_13 eigenvalue groups within degenerate subspace:")
    for val, indices in sorted(V_13_groups.items()):
        print(f"    λ_13 = {val:+.4f}: {len(indices)} states")
    
    print()
    
    if len(V_13_groups) >= 2:
        # Pick two groups with maximally different V_13 eigenvalues
        vals = sorted(V_13_groups.keys())
        val_0 = vals[-1]  # Maximum
        val_1 = vals[0]   # Minimum
        
        # Take first state from each group
        idx_0 = V_13_groups[val_0][0]
        idx_1 = V_13_groups[val_1][0]
        
        r_0 = residues[idx_0]
        r_1 = residues[idx_1]
        
        print(f"  Qubit states:")
        print(f"    |0⟩ = |r={r_0}⟩, V_13 eigenvalue = {val_0:+.4f}")
        print(f"    |1⟩ = |r={r_1}⟩, V_13 eigenvalue = {val_1:+.4f}")
        print(f"    V_13 gap = {val_0 - val_1:.4f}")
        print()
        
        # Verify they have SAME IN-prime eigenvalues
        print(f"  Verification — same IN-prime eigenvalues:")
        for p in [5, 7, 11]:
            v0 = np.cos(2*np.pi*r_0/p)
            v1 = np.cos(2*np.pi*r_1/p)
            print(f"    V_{p}: |0⟩ = {v0:+.4f}, |1⟩ = {v1:+.4f}, diff = {abs(v0-v1):.6f}")
        
        print()
        
        # ==========================================================================
        # STEP 5: Test this qubit
        # ==========================================================================
        
        print("=" * 80)
        print("  STEP 5: Testing the degenerate subspace qubit")
        print("=" * 80)
        print()
        
        # Create basis states
        psi_0 = np.zeros(phi, dtype=complex)
        psi_0[idx_0] = 1.0
        
        psi_1 = np.zeros(phi, dtype=complex)
        psi_1[idx_1] = 1.0
        
        # Cavity state
        psi_cav = (psi_0 + psi_1) / sqrt(2)
        
        # Test evolution
        t_max = 50.0
        n_points = 1000
        gate_amplitude = 0.5
        
        # Under V_5 (should be screened — both have same eigenvalue)
        H_V5 = H + gate_amplitude * V_5
        evals_5, evecs_5 = np.linalg.eigh(H_V5)
        psi_5 = evecs_5.conj().T @ psi_cav
        
        times = np.linspace(0, t_max, n_points)
        F_5 = np.zeros(n_points)
        for i, t in enumerate(times):
            phases = np.exp(-1j * evals_5 * t)
            psi_t = evecs_5 @ (phases * psi_5)
            F_5[i] = np.abs(np.vdot(psi_cav, psi_t))**2
        
        swing_5 = F_5.max() - F_5.min()
        
        # Under V_13 (should rotate — different eigenvalues)
        H_V13 = H + gate_amplitude * V_13
        evals_13, evecs_13 = np.linalg.eigh(H_V13)
        psi_13 = evecs_13.conj().T @ psi_cav
        
        F_13 = np.zeros(n_points)
        F_0_trace = np.zeros(n_points)
        F_1_trace = np.zeros(n_points)
        for i, t in enumerate(times):
            phases = np.exp(-1j * evals_13 * t)
            psi_t = evecs_13 @ (phases * psi_13)
            F_13[i] = np.abs(np.vdot(psi_cav, psi_t))**2
            F_0_trace[i] = np.abs(np.vdot(psi_0, psi_t))**2
            F_1_trace[i] = np.abs(np.vdot(psi_1, psi_t))**2
        
        swing_13 = F_13.max() - F_13.min()
        
        print(f"  Evolution under V_5 (IN-prime):")
        print(f"    Fidelity swing: {swing_5:.4f}")
        print(f"    Status: {'✓ SCREENED' if swing_5 < 0.1 else '✗ NOT SCREENED'}")
        print()
        
        print(f"  Evolution under V_13 (OUT-prime):")
        print(f"    Fidelity swing: {swing_13:.4f}")
        print(f"    Status: {'✓ ROTATES' if swing_13 > 0.5 else '~ PARTIAL' if swing_13 > 0.1 else '✗ NO ROTATION'}")
        print()
        
        # Check Rabi oscillation
        print(f"  Rabi oscillation |0⟩ ↔ |1⟩:")
        print(f"    F_0 range: [{F_0_trace.min():.4f}, {F_0_trace.max():.4f}]")
        print(f"    F_1 range: [{F_1_trace.min():.4f}, {F_1_trace.max():.4f}]")
        
        # Find time of maximum transfer
        t_flip = times[np.argmax(F_1_trace)]
        F0_at_flip = F_0_trace[np.argmax(F_1_trace)]
        F1_at_flip = F_1_trace[np.argmax(F_1_trace)]
        
        print(f"    At t = {t_flip:.2f}: F_0 = {F0_at_flip:.4f}, F_1 = {F1_at_flip:.4f}")
        print()
        
        # ==========================================================================
        # STEP 6: Noise immunity test
        # ==========================================================================
        
        print("=" * 80)
        print("  STEP 6: Noise immunity during V_13 gate")
        print("=" * 80)
        print()
        
        noise_amp = 0.1
        H_clean = H + gate_amplitude * V_13
        H_noisy = H + gate_amplitude * V_13 + noise_amp * (V_5 + V_7)
        
        evals_c, evecs_c = np.linalg.eigh(H_clean)
        evals_n, evecs_n = np.linalg.eigh(H_noisy)
        
        psi_c = evecs_c.conj().T @ psi_cav
        psi_n = evecs_n.conj().T @ psi_cav
        
        F_clean = np.zeros(n_points)
        F_noisy = np.zeros(n_points)
        
        for i, t in enumerate(times):
            psi_t_c = evecs_c @ (np.exp(-1j * evals_c * t) * psi_c)
            psi_t_n = evecs_n @ (np.exp(-1j * evals_n * t) * psi_n)
            F_clean[i] = np.abs(np.vdot(psi_cav, psi_t_c))**2
            F_noisy[i] = np.abs(np.vdot(psi_cav, psi_t_n))**2
        
        max_diff = np.max(np.abs(F_clean - F_noisy))
        avg_diff = np.mean(np.abs(F_clean - F_noisy))
        
        print(f"  Max |F_clean - F_noisy|: {max_diff:.4f}")
        print(f"  Avg |F_clean - F_noisy|: {avg_diff:.4f}")
        print()
        
        if swing_5 < 0.1 and swing_13 > 0.5 and max_diff < 0.1:
            print("  ╔════════════════════════════════════════════════════════════════════╗")
            print("  ║  CONNECTION FOUND!                                                 ║")
            print("  ║  Degenerate subspace qubit: screened by IN, rotated by OUT         ║")
            print("  ╚════════════════════════════════════════════════════════════════════╝")
        elif swing_5 < 0.3 and swing_13 > 0.3:
            print("  ╔════════════════════════════════════════════════════════════════════╗")
            print("  ║  PARTIAL CONNECTION — The structure exists but needs refinement    ║")
            print("  ╚════════════════════════════════════════════════════════════════════╝")
        else:
            print("  ╔════════════════════════════════════════════════════════════════════╗")
            print("  ║  CONNECTION NOT CLEAN — Need different encoding                    ║")
            print("  ╚════════════════════════════════════════════════════════════════════╝")
    
    else:
        print("  ✗ No degenerate subspace found with multiple V_13 eigenvalues")
    
    print()
    
    # ==========================================================================
    # STEP 7: The Hamiltonian connection
    # ==========================================================================
    
    print("=" * 80)
    print("  STEP 7: Does the Hamiltonian H(α) preserve the structure?")
    print("=" * 80)
    print()
    
    # The key question: do |r_0⟩ and |r_1⟩ stay in their joint eigenspace under H?
    # This requires checking if H preserves the CRT blocks
    
    # Check commutator [H, V_5]
    comm_H_V5 = H @ V_5 - V_5 @ H
    comm_norm_5 = np.linalg.norm(comm_H_V5)
    
    comm_H_V13 = H @ V_13 - V_13 @ H
    comm_norm_13 = np.linalg.norm(comm_H_V13)
    
    print(f"  ||[H, V_5]|| = {comm_norm_5:.4f}")
    print(f"  ||[H, V_13]|| = {comm_norm_13:.4f}")
    print()
    
    if comm_norm_5 < 1.0:
        print("  H approximately commutes with V_5 — structure preserved!")
    else:
        print("  H does NOT commute with V_5 — H mixes the joint eigenspaces")
        print("  This is why the basis states |r⟩ scatter under H evolution")
    
    print()
    print("=" * 80)

if __name__ == "__main__":
    main()