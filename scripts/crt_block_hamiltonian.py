#!/usr/bin/env python3
"""
THE CRT BLOCK HAMILTONIAN — The Third Option
==============================================

The wall: H scrambles CRT blocks (||[H,V_5]|| = 1.21).
The insight: WITHIN a CRT block, IN-primes are EXACTLY identity.

Strategy:
  1. Identify the 30 CRT joint eigenspaces (dim 16 each)
  2. Project H into a single 16D block: H_block = P_block^† H P_block
  3. Within this block, V_5 = V_7 = V_11 = constant (by construction!)
  4. Project V_13 into same block: V_13_block has 7 distinct eigenvalues
  5. Test: does H_block + V_13_block produce a protected qubit?

If H_block is nontrivial (not proportional to identity), then:
  - It generates dynamics WITHIN the protected subspace
  - IN-prime screening is EXACT (not approximate)
  - V_13 can rotate between H_block eigenstates

This is the Third Option: not asking H to commute with CRT,
but RESTRICTING H to where CRT protection is already guaranteed.

Author: Claude + Tony (Fancyland LLC)
Date: 2026-03-25
"""

import numpy as np
from math import gcd, sqrt
from scipy.linalg import expm
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
    print("  THE CRT BLOCK HAMILTONIAN — The Third Option")
    print("=" * 80)
    print()

    # Setup
    m = 2310
    alpha_c = sqrt(135/88)

    residues = get_coprime_residues(m)
    phi = len(residues)

    print(f"  m = {m}, φ = {phi}")
    print()

    # Build full operators
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
    V_17 = build_V_diag(residues, 17)
    V_19 = build_V_diag(residues, 19)

    # =========================================================================
    # STEP 1: Build CRT joint eigenspaces
    # =========================================================================

    print("=" * 80)
    print("  STEP 1: CRT Joint Eigenspaces of (V_5, V_7, V_11)")
    print("=" * 80)
    print()

    # Group residues by their (cos(2πr/5), cos(2πr/7), cos(2πr/11)) signature
    joint_eigenspaces = {}
    for i, r in enumerate(residues):
        sig = (round(np.cos(2*np.pi*r/5), 8),
               round(np.cos(2*np.pi*r/7), 8),
               round(np.cos(2*np.pi*r/11), 8))
        if sig not in joint_eigenspaces:
            joint_eigenspaces[sig] = []
        joint_eigenspaces[sig].append(i)

    print(f"  Number of joint eigenspaces: {len(joint_eigenspaces)}")
    dims = sorted([len(v) for v in joint_eigenspaces.values()], reverse=True)
    print(f"  Dimensions: {dims}")
    print()

    # =========================================================================
    # STEP 2: Pick the largest block and project H into it
    # =========================================================================

    print("=" * 80)
    print("  STEP 2: Project H into a single CRT block")
    print("=" * 80)
    print()

    # Sort blocks by size, take the largest
    blocks = sorted(joint_eigenspaces.items(), key=lambda x: -len(x[1]))

    # Test ALL blocks, not just one
    print(f"  Analyzing all {len(blocks)} CRT blocks:")
    print()
    print(f"  {'Block':<6} {'Dim':<5} {'IN-vals':<35} {'H_block spread':<16} {'V13 distinct':<14} {'H nontrivial?'}")
    print(f"  {'-'*6} {'-'*5} {'-'*35} {'-'*16} {'-'*14} {'-'*13}")

    best_block_key = None
    best_H_spread = 0
    best_V13_distinct = 0

    for b_idx, (sig, indices) in enumerate(blocks):
        dim = len(indices)

        # Build projection matrix: columns are position-basis vectors |r_i>
        P_block = np.zeros((phi, dim), dtype=complex)
        for j, idx in enumerate(indices):
            P_block[idx, j] = 1.0

        # Project H into block
        H_block = P_block.conj().T @ H @ P_block
        evals_block = np.linalg.eigvalsh(H_block)
        H_spread = evals_block.max() - evals_block.min()

        # Project V_13 into block
        V13_block = P_block.conj().T @ V_13 @ P_block
        evals_V13 = np.linalg.eigvalsh(V13_block)
        v13_unique = len(set(np.round(evals_V13, 6)))

        nontrivial = "YES" if H_spread > 0.01 else "no"

        sig_str = f"({sig[0]:+.3f},{sig[1]:+.3f},{sig[2]:+.3f})"
        print(f"  {b_idx:<6} {dim:<5} {sig_str:<35} {H_spread:<16.6f} {v13_unique:<14} {nontrivial}")

        if H_spread > best_H_spread and v13_unique > 1:
            best_H_spread = H_spread
            best_V13_distinct = v13_unique
            best_block_key = sig

    print()
    print(f"  Best block: H_spread = {best_H_spread:.6f}, V_13 distinct = {best_V13_distinct}")
    print()

    # =========================================================================
    # STEP 3: Deep analysis of best block
    # =========================================================================

    print("=" * 80)
    print("  STEP 3: Deep analysis of best CRT block")
    print("=" * 80)
    print()

    indices = joint_eigenspaces[best_block_key]
    dim = len(indices)
    residue_vals = [residues[i] for i in indices]

    print(f"  Block signature: (cos5, cos7, cos11) = {best_block_key}")
    print(f"  Dimension: {dim}")
    print(f"  Residues: {residue_vals}")
    print()

    # Build projection
    P_block = np.zeros((phi, dim), dtype=complex)
    for j, idx in enumerate(indices):
        P_block[idx, j] = 1.0

    # Project everything into block
    H_block = P_block.conj().T @ H @ P_block
    V5_block = P_block.conj().T @ V_5 @ P_block
    V7_block = P_block.conj().T @ V_7 @ P_block
    V11_block = P_block.conj().T @ V_11 @ P_block
    V13_block = P_block.conj().T @ V_13 @ P_block
    V17_block = P_block.conj().T @ V_17 @ P_block

    # Verify IN-primes are proportional to identity
    print(f"  IN-prime projections (should be proportional to I):")
    for p, V_block in [(5, V5_block), (7, V7_block), (11, V11_block)]:
        evals = np.linalg.eigvalsh(V_block)
        spread = evals.max() - evals.min()
        mean_val = np.mean(evals)
        print(f"    V_{p:2d}: mean = {mean_val:+.6f}, spread = {spread:.2e}")

    print()

    # V_13 in block
    evals_V13, evecs_V13 = np.linalg.eigh(V13_block)
    print(f"  V_13 in block:")
    print(f"    Eigenvalues: {np.round(evals_V13, 4)}")
    print(f"    Spread: {evals_V13.max() - evals_V13.min():.4f}")
    print()

    # H in block
    evals_H_block, evecs_H_block = np.linalg.eigh(H_block)
    print(f"  H in block:")
    print(f"    Eigenvalues: {np.round(evals_H_block, 6)}")
    print(f"    Spread: {evals_H_block.max() - evals_H_block.min():.6f}")
    print()

    # Key question: [H_block, V13_block] = ?
    comm = H_block @ V13_block - V13_block @ H_block
    comm_norm = np.linalg.norm(comm)
    print(f"  ||[H_block, V_13_block]|| = {comm_norm:.6f}")
    print()

    # =========================================================================
    # STEP 4: Construct qubit from H_block eigenstates
    # =========================================================================

    print("=" * 80)
    print("  STEP 4: Qubit from H_block eigenstates")
    print("=" * 80)
    print()

    # Find H_block eigenstates with DIFFERENT V_13 expectation values
    v13_expectations = np.zeros(dim)
    for i in range(dim):
        psi = evecs_H_block[:, i]
        v13_expectations[i] = np.real(np.vdot(psi, V13_block @ psi))

    print(f"  H_block eigenstate V_13 expectations:")
    for i in range(dim):
        print(f"    |{i}⟩: λ_H = {evals_H_block[i]:+.6f}, ⟨V_13⟩ = {v13_expectations[i]:+.6f}")
    print()

    # Find pair with maximal V_13 splitting
    best_pair = None
    best_split = 0
    for i in range(dim):
        for j in range(i+1, dim):
            split = abs(v13_expectations[i] - v13_expectations[j])
            if split > best_split:
                best_split = split
                best_pair = (i, j)

    i0, i1 = best_pair
    print(f"  Best qubit pair: |{i0}⟩ and |{i1}⟩")
    print(f"    V_13 split: {best_split:.6f}")
    print(f"    H gap: {abs(evals_H_block[i0] - evals_H_block[i1]):.6f}")
    print()

    # =========================================================================
    # STEP 5: Evolution within the CRT block
    # =========================================================================

    print("=" * 80)
    print("  STEP 5: Evolution within the CRT block")
    print("=" * 80)
    print()

    # Qubit states in block coordinates
    psi_0_block = evecs_H_block[:, i0]
    psi_1_block = evecs_H_block[:, i1]
    psi_cav_block = (psi_0_block + psi_1_block) / sqrt(2)
    psi_cav_block = psi_cav_block / np.linalg.norm(psi_cav_block)

    # Lift to full space
    psi_0_full = P_block @ psi_0_block
    psi_1_full = P_block @ psi_1_block
    psi_cav_full = P_block @ psi_cav_block

    t_max = 100.0
    n_points = 2000
    times = np.linspace(0, t_max, n_points)
    gate_amplitude = 0.5

    print(f"  === BLOCK-INTERNAL evolution (16D projected dynamics) ===")
    print()

    # Evolve within the block using H_block + gate * V13_block
    for label, H_evo in [("H_block only", H_block),
                          ("H_block + 0.5*V13", H_block + gate_amplitude * V13_block),
                          ("V13_block only", gate_amplitude * V13_block)]:
        evals_evo, evecs_evo = np.linalg.eigh(H_evo)
        psi_evo = evecs_evo.conj().T @ psi_cav_block

        F = np.zeros(n_points)
        F0 = np.zeros(n_points)
        F1 = np.zeros(n_points)
        for k, t in enumerate(times):
            phases = np.exp(-1j * evals_evo * t)
            psi_t = evecs_evo @ (phases * psi_evo)
            F[k] = np.abs(np.vdot(psi_cav_block, psi_t))**2
            F0[k] = np.abs(np.vdot(psi_0_block, psi_t))**2
            F1[k] = np.abs(np.vdot(psi_1_block, psi_t))**2

        swing = F.max() - F.min()
        print(f"  {label}:")
        print(f"    Cavity fidelity swing: {swing:.4f}")
        print(f"    F0 range: [{F0.min():.4f}, {F0.max():.4f}]")
        print(f"    F1 range: [{F1.min():.4f}, {F1.max():.4f}]")
        print()

    # =========================================================================
    # STEP 6: FULL-SPACE evolution with block-projected initial state
    # =========================================================================

    print("=" * 80)
    print("  STEP 6: Full-space evolution (leakage test)")
    print("=" * 80)
    print()

    # Start from psi_cav_full, evolve under FULL H + V_13
    H_gate_full = H + gate_amplitude * V_13
    evals_gf, evecs_gf = np.linalg.eigh(H_gate_full)
    psi_gf = evecs_gf.conj().T @ psi_cav_full

    F_full = np.zeros(n_points)
    F0_full = np.zeros(n_points)
    F1_full = np.zeros(n_points)
    block_content = np.zeros(n_points)

    # P_block projector
    P_proj = P_block @ P_block.conj().T

    for k, t in enumerate(times):
        phases = np.exp(-1j * evals_gf * t)
        psi_t = evecs_gf @ (phases * psi_gf)
        F_full[k] = np.abs(np.vdot(psi_cav_full, psi_t))**2
        F0_full[k] = np.abs(np.vdot(psi_0_full, psi_t))**2
        F1_full[k] = np.abs(np.vdot(psi_1_full, psi_t))**2
        block_content[k] = np.real(np.vdot(psi_t, P_proj @ psi_t))

    swing_full = F_full.max() - F_full.min()

    print(f"  Full H + 0.5*V_13 evolution:")
    print(f"    Cavity fidelity swing: {swing_full:.4f}")
    print(f"    F0 range: [{F0_full.min():.4f}, {F0_full.max():.4f}]")
    print(f"    F1 range: [{F1_full.min():.4f}, {F1_full.max():.4f}]")
    print(f"    Block content range: [{block_content.min():.4f}, {block_content.max():.4f}]")
    print(f"    Min block content: {block_content.min():.4f} (leakage = {1-block_content.min():.4f})")
    print()

    # =========================================================================
    # STEP 7: THE CRITICAL TEST — Noise immunity within block
    # =========================================================================

    print("=" * 80)
    print("  STEP 7: Noise immunity — clean vs noisy V_13 gate")
    print("=" * 80)
    print()

    noise_amp = 0.1

    # Test A: Block-internal (projected dynamics, no leakage)
    print(f"  === A: Block-internal dynamics (16D) ===")
    H_clean_block = H_block + gate_amplitude * V13_block
    H_noisy_block = H_block + gate_amplitude * V13_block + noise_amp * (V5_block + V7_block)

    ec, vc = np.linalg.eigh(H_clean_block)
    en, vn = np.linalg.eigh(H_noisy_block)
    pc = vc.conj().T @ psi_cav_block
    pn = vn.conj().T @ psi_cav_block

    Fc = np.zeros(n_points)
    Fn = np.zeros(n_points)
    for k, t in enumerate(times):
        psi_c = vc @ (np.exp(-1j * ec * t) * pc)
        psi_n = vn @ (np.exp(-1j * en * t) * pn)
        Fc[k] = np.abs(np.vdot(psi_cav_block, psi_c))**2
        Fn[k] = np.abs(np.vdot(psi_cav_block, psi_n))**2

    max_diff_block = np.max(np.abs(Fc - Fn))
    avg_diff_block = np.mean(np.abs(Fc - Fn))

    print(f"    Clean swing: {Fc.max()-Fc.min():.4f}")
    print(f"    Noisy swing: {Fn.max()-Fn.min():.4f}")
    print(f"    Max |F_clean - F_noisy|: {max_diff_block:.6f}")
    print(f"    Avg |F_clean - F_noisy|: {avg_diff_block:.6f}")
    print()

    # WHY this should be small: V5_block and V7_block are proportional to identity!
    # So the noise term is just a constant shift and doesn't affect dynamics
    print(f"    V5_block proportional to I? spread = {np.ptp(np.linalg.eigvalsh(V5_block)):.2e}")
    print(f"    V7_block proportional to I? spread = {np.ptp(np.linalg.eigvalsh(V7_block)):.2e}")
    print()

    # Test B: Full-space dynamics
    print(f"  === B: Full-space dynamics (480D) ===")
    H_clean_full = H + gate_amplitude * V_13
    H_noisy_full = H + gate_amplitude * V_13 + noise_amp * (V_5 + V_7)

    ecf, vcf = np.linalg.eigh(H_clean_full)
    enf, vnf = np.linalg.eigh(H_noisy_full)
    pcf = vcf.conj().T @ psi_cav_full
    pnf = vnf.conj().T @ psi_cav_full

    Fcf = np.zeros(n_points)
    Fnf = np.zeros(n_points)
    for k, t in enumerate(times):
        psi_c = vcf @ (np.exp(-1j * ecf * t) * pcf)
        psi_n = vnf @ (np.exp(-1j * enf * t) * pnf)
        Fcf[k] = np.abs(np.vdot(psi_cav_full, psi_c))**2
        Fnf[k] = np.abs(np.vdot(psi_cav_full, psi_n))**2

    max_diff_full = np.max(np.abs(Fcf - Fnf))
    avg_diff_full = np.mean(np.abs(Fcf - Fnf))

    print(f"    Clean swing: {Fcf.max()-Fcf.min():.4f}")
    print(f"    Noisy swing: {Fnf.max()-Fnf.min():.4f}")
    print(f"    Max |F_clean - F_noisy|: {max_diff_full:.6f}")
    print(f"    Avg |F_clean - F_noisy|: {avg_diff_full:.6f}")
    print()

    # =========================================================================
    # STEP 8: Rabi |0⟩ → |1⟩ within block
    # =========================================================================

    print("=" * 80)
    print("  STEP 8: Rabi |0⟩ → |1⟩ within the CRT block")
    print("=" * 80)
    print()

    # Start from |0⟩, evolve under V_13_block only (pure gate)
    psi_init = psi_0_block.copy()

    for label, H_evo in [("V13_block only", gate_amplitude * V13_block),
                          ("H_block + V13_block", H_block + gate_amplitude * V13_block)]:
        evals_evo, evecs_evo = np.linalg.eigh(H_evo)
        psi_evo = evecs_evo.conj().T @ psi_init

        F0 = np.zeros(n_points)
        F1 = np.zeros(n_points)
        for k, t in enumerate(times):
            phases = np.exp(-1j * evals_evo * t)
            psi_t = evecs_evo @ (phases * psi_evo)
            F0[k] = np.abs(np.vdot(psi_0_block, psi_t))**2
            F1[k] = np.abs(np.vdot(psi_1_block, psi_t))**2

        max_transfer = F1.max()
        t_transfer = times[np.argmax(F1)]
        F0_at_transfer = F0[np.argmax(F1)]
        leakage = 1 - F0_at_transfer - max_transfer

        print(f"  {label}:")
        print(f"    Max |1⟩ population: {max_transfer:.4f} at t = {t_transfer:.2f}")
        print(f"    |0⟩ at that time: {F0_at_transfer:.4f}")
        print(f"    Leakage: {leakage:.4f}")
        print()

    # =========================================================================
    # STEP 9: Scan ALL blocks for best protected qubit
    # =========================================================================

    print("=" * 80)
    print("  STEP 9: Scan all blocks for best noise immunity")
    print("=" * 80)
    print()

    print(f"  {'Block':<6} {'Dim':<5} {'H spread':<12} {'V13 split':<12} {'Block noise':<14} {'Full noise':<14}")
    print(f"  {'-'*6} {'-'*5} {'-'*12} {'-'*12} {'-'*14} {'-'*14}")

    best_immunity = 999
    best_block_idx = -1

    for b_idx, (sig, block_indices) in enumerate(blocks[:10]):  # Top 10
        dim_b = len(block_indices)

        P_b = np.zeros((phi, dim_b), dtype=complex)
        for j, idx in enumerate(block_indices):
            P_b[idx, j] = 1.0

        H_b = P_b.conj().T @ H @ P_b
        V13_b = P_b.conj().T @ V_13 @ P_b
        V5_b = P_b.conj().T @ V_5 @ P_b
        V7_b = P_b.conj().T @ V_7 @ P_b

        evals_Hb, evecs_Hb = np.linalg.eigh(H_b)
        H_spread_b = evals_Hb.max() - evals_Hb.min()

        # V_13 expectations on H_block eigenstates
        v13_exp = np.array([np.real(np.vdot(evecs_Hb[:,i], V13_b @ evecs_Hb[:,i]))
                           for i in range(dim_b)])
        v13_split_b = v13_exp.max() - v13_exp.min()

        # Find best qubit pair
        bp = None
        bs = 0
        for i in range(dim_b):
            for j in range(i+1, dim_b):
                s = abs(v13_exp[i] - v13_exp[j])
                if s > bs:
                    bs = s
                    bp = (i, j)

        if bp is None:
            continue

        i0b, i1b = bp
        psi_0b = evecs_Hb[:, i0b]
        psi_1b = evecs_Hb[:, i1b]
        psi_cb = (psi_0b + psi_1b) / sqrt(2)
        psi_cb = psi_cb / np.linalg.norm(psi_cb)

        # Block noise test
        H_c = H_b + gate_amplitude * V13_b
        H_n = H_b + gate_amplitude * V13_b + noise_amp * (V5_b + V7_b)

        ec_b, vc_b = np.linalg.eigh(H_c)
        en_b, vn_b = np.linalg.eigh(H_n)
        pc_b = vc_b.conj().T @ psi_cb
        pn_b = vn_b.conj().T @ psi_cb

        n_test = 500
        t_test = np.linspace(0, 50, n_test)
        Fc_b = np.zeros(n_test)
        Fn_b = np.zeros(n_test)
        for k, t in enumerate(t_test):
            psi_ct = vc_b @ (np.exp(-1j * ec_b * t) * pc_b)
            psi_nt = vn_b @ (np.exp(-1j * en_b * t) * pn_b)
            Fc_b[k] = np.abs(np.vdot(psi_cb, psi_ct))**2
            Fn_b[k] = np.abs(np.vdot(psi_cb, psi_nt))**2

        md_block = np.max(np.abs(Fc_b - Fn_b))

        # Full-space noise test
        psi_0f = P_b @ psi_0b
        psi_1f = P_b @ psi_1b
        psi_cf = P_b @ psi_cb

        ecf_b, vcf_b = np.linalg.eigh(H + gate_amplitude * V_13)
        enf_b, vnf_b = np.linalg.eigh(H + gate_amplitude * V_13 + noise_amp * (V_5 + V_7))
        pcf_b = vcf_b.conj().T @ psi_cf
        pnf_b = vnf_b.conj().T @ psi_cf

        Fcf_b = np.zeros(n_test)
        Fnf_b = np.zeros(n_test)
        for k, t in enumerate(t_test):
            psi_ct = vcf_b @ (np.exp(-1j * ecf_b * t) * pcf_b)
            psi_nt = vnf_b @ (np.exp(-1j * enf_b * t) * pnf_b)
            Fcf_b[k] = np.abs(np.vdot(psi_cf, psi_ct))**2
            Fnf_b[k] = np.abs(np.vdot(psi_cf, psi_nt))**2

        md_full = np.max(np.abs(Fcf_b - Fnf_b))

        print(f"  {b_idx:<6} {dim_b:<5} {H_spread_b:<12.6f} {bs:<12.6f} {md_block:<14.6f} {md_full:<14.6f}")

        if md_full < best_immunity:
            best_immunity = md_full
            best_block_idx = b_idx

    print()
    print(f"  Best full-space noise immunity: block {best_block_idx}, max_diff = {best_immunity:.6f}")
    print()

    # =========================================================================
    # FINAL VERDICT
    # =========================================================================

    print("=" * 80)
    print("  FINAL VERDICT")
    print("=" * 80)
    print()

    if max_diff_block < 1e-6:
        print("  BLOCK-INTERNAL RESULT:")
        print(f"    Within a CRT block, V_5 and V_7 are proportional to identity.")
        print(f"    Noise immunity is PERFECT by construction: max_diff = {max_diff_block:.2e}")
        print(f"    This is trivially true — the 'noise' is just a constant phase shift.")
        print()

    if max_diff_full < 0.05:
        print("  ╔════════════════════════════════════════════════════════════════════╗")
        print("  ║  ✓ THE THIRD OPTION: CRT-PROTECTED QUBIT IN FULL SPACE             ║")
        print("  ║  Block-projected states maintain protection under full H evolution  ║")
        print("  ╚════════════════════════════════════════════════════════════════════╝")
    elif max_diff_full < 0.15:
        print("  ╔════════════════════════════════════════════════════════════════════╗")
        print("  ║  ~ PARTIAL PROTECTION                                               ║")
        print(f"  ║  Full-space noise: {max_diff_full:.4f} (H leaks out of CRT block)    ║")
        print("  ╚════════════════════════════════════════════════════════════════════╝")
    else:
        print("  ╔════════════════════════════════════════════════════════════════════╗")
        print("  ║  THE WALL STANDS                                                    ║")
        print(f"  ║  Block-internal: {max_diff_block:.6f} (trivially perfect)              ║")
        print(f"  ║  Full-space:     {max_diff_full:.4f} (H scrambles out of block)        ║")
        print("  ║                                                                    ║")
        print("  ║  CRT protection is exact in position basis but H doesn't stay.     ║")
        print("  ║  The 6:1 selectivity ratio IS the real physics.                     ║")
        print("  ╚════════════════════════════════════════════════════════════════════╝")

    print()

    # =========================================================================
    # STEP 10: The Leakage Rate — How fast does H scramble the block?
    # =========================================================================

    print("=" * 80)
    print("  STEP 10: H leakage rate from CRT block")
    print("=" * 80)
    print()

    # Start in block, evolve under H only (no gate), measure block content vs time
    psi_init_full = psi_cav_full / np.linalg.norm(psi_cav_full)

    evals_H_only, evecs_H_only = np.linalg.eigh(H)
    psi_H = evecs_H_only.conj().T @ psi_init_full

    P_proj = P_block @ P_block.conj().T

    sample_times = np.linspace(0, 200, 4000)
    bc = np.zeros(len(sample_times))

    for k, t in enumerate(sample_times):
        phases = np.exp(-1j * evals_H_only * t)
        psi_t = evecs_H_only @ (phases * psi_H)
        bc[k] = np.real(np.vdot(psi_t, P_proj @ psi_t))

    print(f"  Block content under H-only evolution:")
    print(f"    t=0:   {bc[0]:.4f}")
    print(f"    t=10:  {bc[np.argmin(np.abs(sample_times-10))]:.4f}")
    print(f"    t=25:  {bc[np.argmin(np.abs(sample_times-25))]:.4f}")
    print(f"    t=50:  {bc[np.argmin(np.abs(sample_times-50))]:.4f}")
    print(f"    t=100: {bc[np.argmin(np.abs(sample_times-100))]:.4f}")
    print(f"    t=200: {bc[np.argmin(np.abs(sample_times-200))]:.4f}")
    print(f"    Min:   {bc.min():.4f} at t = {sample_times[np.argmin(bc)]:.1f}")
    print()

    # Half-life: when does block content first drop below 0.5?
    below_half = np.where(bc < 0.5)[0]
    if len(below_half) > 0:
        t_half = sample_times[below_half[0]]
        print(f"  Block half-life: t = {t_half:.1f}")
    else:
        print(f"  Block content never drops below 0.5! (min = {bc.min():.4f})")

    print()
    print("=" * 80)

if __name__ == "__main__":
    main()
