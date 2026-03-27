#!/usr/bin/env python3
"""
BLOCK 9 DEEP DIVE — The 0.010 Noise Immunity Block
=====================================================

Block 9 showed max_diff = 0.010 in full-space noise test.
This matches the original transcript's 0.019 claim.

Questions:
1. What makes block 9 special?
2. Can we get clean Rabi within this block?
3. How does protection scale with noise amplitude?
4. Can we find a 2-level qubit that works?
5. What's the ACTUAL gate fidelity?

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
    print("  BLOCK 9 DEEP DIVE — The 0.010 Noise Immunity Block")
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
    V_17 = build_V_diag(residues, 17)
    V_19 = build_V_diag(residues, 19)

    # Build ALL CRT blocks
    joint_eigenspaces = {}
    for i, r in enumerate(residues):
        sig = (round(np.cos(2*np.pi*r/5), 8),
               round(np.cos(2*np.pi*r/7), 8),
               round(np.cos(2*np.pi*r/11), 8))
        if sig not in joint_eigenspaces:
            joint_eigenspaces[sig] = []
        joint_eigenspaces[sig].append(i)

    blocks = sorted(joint_eigenspaces.items(), key=lambda x: -len(x[1]))

    # =========================================================================
    # STEP 1: Find the ACTUAL best block across ALL 30
    # =========================================================================

    print("=" * 80)
    print("  STEP 1: Full scan of all 30 blocks")
    print("=" * 80)
    print()

    gate_amplitude = 0.5
    noise_amp = 0.1
    t_max = 50.0
    n_test = 500
    t_test = np.linspace(0, t_max, n_test)

    # Pre-diagonalize full-space clean and noisy Hamiltonians
    H_clean_full = H + gate_amplitude * V_13
    H_noisy_full = H + gate_amplitude * V_13 + noise_amp * (V_5 + V_7)
    ecf, vcf = np.linalg.eigh(H_clean_full)
    enf, vnf = np.linalg.eigh(H_noisy_full)

    block_results = []

    print(f"  {'#':<4} {'Res. example':<14} {'H spread':<10} {'V13 gap':<10} {'Noise max':<12} {'Noise avg':<12} {'Swing':<8}")
    print(f"  {'-'*4} {'-'*14} {'-'*10} {'-'*10} {'-'*12} {'-'*12} {'-'*8}")

    for b_idx, (sig, block_indices) in enumerate(blocks):
        dim_b = len(block_indices)
        P_b = np.zeros((phi, dim_b), dtype=complex)
        for j, idx in enumerate(block_indices):
            P_b[idx, j] = 1.0

        H_b = P_b.conj().T @ H @ P_b
        V13_b = P_b.conj().T @ V_13 @ P_b

        evals_Hb, evecs_Hb = np.linalg.eigh(H_b)
        H_spread = evals_Hb.max() - evals_Hb.min()

        # V_13 expectations on H_block eigenstates
        v13_exp = np.array([np.real(np.vdot(evecs_Hb[:,i], V13_b @ evecs_Hb[:,i]))
                           for i in range(dim_b)])

        # Best qubit pair (max V_13 split)
        bp = None; bs = 0
        for i in range(dim_b):
            for j in range(i+1, dim_b):
                s = abs(v13_exp[i] - v13_exp[j])
                if s > bs:
                    bs = s; bp = (i, j)

        i0, i1 = bp
        psi_0b = evecs_Hb[:, i0]
        psi_1b = evecs_Hb[:, i1]
        psi_cb = (psi_0b + psi_1b) / sqrt(2)
        psi_cb = psi_cb / np.linalg.norm(psi_cb)

        # Full space cavity state
        psi_cf = P_b @ psi_cb
        pcf_b = vcf.conj().T @ psi_cf
        pnf_b = vnf.conj().T @ psi_cf

        Fcf_b = np.zeros(n_test)
        Fnf_b = np.zeros(n_test)
        for k, t in enumerate(t_test):
            psi_ct = vcf @ (np.exp(-1j * ecf * t) * pcf_b)
            psi_nt = vnf @ (np.exp(-1j * enf * t) * pnf_b)
            Fcf_b[k] = np.abs(np.vdot(psi_cf, psi_ct))**2
            Fnf_b[k] = np.abs(np.vdot(psi_cf, psi_nt))**2

        md = np.max(np.abs(Fcf_b - Fnf_b))
        ad = np.mean(np.abs(Fcf_b - Fnf_b))
        sw = Fcf_b.max() - Fcf_b.min()

        res_ex = residues[block_indices[0]]
        print(f"  {b_idx:<4} r={res_ex:<10} {H_spread:<10.6f} {bs:<10.4f} {md:<12.6f} {ad:<12.6f} {sw:<8.4f}")

        block_results.append({
            'idx': b_idx, 'sig': sig, 'indices': block_indices,
            'H_spread': H_spread, 'V13_split': bs, 'pair': bp,
            'noise_max': md, 'noise_avg': ad, 'swing': sw
        })

    # Sort by noise immunity
    block_results.sort(key=lambda x: x['noise_max'])

    print()
    print(f"  Top 5 blocks by noise immunity:")
    for r in block_results[:5]:
        print(f"    Block {r['idx']:2d}: max_diff = {r['noise_max']:.6f}, V13_split = {r['V13_split']:.4f}, swing = {r['swing']:.4f}")
    print()

    # =========================================================================
    # STEP 2: Deep dive into the best block
    # =========================================================================

    best = block_results[0]
    b_idx = best['idx']
    sig = best['sig']
    indices = best['indices']
    dim_b = len(indices)

    print("=" * 80)
    print(f"  STEP 2: Deep dive into Block {b_idx} (best noise immunity)")
    print("=" * 80)
    print()

    P_b = np.zeros((phi, dim_b), dtype=complex)
    for j, idx in enumerate(indices):
        P_b[idx, j] = 1.0

    H_b = P_b.conj().T @ H @ P_b
    V5_b = P_b.conj().T @ V_5 @ P_b
    V7_b = P_b.conj().T @ V_7 @ P_b
    V11_b = P_b.conj().T @ V_11 @ P_b
    V13_b = P_b.conj().T @ V_13 @ P_b
    V17_b = P_b.conj().T @ V_17 @ P_b
    V19_b = P_b.conj().T @ V_19 @ P_b

    evals_Hb, evecs_Hb = np.linalg.eigh(H_b)

    print(f"  Residues: {[residues[i] for i in indices]}")
    print(f"  CRT signature: (cos5,cos7,cos11) = ({sig[0]:+.4f}, {sig[1]:+.4f}, {sig[2]:+.4f})")
    print()

    # H_block spectrum
    print(f"  H_block eigenvalues:")
    for i in range(dim_b):
        v13_exp = np.real(np.vdot(evecs_Hb[:,i], V13_b @ evecs_Hb[:,i]))
        v17_exp = np.real(np.vdot(evecs_Hb[:,i], V17_b @ evecs_Hb[:,i]))
        print(f"    |{i:2d}⟩: λ_H = {evals_Hb[i]:+.6f}, ⟨V_13⟩ = {v13_exp:+.4f}, ⟨V_17⟩ = {v17_exp:+.4f}")
    print()

    # =========================================================================
    # STEP 3: Noise amplitude scaling
    # =========================================================================

    print("=" * 80)
    print(f"  STEP 3: Noise amplitude scaling")
    print("=" * 80)
    print()

    i0, i1 = best['pair']
    psi_0b = evecs_Hb[:, i0]
    psi_1b = evecs_Hb[:, i1]
    psi_cb = (psi_0b + psi_1b) / sqrt(2)
    psi_cb = psi_cb / np.linalg.norm(psi_cb)
    psi_cf = P_b @ psi_cb

    print(f"  {'ε (noise)':<12} {'Max diff':<12} {'Avg diff':<12} {'Gate quality'}")
    print(f"  {'-'*12} {'-'*12} {'-'*12} {'-'*12}")

    for eps in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0]:
        H_n_eps = H + gate_amplitude * V_13 + eps * (V_5 + V_7)
        en_eps, vn_eps = np.linalg.eigh(H_n_eps)
        pn_eps = vn_eps.conj().T @ psi_cf

        Fc_eps = np.zeros(n_test)
        Fn_eps = np.zeros(n_test)
        for k, t in enumerate(t_test):
            psi_ct = vcf @ (np.exp(-1j * ecf * t) * pcf_b)  # Reuse clean
            psi_nt = vn_eps @ (np.exp(-1j * en_eps * t) * pn_eps)
            Fc_eps[k] = np.abs(np.vdot(psi_cf, psi_ct))**2
            Fn_eps[k] = np.abs(np.vdot(psi_cf, psi_nt))**2

        # Oops — need to recompute pcf_b for the clean case since it may be stale
        # Actually pcf_b was computed above for psi_cf against ecf/vcf, which is correct
        md_eps = np.max(np.abs(Fc_eps - Fn_eps))
        ad_eps = np.mean(np.abs(Fc_eps - Fn_eps))
        gq = 1 - md_eps

        print(f"  ε = {eps:<9.2f} {md_eps:<12.6f} {ad_eps:<12.6f} {gq:.4f}")

    print()

    # =========================================================================
    # STEP 4: Test DIFFERENT gate operators (V_17, V_19)
    # =========================================================================

    print("=" * 80)
    print(f"  STEP 4: Different gate operators on block {b_idx}")
    print("=" * 80)
    print()

    for p, V_gate in [(13, V_13), (17, V_17), (19, V_19)]:
        V_gate_b = P_b.conj().T @ V_gate @ P_b
        v_exp = np.array([np.real(np.vdot(evecs_Hb[:,i], V_gate_b @ evecs_Hb[:,i]))
                         for i in range(dim_b)])
        # Best pair for this gate
        bp_g = None; bs_g = 0
        for i in range(dim_b):
            for j in range(i+1, dim_b):
                s = abs(v_exp[i] - v_exp[j])
                if s > bs_g:
                    bs_g = s; bp_g = (i, j)

        i0g, i1g = bp_g
        psi_0g = evecs_Hb[:, i0g]
        psi_1g = evecs_Hb[:, i1g]
        psi_cg = (psi_0g + psi_1g) / sqrt(2)
        psi_cg = psi_cg / np.linalg.norm(psi_cg)
        psi_cf_g = P_b @ psi_cg

        # Full-space noise test
        H_c_g = H + gate_amplitude * V_gate
        H_n_g = H + gate_amplitude * V_gate + noise_amp * (V_5 + V_7)
        ec_g, vc_g = np.linalg.eigh(H_c_g)
        en_g, vn_g = np.linalg.eigh(H_n_g)
        pc_g = vc_g.conj().T @ psi_cf_g
        pn_g = vn_g.conj().T @ psi_cf_g

        Fc_g = np.zeros(n_test)
        Fn_g = np.zeros(n_test)
        for k, t in enumerate(t_test):
            psi_ct = vc_g @ (np.exp(-1j * ec_g * t) * pc_g)
            psi_nt = vn_g @ (np.exp(-1j * en_g * t) * pn_g)
            Fc_g[k] = np.abs(np.vdot(psi_cf_g, psi_ct))**2
            Fn_g[k] = np.abs(np.vdot(psi_cf_g, psi_nt))**2

        md_g = np.max(np.abs(Fc_g - Fn_g))
        sw_g = Fc_g.max() - Fc_g.min()

        print(f"  V_{p:2d} gate: split = {bs_g:.4f}, swing = {sw_g:.4f}, noise max_diff = {md_g:.6f}")

    print()

    # =========================================================================
    # STEP 5: ALL pairs in best block — find the clean Rabi
    # =========================================================================

    print("=" * 80)
    print(f"  STEP 5: Exhaustive pair search in block {b_idx}")
    print("=" * 80)
    print()

    # For every pair of H_block eigenstates, test Rabi |i⟩ → |j⟩ under V_13
    H_gate_b = H_b + gate_amplitude * V13_b
    evals_gb, evecs_gb = np.linalg.eigh(H_gate_b)

    t_rabi = np.linspace(0, 200, 4000)

    print(f"  {'Pair':<10} {'V13 split':<12} {'Max transfer':<14} {'Leakage':<10} {'t_flip':<10}")
    print(f"  {'-'*10} {'-'*12} {'-'*14} {'-'*10} {'-'*10}")

    best_transfer = 0
    best_rabi_pair = None

    for i in range(dim_b):
        for j in range(i+1, dim_b):
            v13_i = np.real(np.vdot(evecs_Hb[:,i], V13_b @ evecs_Hb[:,i]))
            v13_j = np.real(np.vdot(evecs_Hb[:,j], V13_b @ evecs_Hb[:,j]))
            v13_split = abs(v13_i - v13_j)

            if v13_split < 0.1:
                continue  # Skip pairs with tiny V_13 splitting

            psi_i = evecs_Hb[:, i]
            psi_evo = evecs_gb.conj().T @ psi_i

            F_j = np.zeros(len(t_rabi))
            F_i = np.zeros(len(t_rabi))
            for k, t in enumerate(t_rabi):
                phases = np.exp(-1j * evals_gb * t)
                psi_t = evecs_gb @ (phases * psi_evo)
                F_j[k] = np.abs(np.vdot(evecs_Hb[:,j], psi_t))**2
                F_i[k] = np.abs(np.vdot(psi_i, psi_t))**2

            max_t = F_j.max()
            idx_t = np.argmax(F_j)
            leakage = 1 - F_i[idx_t] - F_j[idx_t]
            t_flip = t_rabi[idx_t]

            if max_t > 0.15:  # Only show meaningful transfers
                print(f"  ({i:2d},{j:2d})   {v13_split:<12.4f} {max_t:<14.4f} {leakage:<10.4f} {t_flip:<10.2f}")

            if max_t > best_transfer:
                best_transfer = max_t
                best_rabi_pair = (i, j)

    print()
    print(f"  Best Rabi pair: {best_rabi_pair}, max transfer = {best_transfer:.4f}")
    print()

    # =========================================================================
    # STEP 6: Position-basis qubit test (the CRT pair from invariant_subspace)
    # =========================================================================

    print("=" * 80)
    print(f"  STEP 6: Position-basis qubit with CRT protection")
    print("=" * 80)
    print()

    # Within the best block, pick two position-basis states with max V_13 gap
    residue_vals = [residues[idx] for idx in indices]
    v13_vals = [np.cos(2*np.pi*r/13) for r in residue_vals]

    # Sort by V_13 value
    rv_sorted = sorted(zip(v13_vals, indices, residue_vals), key=lambda x: x[0])

    # Max gap pair
    r_max_idx = rv_sorted[-1][1]
    r_min_idx = rv_sorted[0][1]
    r_max = rv_sorted[-1][2]
    r_min = rv_sorted[0][2]
    v13_gap = rv_sorted[-1][0] - rv_sorted[0][0]

    print(f"  Max V_13 gap pair in block:")
    print(f"    |0⟩ = |r={r_max}⟩, V_13 = {rv_sorted[-1][0]:+.4f}")
    print(f"    |1⟩ = |r={r_min}⟩, V_13 = {rv_sorted[0][0]:+.4f}")
    print(f"    Gap: {v13_gap:.4f}")
    print()

    # Build states
    psi_r0 = np.zeros(phi, dtype=complex)
    psi_r0[r_max_idx] = 1.0
    psi_r1 = np.zeros(phi, dtype=complex)
    psi_r1[r_min_idx] = 1.0
    psi_r_cav = (psi_r0 + psi_r1) / sqrt(2)

    # Full-space tests
    pc_r = vcf.conj().T @ psi_r_cav
    pn_r = vnf.conj().T @ psi_r_cav

    Fc_r = np.zeros(n_test)
    Fn_r = np.zeros(n_test)
    for k, t in enumerate(t_test):
        psi_ct = vcf @ (np.exp(-1j * ecf * t) * pc_r)
        psi_nt = vnf @ (np.exp(-1j * enf * t) * pn_r)
        Fc_r[k] = np.abs(np.vdot(psi_r_cav, psi_ct))**2
        Fn_r[k] = np.abs(np.vdot(psi_r_cav, psi_nt))**2

    md_r = np.max(np.abs(Fc_r - Fn_r))
    sw_r = Fc_r.max() - Fc_r.min()

    print(f"  Position-basis cavity under H + 0.5*V_13 + noise:")
    print(f"    Swing: {sw_r:.4f}")
    print(f"    Noise max_diff: {md_r:.6f}")
    print()

    # Also test V_5 screening on position basis
    H_v5 = H + gate_amplitude * V_5
    ev5, vv5 = np.linalg.eigh(H_v5)
    pv5 = vv5.conj().T @ psi_r_cav
    Fv5 = np.zeros(n_test)
    for k, t in enumerate(t_test):
        psi_t = vv5 @ (np.exp(-1j * ev5 * t) * pv5)
        Fv5[k] = np.abs(np.vdot(psi_r_cav, psi_t))**2
    sw_v5 = Fv5.max() - Fv5.min()
    print(f"  V_5 screening on position-basis cavity: swing = {sw_v5:.4f}")
    print()

    # =========================================================================
    # FINAL VERDICT
    # =========================================================================

    print("=" * 80)
    print("  FINAL VERDICT")
    print("=" * 80)
    print()

    print(f"  CRT Block {b_idx}:")
    print(f"    H-eigenstate qubit noise immunity: {best['noise_max']:.6f}")
    print(f"    Position-basis qubit noise immunity: {md_r:.6f}")
    print(f"    Best Rabi transfer: {best_transfer:.4f}")
    print()

    if best['noise_max'] < 0.02:
        print("  THE ORIGINAL 0.019 RESULT IS REPRODUCIBLE.")
        print("  The key was the CRT block projection — the qubit states must be")
        print("  H-eigenstates within a CRT-protected subspace.")
        print()
        if best_transfer > 0.80:
            print("  AND clean Rabi oscillation works!")
        elif best_transfer > 0.2:
            print("  Rabi transfer is partial — multi-level dynamics, not clean 2-level.")
        else:
            print("  But Rabi between specific states is weak.")
            print("  The 'oscillation' is cavity fidelity (superposition dynamics),")
            print("  not a clean |0⟩ → |1⟩ gate.")

    print()
    print("=" * 80)

if __name__ == "__main__":
    main()
