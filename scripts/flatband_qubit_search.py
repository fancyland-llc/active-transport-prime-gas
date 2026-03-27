#!/usr/bin/env python3
"""
FLAT-BAND QUBIT SEARCH
=======================
The cavity state (Perron/Anti-Perron) Rabi was free precession, not a gate.
The screening lives in the FLAT BAND, not the edge modes.

Question: Is there a qubit INSIDE the flat band where:
  - IN-primes (2,3,5,7,11) have NO eigenvalue splitting → screened
  - OUT-primes (13,17,19,...) HAVE eigenvalue splitting → addressable

If yes → quantum computer lives inside the flat band
If no  → the screening is real but not computationally useful

DIAGNOSTIC:
  1. Extract flat-band projector (465D subspace near λ ≈ 0)
  2. Project every V_p into flat-band subspace
  3. Compute eigenvalue spread of each projected V_p
  4. Compare IN-prime vs OUT-prime internal structure
  5. Find optimal qubit pair (max OUT-prime splitting, min IN-prime splitting)

Author: Claude + Tony
Date: 2026-03-25
"""

import numpy as np
from math import gcd, sqrt
import time

# =============================================================================
# INFRASTRUCTURE (same as before)
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

def build_V_p(residues, p):
    """V_p = diag(cos(2πr/p)), Frobenius-normalized."""
    diag = np.array([np.cos(2 * np.pi * r / p) for r in residues])
    V = np.diag(diag)
    norm = np.linalg.norm(V, 'fro')
    if norm > 0:
        V = V / norm
    return V

def build_V_p_raw(residues, p):
    """V_p = diag(cos(2πr/p)), NOT normalized."""
    diag = np.array([np.cos(2 * np.pi * r / p) for r in residues])
    return np.diag(diag)

def get_perron_eigenvalue(D_sym, tol=1e-10, max_iter=1000):
    n = D_sym.shape[0]
    v = np.random.RandomState(42).rand(n)
    v = v / np.linalg.norm(v)
    lam = 0
    for _ in range(max_iter):
        w = D_sym @ v
        lam_new = np.dot(v, w)
        v = w / np.linalg.norm(w)
        if abs(lam_new - lam) < tol:
            break
        lam = lam_new
    return lam_new

# =============================================================================
# BUILD SYSTEM
# =============================================================================

print("=" * 72)
print("FLAT-BAND QUBIT SEARCH")
print("Can we find a screened qubit inside the 465D flat band?")
print("=" * 72)

m = 2310
alpha_c = sqrt(135/88)
t0 = time.time()

residues = get_coprime_residues(m)
n = len(residues)
print(f"\nm = {m}, φ = {n}, α_c = {alpha_c:.6f}")

print("Building system...", end=" ", flush=True)
D_sym = build_D_sym(residues, m)
P_tau = build_P_tau(residues, m)
lambda_P = get_perron_eigenvalue(D_sym)
H = build_H_alpha(D_sym, P_tau, alpha_c, lambda_P)
evals_H, evecs_H = np.linalg.eigh(H)
print("done.")

# =============================================================================
# SECTION 1: IDENTIFY THE FLAT BAND
# =============================================================================

print("\n" + "=" * 72)
print("SECTION 1: FLAT BAND IDENTIFICATION")
print("=" * 72)

# Count eigenvalues near zero
thresholds = [0.001, 0.005, 0.01, 0.05, 0.1]
for thr in thresholds:
    count = np.sum(np.abs(evals_H) < thr)
    print(f"  |λ| < {thr}: {count}/{n} eigenvalues ({100*count/n:.1f}%)")

# Use threshold that captures the flat band cleanly
flat_threshold = 0.01
flat_mask = np.abs(evals_H) < flat_threshold
n_flat = np.sum(flat_mask)
flat_indices = np.where(flat_mask)[0]

# Also identify the edge modes
idx_P = np.argmax(evals_H)
idx_AP = np.argmin(evals_H)

print(f"\nFlat band: {n_flat} states with |λ| < {flat_threshold}")
print(f"Edge modes: Perron (λ={evals_H[idx_P]:.6f}), Anti-Perron (λ={evals_H[idx_AP]:.6f})")
print(f"Gap modes: {n - n_flat - 2} states between flat band and edges")

# Extract flat-band subspace
P_flat = evecs_H[:, flat_indices]  # n × n_flat matrix (projector basis)
print(f"Flat-band projector: {n} × {n_flat}")

# =============================================================================
# SECTION 2: PROJECT V_p INTO FLAT BAND
# =============================================================================

print("\n" + "=" * 72)
print("SECTION 2: V_p PROJECTED INTO FLAT BAND")
print("Comparing eigenvalue spread: IN-primes vs OUT-primes")
print("=" * 72)

in_primes = [2, 3, 5, 7, 11]
out_primes = [13, 17, 19, 23, 29, 31, 37, 41, 43]
all_test_primes = in_primes + out_primes

print(f"\n{'Prime':>5} | {'Type':>4} | {'λ_max':>10} | {'λ_min':>10} | {'Spread':>10} | {'Std':>10} | {'IPR':>10}")
print("-" * 80)

spreads_in = []
spreads_out = []
flat_evals_dict = {}

for p in all_test_primes:
    ptype = "IN" if p in in_primes else "OUT"
    
    # Build V_p (both normalized and raw)
    V_p = build_V_p(residues, p)
    
    # Project into flat band: V_flat = P_flat^† V_p P_flat
    V_flat = P_flat.conj().T @ V_p @ P_flat  # n_flat × n_flat
    
    # Diagonalize within flat band
    evals_flat, evecs_flat = np.linalg.eigh(V_flat.real)  # Should be real for diagonal V
    
    spread = np.max(evals_flat) - np.min(evals_flat)
    std = np.std(evals_flat)
    
    # IPR: how concentrated is the eigenvalue distribution?
    # If all eigenvalues are identical, IPR → 1/n_flat
    # If spread is broad, IPR is larger
    p_vals = np.abs(evals_flat)**2
    if np.sum(p_vals) > 0:
        p_norm = p_vals / np.sum(p_vals)
        ipr = np.sum(p_norm**2)
    else:
        ipr = 0
    
    flat_evals_dict[p] = evals_flat
    
    if ptype == "IN":
        spreads_in.append(spread)
    else:
        spreads_out.append(spread)
    
    print(f"  V_{p:>2} | {ptype:>4} | {np.max(evals_flat):>+10.6f} | {np.min(evals_flat):>+10.6f} | {spread:>10.6f} | {std:>10.6f} | {ipr:>10.6f}")

print(f"\n--- Summary ---")
print(f"  IN-prime  mean spread: {np.mean(spreads_in):.6f} ± {np.std(spreads_in):.6f}")
print(f"  OUT-prime mean spread: {np.mean(spreads_out):.6f} ± {np.std(spreads_out):.6f}")
ratio = np.mean(spreads_out) / np.mean(spreads_in) if np.mean(spreads_in) > 0 else float('inf')
print(f"  Ratio OUT/IN: {ratio:.4f}")

if ratio > 5:
    print(f"  >>> STRONG SELECTIVITY: OUT-primes have {ratio:.1f}× more structure in flat band")
elif ratio > 2:
    print(f"  >>> MODERATE SELECTIVITY: OUT-primes have {ratio:.1f}× more structure")
elif ratio > 1.2:
    print(f"  >>> WEAK SELECTIVITY: OUT-primes have slightly more structure ({ratio:.2f}×)")
else:
    print(f"  >>> NO SELECTIVITY: IN and OUT primes have similar flat-band structure")

# =============================================================================
# SECTION 3: EIGENVALUE HISTOGRAMS (text-based)
# =============================================================================

print("\n" + "=" * 72)
print("SECTION 3: FLAT-BAND EIGENVALUE DISTRIBUTIONS")
print("=" * 72)

for p in [5, 7, 13, 17]:
    evals = flat_evals_dict[p]
    ptype = "IN" if p in in_primes else "OUT"
    print(f"\n  V_{p} ({ptype}) — {n_flat} eigenvalues in flat band:")
    
    # Simple histogram with 20 bins
    n_bins = 20
    hist, bin_edges = np.histogram(evals, bins=n_bins)
    max_count = max(hist) if max(hist) > 0 else 1
    
    for i in range(n_bins):
        bar = "█" * int(40 * hist[i] / max_count)
        lo, hi = bin_edges[i], bin_edges[i+1]
        print(f"    [{lo:+.4f},{hi:+.4f}] {hist[i]:>4d} {bar}")

# =============================================================================
# SECTION 4: PAIRWISE COMMUTATORS WITHIN FLAT BAND
# =============================================================================

print("\n" + "=" * 72)
print("SECTION 4: COMMUTATORS [V_p, V_q] WITHIN FLAT BAND")
print("Do OUT-primes commute with each other? With IN-primes?")
print("=" * 72)

# Project all V_p into flat band
V_flat_dict = {}
for p in all_test_primes:
    V_p = build_V_p(residues, p)
    V_flat_dict[p] = P_flat.conj().T @ V_p @ P_flat

# Compute pairwise commutator norms
test_pairs = [
    (5, 7, "IN×IN"),
    (5, 11, "IN×IN"),
    (7, 11, "IN×IN"),
    (13, 17, "OUT×OUT"),
    (13, 19, "OUT×OUT"),
    (17, 19, "OUT×OUT"),
    (13, 23, "OUT×OUT"),
    (5, 13, "IN×OUT"),
    (7, 13, "IN×OUT"),
    (11, 13, "IN×OUT"),
    (5, 17, "IN×OUT"),
    (7, 17, "IN×OUT"),
]

print(f"\n{'Pair':>12} | {'Type':>8} | {'||[V_p,V_q]||':>14} | {'Interpretation'}")
print("-" * 65)

for p, q, ptype in test_pairs:
    Vp = V_flat_dict[p]
    Vq = V_flat_dict[q]
    comm = Vp @ Vq - Vq @ Vp
    cn = np.linalg.norm(comm, 'fro')
    
    if cn < 1e-10:
        interp = "COMMUTE (same axis)"
    elif cn < 1e-6:
        interp = "Barely non-commuting"
    elif cn < 1e-3:
        interp = "Weakly non-commuting"
    else:
        interp = "NON-COMMUTING (different axes!)"
    
    print(f"  [V_{p},V_{q}] | {ptype:>8} | {cn:>14.8f} | {interp}")

# =============================================================================
# SECTION 5: SEARCH FOR OPTIMAL QUBIT PAIR
# =============================================================================

print("\n" + "=" * 72)
print("SECTION 5: OPTIMAL QUBIT SEARCH")
print("Find two flat-band states with max V_13 splitting + min V_5 splitting")
print("=" * 72)

# Diagonalize V_13 restricted to flat band
V13_flat = V_flat_dict[13]
evals_13, evecs_13 = np.linalg.eigh(V13_flat.real)

# The best qubit from V_13's perspective: the two states with largest eigenvalue gap
idx_max = np.argmax(evals_13)
idx_min = np.argmin(evals_13)
gap_13 = evals_13[idx_max] - evals_13[idx_min]

print(f"\nV_13 in flat band:")
print(f"  Max eigenvalue: {evals_13[idx_max]:+.8f} (state #{idx_max})")
print(f"  Min eigenvalue: {evals_13[idx_min]:+.8f} (state #{idx_min})")
print(f"  Gap: {gap_13:.8f}")

# Now check: for these two states, what does V_5 look like?
psi_0 = evecs_13[:, idx_max]  # Logical |0⟩ candidate
psi_1 = evecs_13[:, idx_min]  # Logical |1⟩ candidate

# Project V_5 into this 2D subspace
qubit_basis = np.column_stack([psi_0, psi_1])
for p in all_test_primes:
    Vp_2x2 = qubit_basis.conj().T @ V_flat_dict[p] @ qubit_basis
    diag_split = abs(Vp_2x2[0,0] - Vp_2x2[1,1])
    off_diag = abs(Vp_2x2[0,1])
    ptype = "IN" if p in in_primes else "OUT"
    print(f"  V_{p:>2} ({ptype:>3}) in V_13-qubit: diag_split={diag_split:.8f}, off_diag={off_diag:.8f}")

# =============================================================================
# SECTION 6: THE REAL TEST — TIME EVOLUTION IN FLAT BAND
# =============================================================================

print("\n" + "=" * 72)
print("SECTION 6: TIME EVOLUTION WITHIN FLAT BAND")
print("Does V_13 rotate the flat-band qubit? Does V_5 leave it alone?")
print("=" * 72)

# Construct full-space versions of the flat-band qubit states
psi_0_full = P_flat @ psi_0  # Back to 480D
psi_1_full = P_flat @ psi_1
psi_fb_cavity = (psi_0_full + psi_1_full) / np.sqrt(2)
psi_fb_cavity = psi_fb_cavity / np.linalg.norm(psi_fb_cavity)

print(f"\nFlat-band qubit states (in full 480D space):")
print(f"  |0_fb⟩ norm: {np.linalg.norm(psi_0_full):.6f}")
print(f"  |1_fb⟩ norm: {np.linalg.norm(psi_1_full):.6f}")
print(f"  Overlap ⟨0|1⟩: {abs(np.dot(psi_0_full.conj(), psi_1_full)):.8f}")

# Time evolve under H + A*V_p for various primes
# Use eigendecomposition for efficiency
def measure_rabi_flatband(H_total, psi_init, psi_ref, t_max=200.0, dt=0.05):
    """Measure fidelity oscillation of psi_init under H_total, referenced against psi_ref."""
    evals, evecs = np.linalg.eigh(H_total)
    coeffs = evecs.conj().T @ psi_init
    
    ts = np.arange(0, t_max, dt)
    F = np.zeros(len(ts))
    for k, t in enumerate(ts):
        phases = np.exp(-1j * evals * t)
        psi_t = evecs @ (coeffs * phases)
        F[k] = abs(np.dot(psi_ref.conj(), psi_t))**2
    
    return ts, F

A_gate = 0.5

print(f"\nTime evolution of flat-band cavity state under H + {A_gate}·V_p:")
print(f"  (Measuring F(t) = |⟨ψ_fb_cav|ψ(t)⟩|²)")
print(f"\n{'Prime':>5} | {'Type':>4} | {'F_min':>8} | {'F_max':>8} | {'Swing':>8} | {'Status'}")
print("-" * 60)

for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]:
    V_p = build_V_p(residues, p)
    H_total = H + A_gate * V_p
    
    ts, F = measure_rabi_flatband(H_total, psi_fb_cavity, psi_fb_cavity, t_max=200.0, dt=0.05)
    F_min, F_max = np.min(F), np.max(F)
    swing = F_max - F_min
    ptype = "IN" if p in in_primes else "OUT"
    
    if swing < 0.05:
        status = "FROZEN ✓"
    elif swing > 0.5:
        status = "ROTATES ✗" if ptype == "IN" else "GATE ✓"
    else:
        status = "WOBBLE"
    
    print(f"  V_{p:>2} | {ptype:>4} |  {F_min:>6.4f} |  {F_max:>6.4f} |  {swing:>6.4f} | {status}")

# =============================================================================
# SECTION 7: RAW V_p TEST (unnormalized, to see if normalization matters)
# =============================================================================

print("\n" + "=" * 72)
print("SECTION 7: REPEAT WITH UNNORMALIZED V_p")
print("(Testing if Frobenius normalization is hiding the selectivity)")
print("=" * 72)

print(f"\n{'Prime':>5} | {'Type':>4} | {'||V||_F':>8} | {'F_min':>8} | {'F_max':>8} | {'Swing':>8} | {'Status'}")
print("-" * 72)

for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]:
    V_p_raw = build_V_p_raw(residues, p)
    frob = np.linalg.norm(V_p_raw, 'fro')
    H_total = H + A_gate * V_p_raw
    
    ts, F = measure_rabi_flatband(H_total, psi_fb_cavity, psi_fb_cavity, t_max=200.0, dt=0.05)
    F_min, F_max = np.min(F), np.max(F)
    swing = F_max - F_min
    ptype = "IN" if p in in_primes else "OUT"
    
    if swing < 0.05:
        status = "FROZEN ✓"
    elif swing > 0.5:
        status = "ROTATES ✗" if ptype == "IN" else "GATE ✓"
    else:
        status = "WOBBLE"
    
    print(f"  V_{p:>2} | {ptype:>4} | {frob:>8.2f} | {F_min:>6.4f} | {F_max:>6.4f} | {swing:>6.4f} | {status}")

# =============================================================================
# SECTION 8: THE DEFINITIVE IN vs OUT TEST
# =============================================================================

print("\n" + "=" * 72)
print("SECTION 8: DEFINITIVE SELECTIVITY TEST")
print("Evolve under ONLY V_p (no H base) in the flat band subspace")
print("This isolates the perturbation's effect from the base Hamiltonian")
print("=" * 72)

# The flat band has eigenvalues ≈ 0, so H ≈ 0 in this subspace.
# If we evolve under V_p alone (projected to flat band), we see ONLY the perturbation.

print(f"\nEvolution under V_p ONLY (no H), within flat band:")
print(f"Starting state: flat-band V_13-qubit cavity = (|0_fb⟩ + |1_fb⟩)/√2")

fb_cavity = (psi_0 + psi_1) / np.sqrt(2)  # In flat-band coordinates
fb_cavity = fb_cavity / np.linalg.norm(fb_cavity)

print(f"\n{'Prime':>5} | {'Type':>4} | {'F_min':>8} | {'F_max':>8} | {'Swing':>8} | {'Period':>8} | {'Status'}")
print("-" * 68)

for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]:
    Vp_fb = V_flat_dict[p].real  # n_flat × n_flat
    
    # Diagonalize V_p in flat band
    ev, ew = np.linalg.eigh(Vp_fb)
    coeffs = ew.conj().T @ fb_cavity
    
    # Time evolution under V_p alone at amplitude A
    A = 0.5
    t_max = 500.0
    dt = 0.1
    ts = np.arange(0, t_max, dt)
    F = np.zeros(len(ts))
    
    for k, t in enumerate(ts):
        phases = np.exp(-1j * A * ev * t)
        psi_t = ew @ (coeffs * phases)
        F[k] = abs(np.dot(fb_cavity.conj(), psi_t))**2
    
    F_min, F_max = np.min(F), np.max(F)
    swing = F_max - F_min
    ptype = "IN" if p in in_primes else "OUT"
    
    # Find period
    period = np.nan
    if swing > 0.1:
        for k in range(1, len(F) - 1):
            if F[k] < F[k-1] and F[k] < F[k+1] and F[k] < 0.5:
                period = 2 * ts[k]
                break
    
    if swing < 0.05:
        status = "SCREENED ✓"
    elif swing > 0.5:
        status = "ROTATES"
    else:
        status = f"WOBBLE ({swing:.3f})"
    
    period_str = f"{period:.1f}" if not np.isnan(period) else "N/A"
    print(f"  V_{p:>2} | {ptype:>4} | {F_min:>6.4f} | {F_max:>6.4f} | {swing:>6.4f} | {period_str:>8} | {status}")

# =============================================================================
# SECTION 9: THE RATIO THAT MATTERS
# =============================================================================

print("\n" + "=" * 72)
print("SECTION 9: V_p COUPLING STRENGTH IN FLAT-BAND QUBIT")
print("How strongly does each V_p couple |0_fb⟩ ↔ |1_fb⟩ ?")
print("=" * 72)

# The coupling strength is the off-diagonal element ⟨0_fb|V_p|1_fb⟩
# AND the detuning is the diagonal difference ⟨0|V|0⟩ - ⟨1|V|1⟩
# Rabi frequency = sqrt(detuning² + 4*coupling²)

print(f"\n{'Prime':>5} | {'Type':>4} | {'Coupling |g|':>12} | {'Detuning δ':>12} | {'Ω_Rabi':>12} | {'T_Rabi':>10}")
print("-" * 72)

for p in all_test_primes:
    Vp_fb = V_flat_dict[p]
    g = abs(psi_0.conj() @ Vp_fb @ psi_1)           # Off-diagonal coupling
    delta = abs(psi_0.conj() @ Vp_fb @ psi_0 - psi_1.conj() @ Vp_fb @ psi_1)  # Detuning
    omega = 0.5 * np.sqrt(delta**2 + 4*g**2)  # Rabi frequency at A=0.5
    T_rabi = 2*np.pi/omega if omega > 1e-12 else float('inf')
    ptype = "IN" if p in in_primes else "OUT"
    
    T_str = f"{T_rabi:.1f}" if T_rabi < 1e6 else "∞"
    print(f"  V_{p:>2} | {ptype:>4} | {g:>12.8f} | {delta:>12.8f} | {omega:>12.8f} | {T_str:>10}")


# =============================================================================
# FINAL ASSESSMENT
# =============================================================================

print("\n" + "=" * 72)
print("FINAL ASSESSMENT")
print("=" * 72)

# Compare average coupling for IN vs OUT
in_couplings = []
out_couplings = []
for p in all_test_primes:
    Vp_fb = V_flat_dict[p]
    g = abs(psi_0.conj() @ Vp_fb @ psi_1)
    delta = abs(psi_0.conj() @ Vp_fb @ psi_0 - psi_1.conj() @ Vp_fb @ psi_1)
    strength = np.sqrt(delta**2 + 4*g**2)
    if p in in_primes:
        in_couplings.append(strength)
    else:
        out_couplings.append(strength)

avg_in = np.mean(in_couplings)
avg_out = np.mean(out_couplings)
selectivity = avg_out / avg_in if avg_in > 0 else float('inf')

print(f"\n  Average IN-prime coupling in flat-band qubit:  {avg_in:.8f}")
print(f"  Average OUT-prime coupling in flat-band qubit: {avg_out:.8f}")
print(f"  Selectivity ratio (OUT/IN): {selectivity:.4f}")

if selectivity > 10:
    print(f"\n  >>> STRONG SELECTIVITY ({selectivity:.1f}×)")
    print(f"  >>> A screened qubit exists inside the flat band!")
    print(f"  >>> IN-primes are {selectivity:.1f}× weaker — arithmetic immunity preserved")
elif selectivity > 3:
    print(f"\n  >>> MODERATE SELECTIVITY ({selectivity:.1f}×)")
    print(f"  >>> Partial screening — may be useful with error correction")
elif selectivity > 1.5:
    print(f"\n  >>> WEAK SELECTIVITY ({selectivity:.1f}×)")
    print(f"  >>> Marginal advantage — probably not computationally useful")
else:
    print(f"\n  >>> NO SELECTIVITY ({selectivity:.2f}×)")
    print(f"  >>> The flat-band V_13 qubit has no arithmetic immunity")
    print(f"  >>> IN and OUT primes couple equally to these states")

print(f"\nDone in {time.time()-t0:.1f}s")
