#!/usr/bin/env python3
"""
micro_sweep.py — Ultra-high resolution probe of the crash point
================================================================

Tests whether the "plateau" between m=2310 and m=30030 is:
  (a) real convergence to α_c ≈ 1.2386, or
  (b) a grid-resolution artifact (both crash between same grid points)

Strategy:
  1. Micro-sweep 2000 points in [1.230, 1.245] for m=2310 (fast: ~5 min)
  2. Resolution: Δα = 0.0000075 — 300x finer than the original sweep
  3. Compare crash point to √(135/88) = 1.238584...
  4. If time/patience allows: repeat for m=30030 (slow: ~9 hours)
     (use --full flag for m=30030)

Usage:
  python micro_sweep.py               # m=2310 only (~5 min)
  python micro_sweep.py --both        # m=2310 + m=30030 (~9 hours)

Author: Antonio P. Matos / Fancyland LLC
Date: March 25, 2026
"""

import numpy as np
from math import gcd
import sys
import os
import time
import json

# ================================================================
SQRT_3_2    = np.sqrt(3.0 / 2.0)     # 1.224744871391589
SQRT_2_3    = np.sqrt(2.0 / 3.0)     # 0.816496580927726
SQRT_135_88 = np.sqrt(135.0 / 88.0)  # 1.238584235767...
# ================================================================


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


def gap_position(eigs):
    consec = np.diff(eigs)
    j_max = int(np.argmax(consec))
    return j_max / max(len(consec) - 1, 1)


def find_crash_alpha(alphas, positions, threshold=0.5):
    for i in range(len(positions) - 1):
        if positions[i] > threshold and positions[i + 1] <= threshold:
            frac = (threshold - positions[i]) / (positions[i + 1] - positions[i])
            return alphas[i] + frac * (alphas[i + 1] - alphas[i])
    return None


def sweep_modulus(m, alphas, label=""):
    res = coprime_residues(m)
    n = len(res)
    
    print(f"\n{'='*72}")
    print(f"  MICRO-SWEEP: m={m}, φ={n}  {label}")
    print(f"  α range: [{alphas[0]:.10f}, {alphas[-1]:.10f}]")
    print(f"  Points: {len(alphas)}, spacing: {alphas[1]-alphas[0]:.10f}")
    print(f"{'='*72}")
    
    t0 = time.time()
    print(f"  Building D_sym ({n}x{n})...", end=" ", flush=True)
    D = build_D_sym(res, m)
    print(f"{time.time()-t0:.1f}s")
    
    print(f"  Perron eigenvalue...", end=" ", flush=True)
    t0 = time.time()
    lam_P = perron_eigenvalue(D)
    print(f"λ_P = {lam_P:.6f} ({time.time()-t0:.1f}s)")
    
    perm = inversion_perm(res, m)
    C = commutator_fast(D, perm)
    A = D / lam_P
    B = C / lam_P
    del D, C, perm
    
    n_alpha = len(alphas)
    gap_pos_arr = np.empty(n_alpha)
    perron_gap_arr = np.empty(n_alpha)
    
    t0 = time.time()
    print(f"  Sweeping {n_alpha} α values...")
    
    for i, alpha in enumerate(alphas):
        H = A + (1j * alpha) * B
        eigs = np.linalg.eigvalsh(H)
        gap_pos_arr[i] = gap_position(eigs)
        perron_gap_arr[i] = eigs[-1] - eigs[-2]
        
        if n > 200 and ((i + 1) % 200 == 0 or i == 0):
            elapsed = time.time() - t0
            rate = elapsed / (i + 1)
            eta = rate * (n_alpha - i - 1)
            print(f"    [{i+1:>5d}/{n_alpha}] {elapsed:.1f}s elapsed, ~{eta:.0f}s ETA")
    
    elapsed = time.time() - t0
    print(f"  Done: {elapsed:.1f}s ({elapsed/n_alpha:.3f}s/α)")
    
    # Find crash at multiple thresholds for robustness
    crash_50 = find_crash_alpha(alphas, gap_pos_arr, 0.5)
    crash_25 = find_crash_alpha(alphas, gap_pos_arr, 0.25)
    crash_75 = find_crash_alpha(alphas, gap_pos_arr, 0.75)
    
    # Steepest descent
    derivs = np.diff(gap_pos_arr) / np.diff(alphas)
    j_min = int(np.argmin(derivs))
    steepest_alpha = 0.5 * (alphas[j_min] + alphas[j_min + 1])
    steepest_deriv = derivs[j_min]
    
    print(f"\n  RESULTS:")
    print(f"  {'Candidate':>20s}  {'Value':>14s}  {'Δ from crash':>14s}")
    print(f"  {'-'*20}  {'-'*14}  {'-'*14}")
    
    candidates = [
        ("√(3/2)", SQRT_3_2),
        ("√(135/88)", SQRT_135_88),
        ("Measured crash", crash_50),
        ("Steepest descent", steepest_alpha),
    ]
    
    for name, val in candidates:
        if val is not None and crash_50 is not None:
            delta = val - crash_50
            print(f"  {name:>20s}  {val:>14.10f}  {delta:>+14.10f}")
        elif val is not None:
            print(f"  {name:>20s}  {val:>14.10f}  {'N/A':>14s}")
    
    if crash_50:
        print(f"\n  Crash thresholds:")
        if crash_75: print(f"    gap_pos < 0.75 at α = {crash_75:.10f}")
        print(f"    gap_pos < 0.50 at α = {crash_50:.10f}")
        if crash_25: print(f"    gap_pos < 0.25 at α = {crash_25:.10f}")
        if crash_75 and crash_25:
            width = crash_25 - crash_75
            center = 0.5 * (crash_75 + crash_25)
            print(f"    Transition width (75→25): {width:.10f}")
            print(f"    Transition center: {center:.10f}")
        
        print(f"\n  Error analysis (crash @ 0.5):")
        err_32 = crash_50 - SQRT_3_2
        err_135 = crash_50 - SQRT_135_88
        print(f"    From √(3/2):    {err_32:+.12f} ({100*err_32/SQRT_3_2:+.8f}%)")
        print(f"    From √(135/88): {err_135:+.12f} ({100*err_135/SQRT_135_88:+.8f}%)")
    
    print(f"\n  Steepest descent:")
    print(f"    α = {steepest_alpha:.10f}")
    print(f"    d(gap_pos)/dα = {steepest_deriv:.2f}")
    err_steep_32 = steepest_alpha - SQRT_3_2
    err_steep_135 = steepest_alpha - SQRT_135_88
    print(f"    From √(3/2):    {err_steep_32:+.12f} ({100*err_steep_32/SQRT_3_2:+.8f}%)")
    print(f"    From √(135/88): {err_steep_135:+.12f} ({100*err_steep_135/SQRT_135_88:+.8f}%)")
    
    return {
        "m": m, "phi": n,
        "crash_50": crash_50, "crash_25": crash_25, "crash_75": crash_75,
        "steepest_alpha": steepest_alpha, "steepest_deriv": steepest_deriv,
        "alphas": alphas.tolist(),
        "gap_pos": gap_pos_arr.tolist(),
    }


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Ultra-fine grid around the transition
    # Center on √(135/88) = 1.238584235767...
    # Range: ±0.0075 → [1.231, 1.246]
    center = SQRT_135_88
    half_width = 0.0075
    n_points = 2000
    alphas = np.linspace(center - half_width, center + half_width, n_points)
    
    # Ensure exact candidates are in the grid
    mandatory = np.array([SQRT_3_2, SQRT_135_88, 1.2385702095, 1.2385680356])
    # Only include mandatory points that fall within range
    mandatory = mandatory[(mandatory >= alphas[0]) & (mandatory <= alphas[-1])]
    alphas = np.sort(np.unique(np.concatenate([alphas, mandatory])))
    
    print("=" * 72)
    print("  MICRO-SWEEP: Ultra-High Resolution Phase Transition Probe")
    print(f"  Grid: {len(alphas)} points in [{alphas[0]:.10f}, {alphas[-1]:.10f}]")
    print(f"  Spacing: Δα ≈ {np.mean(np.diff(alphas)):.10f}")
    print(f"  Candidates:")
    print(f"    √(3/2)    = {SQRT_3_2:.12f}")
    print(f"    √(135/88) = {SQRT_135_88:.12f}")
    print(f"    Old m=2310 crash  = 1.238570209538")
    print(f"    Old m=30030 crash = 1.238568035592")
    print("=" * 72)
    
    results = {}
    t_total = time.time()
    
    # Always do m=2310 (fast)
    results[2310] = sweep_modulus(2310, alphas, label="(FAST)")
    
    if '--both' in sys.argv:
        results[30030] = sweep_modulus(30030, alphas, label="(SLOW — ~9 hours)")
    
    total = time.time() - t_total
    
    print(f"\n{'='*72}")
    print(f"  Total: {total:.1f}s ({total/60:.1f} min)")
    print(f"{'='*72}")
    
    # Comparison table
    if len(results) >= 1:
        print(f"\n  COMPARISON TABLE:")
        print(f"  {'m':>6s}  {'crash_α':>16s}  {'err(√3/2)':>14s}  {'err(√135/88)':>14s}")
        print(f"  {'-'*6}  {'-'*16}  {'-'*14}  {'-'*14}")
        for m, r in sorted(results.items()):
            c = r["crash_50"]
            if c:
                e32 = f"{100*(c-SQRT_3_2)/SQRT_3_2:+.6f}%"
                e135 = f"{100*(c-SQRT_135_88)/SQRT_135_88:+.6f}%"
                print(f"  {m:>6d}  {c:>16.12f}  {e32:>14s}  {e135:>14s}")
    
    # Save results
    out_path = os.path.join(script_dir, "micro_sweep_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved: {out_path}")


if __name__ == "__main__":
    main()
