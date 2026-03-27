#!/usr/bin/env python3
"""
data_run_ipr_gap.py — Focused data collection for paper revision
=================================================================

Two tasks:
  Task 1: IPR at m=385 (non-primorial, phi=240) at three alpha values
  Task 2: Dense gap stability Δλ sweep at m=30030 (200 alpha values)

Task 1 output: IPR ratios for m=385 to add as non-primorial control
Task 2 output: Refine "15 sampled alpha values" to 200 for Δλ stability

Usage:
  python data_run_ipr_gap.py              # both tasks
  python data_run_ipr_gap.py --ipr-only   # Task 1 only (< 1 second)
  python data_run_ipr_gap.py --gap-only   # Task 2 only (~10 min parallel)

Author: Antonio P. Matos / Fancyland LLC
Date: March 26, 2026
"""

import os
import sys

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
if hasattr(sys.stderr, 'reconfigure'):
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')

# BLAS thread config — before numpy
_CPU = os.cpu_count() or 8
_N_WORKERS = max(1, min(8, _CPU // 4))
_BLAS = max(1, _CPU // _N_WORKERS)
os.environ["OMP_NUM_THREADS"] = str(_BLAS)
os.environ["MKL_NUM_THREADS"] = str(_BLAS)
os.environ["OPENBLAS_NUM_THREADS"] = str(_BLAS)

import numpy as np
from math import gcd
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import json
import tempfile

# Constants
SQRT_2_3    = np.sqrt(2.0 / 3.0)     # 0.8165
SQRT_3_2    = np.sqrt(3.0 / 2.0)     # 1.2247
SQRT_135_88 = np.sqrt(135.0 / 88.0)  # 1.2386

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


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
    v = np.ones(n) / np.sqrt(n)
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


def build_operator(m):
    res = coprime_residues(m)
    n = len(res)
    D = build_D_sym(res, m)
    lam_P = perron_eigenvalue(D)
    perm = inversion_perm(res, m)
    C = commutator_fast(D, perm)
    A = D / lam_P
    B = C / lam_P
    return A, B, n, lam_P


# ================================================================
# Task 1: IPR at m=385
# ================================================================

def task_ipr_m385():
    print("=" * 72)
    print("  TASK 1: IPR at m=385 (non-primorial, phi=240)")
    print("=" * 72)

    A, B, n, lam_P = build_operator(385)
    print(f"  m=385, phi={n}, lam_P={lam_P:.6f}")

    alpha_values = {
        "sqrt_2_3":    SQRT_2_3,
        "sqrt_3_2":    SQRT_3_2,
        "sqrt_135_88": SQRT_135_88,
        "alpha_c_385":  1.24227,  # crash point from job1
    }

    results = {}
    delocalized = 1.0 / n

    for label, alpha in alpha_values.items():
        t0 = time.time()
        H = A + (1j * alpha) * B
        evals, evecs = np.linalg.eigh(H)
        ipr = np.sum(np.abs(evecs)**4, axis=0)
        mean_ipr = float(np.mean(ipr))
        max_ipr = float(np.max(ipr))
        ratio = mean_ipr / delocalized
        elapsed = time.time() - t0

        results[label] = {
            "alpha": float(alpha),
            "mean_ipr": mean_ipr,
            "max_ipr": max_ipr,
            "delocalized_ref": delocalized,
            "ipr_ratio": ratio,
            "lambda_P": float(evals[-1]),
            "lambda_minus": float(evals[0]),
            "delta_lambda": float(evals[-1] - evals[0]),
            "elapsed_s": elapsed,
        }

        print(f"\n  alpha = {alpha:.10f} ({label})")
        print(f"    Mean IPR:    {mean_ipr:.8f}")
        print(f"    Max IPR:     {max_ipr:.8f}")
        print(f"    Delocal ref: {delocalized:.8f}")
        print(f"    IPR ratio:   {ratio:.4f}x")
        print(f"    lambda_P:    {evals[-1]:.8f}")
        print(f"    lambda_-:    {evals[0]:.8f}")
        print(f"    Delta_lam:   {evals[-1] - evals[0]:.8f}")
        print(f"    Time:        {elapsed:.3f}s")

    return results


# ================================================================
# Task 2: Dense gap stability at m=30030
# ================================================================

_W = {}

def _init(a_path, b_path, n):
    _W['A'] = np.load(a_path)
    _W['B'] = np.load(b_path)
    _W['n'] = n


def _gap_measure(alpha):
    """Eigenvalue solve -> extract Perron, anti-Perron, gap, variance."""
    H = _W['A'] + (1j * alpha) * _W['B']
    eigs = np.linalg.eigvalsh(H)
    lam_P = float(eigs[-1])
    lam_minus = float(eigs[0])
    delta_lam = lam_P - lam_minus
    var = float(np.var(eigs))
    # gap position
    consec = np.diff(eigs)
    j_max = int(np.argmax(consec))
    gap_pos = j_max / max(len(consec) - 1, 1)
    return (float(alpha), lam_P, lam_minus, delta_lam, var, gap_pos)


def task_gap_stability():
    print("\n" + "=" * 72)
    print("  TASK 2: Dense gap stability at m=30030 (200 alpha values)")
    print("=" * 72)

    t_build = time.time()
    A, B, n, lam_P = build_operator(30030)
    print(f"  m=30030, phi={n}, lam_P={lam_P:.6f}")
    print(f"  Build time: {time.time() - t_build:.1f}s")

    # Crash point from overnight batch
    crash_center = 1.2387666070266263
    # Window: ±0.00005 around crash (same as Phase B of overnight)
    half_width = 0.00005
    alphas = np.sort(np.unique(np.concatenate([
        np.linspace(crash_center - half_width, crash_center + half_width, 200),
        [SQRT_3_2, SQRT_135_88, crash_center],
    ])))
    n_a = len(alphas)

    print(f"  {n_a} alpha values in [{alphas[0]:.10f}, {alphas[-1]:.10f}]")
    print(f"  Workers: {_N_WORKERS}, OMP threads: {_BLAS}")

    # Save matrices for workers
    tmpdir = tempfile.mkdtemp(prefix="gap_stability_")
    a_path = os.path.join(tmpdir, "A.npy")
    b_path = os.path.join(tmpdir, "B.npy")
    np.save(a_path, A)
    np.save(b_path, B)
    del A, B  # free memory

    out = np.empty((n_a, 6))
    t0 = time.time()

    try:
        with ProcessPoolExecutor(
            max_workers=_N_WORKERS, initializer=_init,
            initargs=(a_path, b_path, n)
        ) as pool:
            futs = {pool.submit(_gap_measure, a): i for i, a in enumerate(alphas)}
            done = 0
            for f in as_completed(futs):
                idx = futs[f]
                out[idx] = f.result()
                done += 1
                if done % max(1, n_a // 10) == 0 or done == n_a:
                    el = time.time() - t0
                    eta = el / done * (n_a - done) if done < n_a else 0
                    print(f"    [{done:>4d}/{n_a}] {el:.0f}s elapsed, ~{eta:.0f}s ETA")
    finally:
        try:
            os.remove(a_path)
            os.remove(b_path)
            os.rmdir(tmpdir)
        except Exception:
            pass

    elapsed = time.time() - t0
    print(f"  Done: {elapsed:.1f}s ({elapsed/60:.1f} min)")

    # Sort by alpha
    order = np.argsort(out[:, 0])
    out = out[order]

    alphas_out   = out[:, 0]
    lam_P_out    = out[:, 1]
    lam_minus_out = out[:, 2]
    delta_lam_out = out[:, 3]
    var_out      = out[:, 4]
    gap_pos_out  = out[:, 5]

    # Statistics on delta_lambda
    mean_dl = float(np.mean(delta_lam_out))
    std_dl  = float(np.std(delta_lam_out))
    min_dl  = float(np.min(delta_lam_out))
    max_dl  = float(np.max(delta_lam_out))

    print(f"\n  Gap stability results ({n_a} points):")
    print(f"    Delta_lambda mean:   {mean_dl:.10f}")
    print(f"    Delta_lambda std:    {std_dl:.2e}")
    print(f"    Delta_lambda range:  [{min_dl:.10f}, {max_dl:.10f}]")
    print(f"    Delta_lambda spread: {max_dl - min_dl:.2e}")

    # Variance stats
    mean_var = float(np.mean(var_out))
    print(f"    Variance mean:       {mean_var:.8f}")
    print(f"    Var * phi:           {mean_var * n:.4f}")

    results = {
        "n_points": n_a,
        "alpha_range": [float(alphas_out[0]), float(alphas_out[-1])],
        "delta_lambda": {
            "mean": mean_dl,
            "std": std_dl,
            "min": min_dl,
            "max": max_dl,
            "spread": float(max_dl - min_dl),
        },
        "variance_mean": mean_var,
        "var_times_phi": float(mean_var * n),
        "elapsed_s": elapsed,
        "raw_data": {
            "alphas": alphas_out.tolist(),
            "delta_lambda": delta_lam_out.tolist(),
            "lam_P": lam_P_out.tolist(),
            "lam_minus": lam_minus_out.tolist(),
            "gap_pos": gap_pos_out.tolist(),
            "variance": var_out.tolist(),
        },
    }

    return results


# ================================================================
# Main
# ================================================================

def main():
    do_ipr = '--gap-only' not in sys.argv
    do_gap = '--ipr-only' not in sys.argv

    all_results = {}

    if do_ipr:
        all_results["ipr_m385"] = task_ipr_m385()

    if do_gap:
        all_results["gap_stability_m30030"] = task_gap_stability()

    # Save results
    out_path = os.path.join(SCRIPT_DIR, "data_run_results.json")
    with open(out_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\n  Results saved to: {out_path}")


if __name__ == "__main__":
    main()
