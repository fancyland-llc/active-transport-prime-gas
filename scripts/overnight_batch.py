#!/usr/bin/env python3
"""
overnight_batch.py — Comprehensive overnight computation batch
==============================================================

Five jobs, parallelized across alpha-values:

  Job 1: m=385   coarse sweep       (non-primorial control, φ=240)
  Job 2: m=2310  post-transition    (α from 1.25 to 2.5)
  Job 3: m=30030 micro-sweep        (THE TEST: √(135/88)?)
  Job 4: m=15015 coarse+micro       (non-primorial φ=5760 control)
  Job 5: Eigenvector structure dump  (at α_c for m=2310 + m=30030)

Uses ProcessPoolExecutor. Each worker gets OMP_NUM_THREADS = cpu//workers
so the full CPU is utilized without oversubscription.

Usage:
  python overnight_batch.py --calibrate         # benchmark (30s)
  python overnight_batch.py --quick             # jobs 1+2 only (~10min)
  python overnight_batch.py --medium            # jobs 1+2+3 (~hours)
  python overnight_batch.py                     # all 5 jobs (overnight)

  set OVERNIGHT_WORKERS=4   (override worker count)

Outputs saved to: overnight_results/  (one .npz + one .json per job)

Author: Antonio P. Matos / Fancyland LLC
Date: March 25, 2026
"""

import os
import sys

# Force UTF-8 output on Windows (cp1252 can't handle √, α, φ, × chars)
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
if hasattr(sys.stderr, 'reconfigure'):
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')

# ================================================================
# BLAS thread configuration — MUST be before numpy import
# ================================================================
_CPU = os.cpu_count() or 8
_N_WORKERS = int(os.environ.get("OVERNIGHT_WORKERS", "0"))
if _N_WORKERS <= 0:
    _N_WORKERS = max(1, min(8, _CPU // 4))
_BLAS = max(1, _CPU // _N_WORKERS)
os.environ["OMP_NUM_THREADS"] = str(_BLAS)
os.environ["MKL_NUM_THREADS"] = str(_BLAS)
os.environ["OPENBLAS_NUM_THREADS"] = str(_BLAS)

import numpy as np
from math import gcd
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
import time
import json
import tempfile
import traceback
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ================================================================
# Constants
# ================================================================
SQRT_3_2    = np.sqrt(3.0 / 2.0)     # 1.224744871391589
SQRT_2_3    = np.sqrt(2.0 / 3.0)     # 0.816496580927726
SQRT_135_88 = np.sqrt(135.0 / 88.0)  # 1.238584235767...

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR    = os.path.join(SCRIPT_DIR, "overnight_results")


# ================================================================
# Utility functions
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


def gap_position(eigs):
    consec = np.diff(eigs)
    j_max = int(np.argmax(consec))
    return j_max / max(len(consec) - 1, 1)


def chirality(m):
    primes = []
    t = m
    for p in range(2, int(t**0.5) + 1):
        if t % p == 0:
            primes.append(p)
            while t % p == 0:
                t //= p
    if t > 1:
        primes.append(t)
    odd = [p for p in primes if p > 2]
    np_ = sum(1 for p in odd if p % 4 == 1)
    nm_ = sum(1 for p in odd if p % 4 == 3)
    return abs(np_ - nm_), np_, nm_, primes


def build_operator(m):
    """Build normalized operator components (A, B) and metadata."""
    t0 = time.time()
    res = coprime_residues(m)
    n = len(res)
    D = build_D_sym(res, m)
    lam_P = perron_eigenvalue(D)
    perm = inversion_perm(res, m)
    C = commutator_fast(D, perm)

    norm_D = float(np.linalg.norm(D, 'fro'))
    norm_C = float(np.linalg.norm(C, 'fro'))

    A = D / lam_P
    B = C / lam_P
    del D, C, perm

    I_val, np_, nm_, primes = chirality(m)
    meta = {
        "m": m, "phi": n, "I": I_val,
        "lam_P": float(lam_P),
        "norm_D": norm_D, "norm_C": norm_C,
        "comm_ratio": norm_C / norm_D,
        "coupling": norm_C / lam_P,
        "primes": primes,
        "build_time": time.time() - t0,
    }
    return A, B, meta


# ================================================================
# Worker functions (for ProcessPoolExecutor)
# ================================================================

_W = {}  # worker-local state


def _init(a_path, b_path, n):
    _W['A'] = np.load(a_path)
    _W['B'] = np.load(b_path)
    _W['n'] = n


def _eigvals(alpha):
    """Eigenvalues only — fast."""
    H = _W['A'] + (1j * alpha) * _W['B']
    eigs = np.linalg.eigvalsh(H)
    gp = gap_position(eigs)
    pgap = float(eigs[-1] - eigs[-2])
    consec = np.diff(eigs)
    mig = float(np.max(consec))
    var = float(np.var(eigs))
    top5 = eigs[-5:][::-1].tolist()
    return (float(alpha), gp, pgap, mig, var, top5)


def _eigvecs(alpha):
    """Full eigenvectors + eigenvalues — expensive."""
    n = _W['n']
    H = _W['A'] + (1j * alpha) * _W['B']
    evals, evecs = np.linalg.eigh(H)
    ipr = np.sum(np.abs(evecs)**4, axis=0)
    return {
        "alpha": float(alpha),
        "eigenvalues": evals.tolist(),
        "ipr_per_mode": ipr.tolist(),
        "mean_ipr": float(np.mean(ipr)),
        "ipr_ratio": float(np.mean(ipr) * n),
        "perron_evec_abs": np.abs(evecs[:, -1]).tolist(),
        "anti_perron_evec_abs": np.abs(evecs[:, 0]).tolist(),
    }


# ================================================================
# Parallel sweep engine
# ================================================================

def save_AB_to_temp(A, B):
    """Save A, B to temp .npy files for worker initialization."""
    tmpdir = tempfile.mkdtemp(prefix="overnight_")
    a_path = os.path.join(tmpdir, "A.npy")
    b_path = os.path.join(tmpdir, "B.npy")
    np.save(a_path, A)
    np.save(b_path, B)
    return a_path, b_path, tmpdir


def cleanup_temp(a_path, b_path, tmpdir):
    try:
        os.remove(a_path)
        os.remove(b_path)
        os.rmdir(tmpdir)
    except Exception:
        pass


def parallel_eigvals(A, B, alphas, workers=_N_WORKERS, label=""):
    """Parallel eigenvalue sweep. Returns structured arrays."""
    n = A.shape[0]
    n_a = len(alphas)
    if label:
        print(f"  [{label}] {n_a} α × {n}×{n}, {workers} workers, "
              f"OMP={os.environ.get('OMP_NUM_THREADS','?')}")
    else:
        print(f"  {n_a} α × {n}×{n}, {workers} workers")

    a_path, b_path, tmpdir = save_AB_to_temp(A, B)
    out_alpha = np.empty(n_a)
    out_gpos  = np.empty(n_a)
    out_pgap  = np.empty(n_a)
    out_mig   = np.empty(n_a)
    out_var   = np.empty(n_a)
    out_top5  = np.empty((n_a, 5))

    t0 = time.time()
    try:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=_init,
            initargs=(a_path, b_path, n)
        ) as pool:
            futs = {pool.submit(_eigvals, a): i for i, a in enumerate(alphas)}
            done = 0
            for f in as_completed(futs):
                idx = futs[f]
                alpha, gp, pg, mg, v, t5 = f.result()
                out_alpha[idx] = alpha
                out_gpos[idx]  = gp
                out_pgap[idx]  = pg
                out_mig[idx]   = mg
                out_var[idx]   = v
                out_top5[idx]  = t5
                done += 1
                if done % max(1, n_a // 10) == 0 or done == n_a:
                    el = time.time() - t0
                    eta = el / done * (n_a - done)
                    print(f"    [{done:>5d}/{n_a}] {el:.0f}s, ~{eta:.0f}s ETA")
    finally:
        cleanup_temp(a_path, b_path, tmpdir)

    # Sort by alpha
    order = np.argsort(out_alpha)
    elapsed = time.time() - t0
    print(f"  Done: {elapsed:.1f}s ({elapsed/60:.1f} min), {elapsed/n_a:.2f}s/α eff")

    return {
        "alphas":     out_alpha[order],
        "gap_pos":    out_gpos[order],
        "perron_gap": out_pgap[order],
        "max_int_gap": out_mig[order],
        "eig_var":    out_var[order],
        "top5":       out_top5[order],
        "var_times_phi": out_var[order] * n,
        "elapsed_s":  elapsed,
    }


def parallel_eigvecs(A, B, alphas, workers=_N_WORKERS, label=""):
    """Parallel eigenvector sweep. Returns list of dicts."""
    n = A.shape[0]
    print(f"  [{label}] {len(alphas)} α (full eigvec), {workers} workers")
    a_path, b_path, tmpdir = save_AB_to_temp(A, B)
    results = {}
    t0 = time.time()
    try:
        with ProcessPoolExecutor(
            max_workers=min(workers, len(alphas)),
            initializer=_init, initargs=(a_path, b_path, n)
        ) as pool:
            futs = {pool.submit(_eigvecs, a): a for a in alphas}
            for f in as_completed(futs):
                r = f.result()
                results[r["alpha"]] = r
                print(f"    α={r['alpha']:.8f}: IPR ratio={r['ipr_ratio']:.3f}x")
    finally:
        cleanup_temp(a_path, b_path, tmpdir)
    print(f"  Done: {time.time()-t0:.1f}s")
    return results


# ================================================================
# Analysis helpers
# ================================================================

def find_crash(alphas, gpos, threshold=0.5):
    for i in range(len(gpos) - 1):
        if gpos[i] > threshold and gpos[i+1] <= threshold:
            frac = (threshold - gpos[i]) / (gpos[i+1] - gpos[i])
            return float(alphas[i] + frac * (alphas[i+1] - alphas[i]))
    return None


def transition_metrics(sweep, phi):
    """Extract transition analysis from a sweep result dict."""
    a = sweep["alphas"]
    g = sweep["gap_pos"]
    c50 = find_crash(a, g, 0.5)
    c75 = find_crash(a, g, 0.75)
    c25 = find_crash(a, g, 0.25)

    derivs = np.diff(g) / np.diff(a)
    j = int(np.argmin(derivs))
    steep_a = float(0.5 * (a[j] + a[j+1]))
    steep_d = float(derivs[j])

    width = float(c25 - c75) if (c25 and c75) else None
    center = float(0.5 * (c25 + c75)) if (c25 and c75) else None

    info = {
        "crash_50": c50, "crash_75": c75, "crash_25": c25,
        "steepest_alpha": steep_a, "steepest_deriv": steep_d,
        "width_75_25": width, "center": center,
    }
    if c50:
        info["err_sqrt32"]    = float(c50 - SQRT_3_2)
        info["err_sqrt32_pct"] = float(100*(c50 - SQRT_3_2)/SQRT_3_2)
        info["err_sqrt135"]    = float(c50 - SQRT_135_88)
        info["err_sqrt135_pct"] = float(100*(c50 - SQRT_135_88)/SQRT_135_88)
    return info


def print_transition(info, phi, label=""):
    if label:
        print(f"\n  === {label} ===")
    c = info.get("crash_50")
    if c is None:
        print("  No crash detected in sweep range")
        return
    print(f"  Crash (gap_pos=0.5): α = {c:.12f}")
    print(f"    From √(3/2):    {info['err_sqrt32']:+.12f} ({info['err_sqrt32_pct']:+.8f}%)")
    print(f"    From √(135/88): {info['err_sqrt135']:+.12f} ({info['err_sqrt135_pct']:+.8f}%)")
    w = info.get("width_75_25")
    if w:
        print(f"  Transition width (75→25): {w:.2e}")
        print(f"  Transition center:        {info['center']:.12f}")
    print(f"  Steepest descent: α = {info['steepest_alpha']:.12f}, "
          f"d/dα = {info['steepest_deriv']:.1f}")


# ================================================================
# Job 1: m=385 coarse (small non-primorial control)
# ================================================================

def job_m385():
    print("\n" + "=" * 72)
    print("  JOB 1: m=385 coarse sweep (non-primorial, φ=240)")
    print("=" * 72)

    A, B, meta = build_operator(385)
    phi = meta["phi"]
    print(f"  m=385, φ={phi}, I={meta['I']}, primes={meta['primes']}")
    print(f"  λ_P = {meta['lam_P']:.6f}, ||[D,τ]||/||D|| = {meta['comm_ratio']:.6f}")

    # Coarse sweep
    alphas = np.sort(np.unique(np.concatenate([
        np.linspace(0.5, 1.45, 222),
        [SQRT_2_3, SQRT_3_2, SQRT_135_88],
    ])))
    sweep = parallel_eigvals(A, B, alphas, label="m=385 coarse")
    info = transition_metrics(sweep, phi)
    print_transition(info, phi, "m=385 RESULT")

    # If crash found, do a quick micro-sweep
    micro_sweep = None
    micro_info = None
    if info["crash_50"]:
        center = info["crash_50"]
        micro_a = np.sort(np.unique(np.concatenate([
            np.linspace(center - 0.003, center + 0.003, 500),
            [SQRT_3_2, SQRT_135_88],
        ])))
        micro_a = micro_a[(micro_a > 0) & (micro_a < 3.0)]
        micro_sweep = parallel_eigvals(A, B, micro_a, label="m=385 micro")
        micro_info = transition_metrics(micro_sweep, phi)
        print_transition(micro_info, phi, "m=385 MICRO RESULT")

    result = {
        "meta": meta, "coarse_sweep": sweep, "coarse_info": info,
        "micro_sweep": micro_sweep, "micro_info": micro_info,
    }
    save_job("job1_m385", result, meta)
    return result


# ================================================================
# Job 2: m=2310 post-transition (α > 1.24)
# ================================================================

def job_post_transition():
    print("\n" + "=" * 72)
    print("  JOB 2: m=2310 post-transition sweep (α = 1.24 to 2.5)")
    print("=" * 72)

    A, B, meta = build_operator(2310)
    phi = meta["phi"]
    print(f"  m=2310, φ={phi}, I={meta['I']}")

    alphas = np.linspace(1.24, 2.5, 300)
    sweep = parallel_eigvals(A, B, alphas, label="m=2310 post-transition")

    # Find any secondary features
    g = sweep["gap_pos"]
    dg = np.diff(g) / np.diff(alphas)
    # Look for inflection points (sign changes in 2nd derivative)
    ddg = np.diff(dg) / np.diff(0.5*(alphas[:-1] + alphas[1:]))
    sign_changes = np.where(np.diff(np.sign(ddg)))[0]

    print(f"\n  Post-transition features:")
    print(f"    Gap pos range: [{np.min(g):.4f}, {np.max(g):.4f}]")
    print(f"    Var·φ range:   [{np.min(sweep['var_times_phi']):.4f}, "
          f"{np.max(sweep['var_times_phi']):.4f}]")
    print(f"    Inflection points: {len(sign_changes)}")
    for sc in sign_changes[:5]:
        print(f"      α ≈ {alphas[sc+1]:.4f}")

    # Perron gap trajectory
    pg = sweep["perron_gap"]
    print(f"    Perron gap: {pg[0]:.6f} → {pg[-1]:.6f}")
    print(f"    Perron gap min: {np.min(pg):.6f} at α={alphas[np.argmin(pg)]:.4f}")

    result = {"meta": meta, "sweep": sweep}
    save_job("job2_post_transition", result, meta)
    return result


# ================================================================
# Job 3: m=30030 micro-sweep (THE TEST)
# ================================================================

def job_m30030_micro():
    print("\n" + "=" * 72)
    print("  JOB 3: m=30030 MICRO-SWEEP — THE TEST")
    print("  Is α_c = √(135/88) = {:.12f} ?".format(SQRT_135_88))
    print("=" * 72)

    A, B, meta = build_operator(30030)
    phi = meta["phi"]
    print(f"  m=30030, φ={phi}, I={meta['I']}")
    print(f"  λ_P = {meta['lam_P']:.6f}, ||[D,τ]||/||D|| = {meta['comm_ratio']:.6f}")
    print(f"  κ = ||[D,τ]||/λ_P = {meta['coupling']:.8f}")

    # Phase A: locate crash at moderate resolution (200 pts in ±0.004)
    print("\n  --- Phase A: locate crash (200 pts) ---")
    phase_a = np.sort(np.unique(np.concatenate([
        np.linspace(SQRT_135_88 - 0.004, SQRT_135_88 + 0.004, 200),
        [SQRT_3_2, SQRT_135_88],
    ])))
    sweep_a = parallel_eigvals(A, B, phase_a, label="Phase A (coarse micro)")
    info_a = transition_metrics(sweep_a, phi)
    print_transition(info_a, phi, "Phase A RESULT")

    crash_est = info_a.get("crash_50", SQRT_135_88)

    # Phase B: ultra-fine around the crash (500 pts in ±0.00005)
    print("\n  --- Phase B: ultra-micro (500 pts around crash) ---")
    phase_b = np.sort(np.unique(np.concatenate([
        np.linspace(crash_est - 0.00005, crash_est + 0.00005, 500),
        [SQRT_3_2, SQRT_135_88, crash_est],
    ])))
    phase_b = phase_b[(phase_b > 0)]
    sweep_b = parallel_eigvals(A, B, phase_b, label="Phase B (ultra-micro)")
    info_b = transition_metrics(sweep_b, phi)
    print_transition(info_b, phi, "Phase B RESULT")

    # Phase C: eigenvectors in the transition window (Gemini's asks)
    evec_results = None
    crash_final = info_b.get("crash_50", crash_est)
    width = info_b.get("width_75_25")
    if width and width > 0:
        # 15 α values spanning the transition window
        evec_alphas = np.linspace(crash_final - width, crash_final + width, 15)
    else:
        # fallback: 15 points in ±1e-5 around crash
        evec_alphas = np.linspace(crash_final - 1e-5, crash_final + 1e-5, 15)

    print(f"\n  --- Phase C: eigenvectors at {len(evec_alphas)} points ---")
    evec_results = parallel_eigvecs(A, B, evec_alphas, label="Phase C eigvecs")

    # Full eigenvalue spectrum dump at the crash point
    print(f"\n  --- Full spectrum at α_c = {crash_final:.10f} ---")
    evec_at_crash = evec_results.get(
        min(evec_results.keys(), key=lambda a: abs(a - crash_final)), None)
    if evec_at_crash:
        evals = np.array(evec_at_crash["eigenvalues"])
        print(f"  Eigenvalue range: [{evals[0]:.6f}, {evals[-1]:.6f}]")
        print(f"  Variance: {np.var(evals):.8f}, Var·φ = {np.var(evals)*phi:.4f}")
        print(f"  IPR ratio: {evec_at_crash['ipr_ratio']:.4f}x")

    result = {
        "meta": meta,
        "phase_a_sweep": sweep_a, "phase_a_info": info_a,
        "phase_b_sweep": sweep_b, "phase_b_info": info_b,
        "evec_results": evec_results,
    }
    save_job("job3_m30030_micro", result, meta)
    return result


# ================================================================
# Job 4: m=15015 coarse + micro (non-primorial control)
# ================================================================

def job_m15015():
    print("\n" + "=" * 72)
    print("  JOB 4: m=15015 (non-primorial, φ=5760)")
    print("  Question: does √(135/88) hold without factor of 2?")
    print("=" * 72)

    A, B, meta = build_operator(15015)
    phi = meta["phi"]
    print(f"  m=15015, φ={phi}, I={meta['I']}, primes={meta['primes']}")
    print(f"  λ_P = {meta['lam_P']:.6f}, ||[D,τ]||/||D|| = {meta['comm_ratio']:.6f}")
    print(f"  κ = {meta['coupling']:.8f}")

    # Phase A: coarse to find IF and WHERE there's a crash
    print("\n  --- Phase A: coarse sweep (222 pts) ---")
    coarse = np.sort(np.unique(np.concatenate([
        np.linspace(0.5, 2.0, 222),
        [SQRT_2_3, SQRT_3_2, SQRT_135_88],
    ])))
    sweep_c = parallel_eigvals(A, B, coarse, label="m=15015 coarse")
    info_c = transition_metrics(sweep_c, phi)
    print_transition(info_c, phi, "m=15015 COARSE RESULT")

    micro_sweep = None
    micro_info = None
    evec_results = None

    # Phase B: micro around crash if found
    if info_c.get("crash_50"):
        crash_coarse = info_c["crash_50"]
        print(f"\n  --- Phase B: micro sweep (500 pts around α={crash_coarse:.6f}) ---")
        micro_a = np.sort(np.unique(np.concatenate([
            np.linspace(crash_coarse - 0.004, crash_coarse + 0.004, 500),
            [SQRT_3_2, SQRT_135_88],
        ])))
        micro_a = micro_a[(micro_a > 0)]
        micro_sweep = parallel_eigvals(A, B, micro_a, label="m=15015 micro")
        micro_info = transition_metrics(micro_sweep, phi)
        print_transition(micro_info, phi, "m=15015 MICRO RESULT")

        # Phase C: eigenvectors at crash
        crash_micro = micro_info.get("crash_50", crash_coarse)
        evec_alphas = np.linspace(crash_micro - 1e-5, crash_micro + 1e-5, 5)
        print(f"\n  --- Phase C: eigenvectors at crash ---")
        evec_results = parallel_eigvecs(A, B, evec_alphas, label="m=15015 eigvecs")
    else:
        print("\n  *** No crash detected in [0.5, 2.0]. Extending to [2.0, 4.0]... ***")
        ext = np.linspace(2.0, 4.0, 100)
        sweep_ext = parallel_eigvals(A, B, ext, label="m=15015 extended")
        info_ext = transition_metrics(sweep_ext, phi)
        print_transition(info_ext, phi, "m=15015 EXTENDED")

    result = {
        "meta": meta,
        "coarse_sweep": sweep_c, "coarse_info": info_c,
        "micro_sweep": micro_sweep, "micro_info": micro_info,
        "evec_results": evec_results,
    }
    save_job("job4_m15015", result, meta)
    return result


# ================================================================
# Job 5: Eigenvector structure at α_c (double helix hunt)
# ================================================================

def job_eigenvectors(crash_alpha_30030):
    print("\n" + "=" * 72)
    print("  JOB 5: Eigenvector Structure at α_c")
    print("=" * 72)

    results = {}

    # --- m=2310 (fast) ---
    print("\n  --- m=2310: 9 eigenvector points around crash ---")
    A2, B2, meta2 = build_operator(2310)
    # Our micro-sweep found m=2310 crash at ~1.238922
    crash_2310 = 1.238922
    evec_a = np.linspace(crash_2310 - 2e-5, crash_2310 + 2e-5, 9)
    evec_2310 = parallel_eigvecs(A2, B2, evec_a, label="m=2310 eigvecs")
    results["m2310"] = {"meta": meta2, "evec": evec_2310}
    del A2, B2

    # --- m=30030 (expensive, single point) ---
    if crash_alpha_30030:
        print(f"\n  --- m=30030: eigenvectors at α_c = {crash_alpha_30030:.10f} ---")
        A3, B3, meta3 = build_operator(30030)
        # Just 3 points: slightly below, at, and slightly above
        evec_a3 = [crash_alpha_30030 - 1e-6, crash_alpha_30030, crash_alpha_30030 + 1e-6]
        evec_30030 = parallel_eigvecs(A3, B3, evec_a3, label="m=30030 eigvecs")
        results["m30030"] = {"meta": meta3, "evec": evec_30030}
        del A3, B3
    else:
        print("  *** Skipping m=30030 eigvecs (no crash detected) ***")

    save_job("job5_eigenvectors", results, {})
    return results


# ================================================================
# Save / checkpoint helpers
# ================================================================

def make_serializable(obj):
    """Convert numpy arrays/types to JSON-serializable types."""
    if isinstance(obj, dict):
        return {str(k): make_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [make_serializable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def save_job(name, result, meta):
    """Save job results as .json checkpoint."""
    os.makedirs(OUT_DIR, exist_ok=True)
    json_path = os.path.join(OUT_DIR, f"{name}.json")
    with open(json_path, "w") as f:
        json.dump(make_serializable(result), f, indent=1)
    print(f"\n  ✓ Saved: {json_path} ({os.path.getsize(json_path):,} bytes)")

    # Also save key arrays as npz for fast numpy loading
    npz_data = {}
    for key, val in result.items():
        if isinstance(val, dict):
            for k2, v2 in val.items():
                if isinstance(v2, np.ndarray):
                    npz_data[f"{key}__{k2}"] = v2
    if npz_data:
        npz_path = os.path.join(OUT_DIR, f"{name}.npz")
        np.savez_compressed(npz_path, **npz_data)
        print(f"  ✓ Saved: {npz_path} ({os.path.getsize(npz_path):,} bytes)")


# ================================================================
# Calibration
# ================================================================

def calibrate():
    print("=" * 72)
    print("  CALIBRATION — Benchmarking eigendecomposition")
    print(f"  CPU cores: {_CPU}")
    print(f"  Workers: {_N_WORKERS}")
    print(f"  OMP_NUM_THREADS: {os.environ.get('OMP_NUM_THREADS', '?')}")
    print("=" * 72)

    for m_test, name in [(2310, "m=2310"), (30030, "m=30030")]:
        print(f"\n  --- {name} ---")
        t0 = time.time()
        A, B, meta = build_operator(m_test)
        n = meta["phi"]
        print(f"  Build: {time.time()-t0:.2f}s, φ={n}")

        H = A + (1j * SQRT_135_88) * B

        # Warm up
        _ = np.linalg.eigvalsh(H)

        # Eigenvalues only (eigvalsh)
        times_valsh = []
        for _ in range(3):
            t0 = time.time()
            _ = np.linalg.eigvalsh(H)
            times_valsh.append(time.time() - t0)
        t_valsh = np.median(times_valsh)

        # Eigenvalues + vectors (eigh) — only for m=2310
        t_eigh = None
        if n <= 1000:
            times_eigh = []
            for _ in range(2):
                t0 = time.time()
                _ = np.linalg.eigh(H)
                times_eigh.append(time.time() - t0)
            t_eigh = np.median(times_eigh)

        print(f"  eigvalsh: {t_valsh:.3f}s")
        if t_eigh:
            print(f"  eigh:     {t_eigh:.3f}s")

        # Estimate parallel sweep time
        for n_pts in [200, 500, 2000]:
            serial = n_pts * t_valsh
            parallel = n_pts * t_valsh / _N_WORKERS
            print(f"  {n_pts} pts: serial={serial/60:.1f}min, "
                  f"{_N_WORKERS} workers={parallel/60:.1f}min")

        del A, B, H

    # Also calibrate m=15015
    print(f"\n  --- m=15015 ---")
    t0 = time.time()
    A, B, meta = build_operator(15015)
    n = meta["phi"]
    print(f"  Build: {time.time()-t0:.2f}s, φ={n}")

    H = A + (1j * SQRT_135_88) * B
    _ = np.linalg.eigvalsh(H)  # warm up
    times = []
    for _ in range(3):
        t0 = time.time()
        _ = np.linalg.eigvalsh(H)
        times.append(time.time() - t0)
    t_m = np.median(times)
    print(f"  eigvalsh: {t_m:.3f}s")
    for n_pts in [222, 500]:
        parallel = n_pts * t_m / _N_WORKERS
        print(f"  {n_pts} pts: {_N_WORKERS} workers={parallel/60:.1f}min")

    print("\n  === TOTAL OVERNIGHT ESTIMATE ===")
    # Rough estimate based on calibration
    # m=30030: 200 (phase A) + 500 (phase B) + 15 eigvec
    # m=15015: 222 (coarse) + 500 (micro) + 5 eigvec
    # m=385:   222 + 500 (trivial)
    # m=2310:  300 + 9 eigvec (trivial)
    # (we don't have actual t_valsh for m=30030 if it errored)


# ================================================================
# Main
# ================================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("=" * 72)
    print("  OVERNIGHT BATCH — Active Transport Phase Transition Hunt")
    print(f"  Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  CPU: {_CPU} cores, Workers: {_N_WORKERS}, OMP: {_BLAS}")
    print(f"  Target: √(135/88) = {SQRT_135_88:.12f}")
    print(f"  Output: {OUT_DIR}")
    print("=" * 72)

    if '--calibrate' in sys.argv:
        calibrate()
        return

    mode = 'full'
    if '--quick' in sys.argv:
        mode = 'quick'
    elif '--medium' in sys.argv:
        mode = 'medium'

    t_total = time.time()
    all_results = {}

    # === Fast jobs first ===
    all_results['m385'] = job_m385()
    all_results['post_transition'] = job_post_transition()

    if mode == 'quick':
        total = time.time() - t_total
        print(f"\n  Quick mode complete: {total:.0f}s ({total/60:.1f} min)")
        save_summary(all_results)
        return

    # === The big test ===
    all_results['m30030'] = job_m30030_micro()

    if mode == 'medium':
        total = time.time() - t_total
        print(f"\n  Medium mode complete: {total:.0f}s ({total/60:.1f} min)")
        save_summary(all_results)
        return

    # === Non-primorial control ===
    all_results['m15015'] = job_m15015()

    # === Eigenvector dump ===
    crash_30030 = None
    m30030_result = all_results.get('m30030', {})
    for phase in ['phase_b_info', 'phase_a_info']:
        info = m30030_result.get(phase, {})
        if info and info.get('crash_50'):
            crash_30030 = info['crash_50']
            break
    all_results['eigenvectors'] = job_eigenvectors(crash_30030)

    total = time.time() - t_total
    print(f"\n{'='*72}")
    print(f"  OVERNIGHT BATCH COMPLETE: {total:.0f}s ({total/3600:.2f} hours)")
    print(f"{'='*72}")
    save_summary(all_results)


def save_summary(all_results):
    """Final summary JSON with key results across all jobs."""
    summary = {
        "timestamp": time.strftime('%Y-%m-%d %H:%M:%S'),
        "workers": _N_WORKERS,
        "omp_threads": _BLAS,
        "target": float(SQRT_135_88),
        "jobs_completed": list(all_results.keys()),
    }

    # Extract crash points for comparison
    crash_table = {}
    for job_key, result in all_results.items():
        if not isinstance(result, dict):
            continue
        for info_key in ['micro_info', 'phase_b_info', 'phase_a_info', 'coarse_info']:
            info = result.get(info_key) if isinstance(result.get(info_key), dict) else None
            if info and info.get('crash_50'):
                meta = result.get('meta', {})
                crash_table[job_key] = {
                    "m": meta.get("m"),
                    "phi": meta.get("phi"),
                    "crash_50": info["crash_50"],
                    "err_sqrt32_pct": info.get("err_sqrt32_pct"),
                    "err_sqrt135_pct": info.get("err_sqrt135_pct"),
                    "width": info.get("width_75_25"),
                    "steepest_deriv": info.get("steepest_deriv"),
                    "source": info_key,
                }
                break

    summary["crash_table"] = crash_table

    # Print comparison
    print("\n  ╔══════════════════════════════════════════════════════════════════╗")
    print("  ║               OVERNIGHT RESULTS — CRASH TABLE                  ║")
    print("  ╠══════════════════════════════════════════════════════════════════╣")
    print("  ║  {:>7s}  {:>5s}  {:>14s}  {:>10s}  {:>10s}  {:>10s} ║".format(
        "m", "φ", "crash_α", "err(√3/2)", "err(135)", "width"))
    print("  ╠══════════════════════════════════════════════════════════════════╣")
    for key, row in sorted(crash_table.items(), key=lambda x: x[1].get("m", 0)):
        m = row.get("m", "?")
        phi = row.get("phi", "?")
        ca = row["crash_50"]
        e32 = row.get("err_sqrt32_pct")
        e135 = row.get("err_sqrt135_pct")
        w = row.get("width")
        print("  ║  {:>7s}  {:>5s}  {:>14.10f}  {:>+9.5f}%  {:>+9.5f}%  {:>10s} ║".format(
            str(m), str(phi), ca,
            e32 if e32 else 0,
            e135 if e135 else 0,
            f"{w:.2e}" if w else "N/A",
        ))
    print("  ╚══════════════════════════════════════════════════════════════════╝")

    summary_path = os.path.join(OUT_DIR, "overnight_summary.json")
    with open(summary_path, "w") as f:
        json.dump(make_serializable(summary), f, indent=2)
    print(f"\n  Final summary: {summary_path}")


if __name__ == '__main__':
    mp.freeze_support()
    main()
