#!/usr/bin/env python3
"""
cache_results.py — Save all locked-in computation results to disk
=================================================================

Reads phase_transition_summary.txt and chiral_quench_summary.txt,
packages everything into a single .npz cache file that downstream
analytics can load in milliseconds.

Run once after computations finish.  Never recompute — just load.

Usage:
  python cache_results.py

Output:
  active_transport_cache.npz   — full sweep data cache
  active_transport_cache.json  — human-readable summary

Author: Antonio P. Matos / Fancyland LLC
Date: March 25, 2026
"""

import numpy as np
import json
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ================================================================
#  LOCKED-IN DATA — from chiral_phase_transition.py (full run)
# ================================================================

# Constants
SQRT_3_2 = np.sqrt(3.0 / 2.0)   # 1.224744871391589
SQRT_2_3 = np.sqrt(2.0 / 3.0)   # 0.816496580927726

# Phase transition crash points (gap_pos = 0.5 interpolation)
crash_data = {
    30:    {"phi": 8,    "I": 0, "crash_alpha": None,             "steepest_alpha": 0.5128947368, "steepest_deriv": 0.00},
    210:   {"phi": 48,   "I": 1, "crash_alpha": 1.2702512563,     "steepest_alpha": 1.2702261307, "steepest_deriv": -432.61},
    2310:  {"phi": 480,  "I": 2, "crash_alpha": 1.2385702095,     "steepest_alpha": 1.2385678392, "steepest_deriv": -441.30},
    30030: {"phi": 5760, "I": 1, "crash_alpha": 1.2385680356,     "steepest_alpha": 1.2385678392, "steepest_deriv": -442.15},
}

# IPR data at alpha = sqrt(3/2)
ipr_data = {
    30:    {"alpha": SQRT_3_2, "mean_ipr": 0.18129892, "ratio": 1.45},
    210:   {"alpha": SQRT_3_2, "mean_ipr": 0.03909618, "ratio": 1.88},
    2310:  {"alpha": SQRT_3_2, "mean_ipr": 0.00475671, "ratio": 2.28},
    30030: {"alpha": SQRT_3_2, "mean_ipr": 0.00045668, "ratio": 2.63},
}

# Flat-band condensation at alpha = sqrt(2/3)
flatband_data = {
    30:    {"phi": 8,    "var": 0.203022, "var_times_phi": 1.624},
    210:   {"phi": 48,   "var": 0.036839, "var_times_phi": 1.768},
    2310:  {"phi": 480,  "var": 0.003702, "var_times_phi": 1.777},
    30030: {"phi": 5760, "var": 0.000309, "var_times_phi": 1.780},
}

# Perron gap convergence at alpha = sqrt(2/3)
perron_gap_data = {
    30:    {"gap": 0.856826, "width": 1.634124, "ratio": 0.52433},
    210:   {"gap": 0.813760, "width": 1.603232, "ratio": 0.50757},
    2310:  {"gap": 0.814621, "width": 1.591708, "ratio": 0.51179},
    30030: {"gap": 0.814642, "width": 1.590752, "ratio": 0.51211},
}

# Commutator ratio convergence
comm_ratio_data = {
    30:    {"ratio": 0.528498, "ratio_sq": 0.279310},
    210:   {"ratio": 0.697550, "ratio_sq": 0.486576},
    2310:  {"ratio": 0.706280, "ratio_sq": 0.498831},
    30030: {"ratio": 0.706982, "ratio_sq": 0.499823},
}

# Spectral measure convergence (anti-Perron eigenvalue)
spectral_measure = {
    30:    {"lambda_max": 1.0000, "lambda_min": -0.6341, "frac_near_zero": 0.375},
    210:   {"lambda_max": 1.0000, "lambda_min": -0.6032, "frac_near_zero": 0.875},
    2310:  {"lambda_max": 1.0000, "lambda_min": -0.5917, "frac_near_zero": 0.983},
    30030: {"lambda_max": 1.0000, "lambda_min": -0.5908, "frac_near_zero": 0.999},
}

# Eigenvalue at alpha = sqrt(3/2) from phase transition sweep
at_sqrt32 = {
    30:    {"gap_pos": 1.000000, "perron_gap": 0.726312, "max_int_gap": 0.726312, "eig_var": 0.24285904},
    210:   {"gap_pos": 1.000000, "perron_gap": 0.663456, "max_int_gap": 0.663456, "eig_var": 0.04811722},
    2310:  {"gap_pos": 1.000000, "perron_gap": 0.666355, "max_int_gap": 0.666355, "eig_var": 0.00485637},
    30030: {"gap_pos": 1.000000, "perron_gap": 0.666526, "max_int_gap": 0.666526, "eig_var": 0.00040503},
}

# ================================================================
#  Derived quantities
# ================================================================

# √(135/88) candidate from Claude Opus 4.5 analysis
SQRT_135_88 = np.sqrt(135.0 / 88.0)  # 1.238584200...

# Error analysis
errors_from_sqrt32 = {}
errors_from_sqrt135_88 = {}
for m, d in crash_data.items():
    if d["crash_alpha"] is not None:
        errors_from_sqrt32[m] = {
            "absolute": d["crash_alpha"] - SQRT_3_2,
            "percent":  100 * (d["crash_alpha"] - SQRT_3_2) / SQRT_3_2,
        }
        errors_from_sqrt135_88[m] = {
            "absolute": d["crash_alpha"] - SQRT_135_88,
            "percent":  100 * (d["crash_alpha"] - SQRT_135_88) / SQRT_135_88,
        }

# ================================================================
#  Save to disk
# ================================================================

# JSON summary (human-readable)
summary = {
    "description": "Active Transport on the Prime Gas — locked-in computation results",
    "date": "2026-03-25",
    "author": "Antonio P. Matos / Fancyland LLC",
    "scripts": ["chiral_quench_test.py", "chiral_phase_transition.py", "chiral_analytics.py"],
    "constants": {
        "sqrt_3_2": float(SQRT_3_2),
        "sqrt_2_3": float(SQRT_2_3),
        "sqrt_135_88": float(SQRT_135_88),
        "product_kappa_alpha": float(SQRT_2_3 * SQRT_3_2),
    },
    "crash_data": {str(k): v for k, v in crash_data.items()},
    "ipr_data": {str(k): v for k, v in ipr_data.items()},
    "flatband_data": {str(k): v for k, v in flatband_data.items()},
    "perron_gap_data": {str(k): v for k, v in perron_gap_data.items()},
    "comm_ratio_data": {str(k): v for k, v in comm_ratio_data.items()},
    "spectral_measure": {str(k): v for k, v in spectral_measure.items()},
    "at_sqrt32": {str(k): v for k, v in at_sqrt32.items()},
    "errors_from_sqrt32": {str(k): v for k, v in errors_from_sqrt32.items()},
    "errors_from_sqrt135_88": {str(k): v for k, v in errors_from_sqrt135_88.items()},
    "falsified": [
        "sqrt(2/3) quench hypothesis",
        "chirality balance controls localization",
        "sqrt(3/2) as exact critical shear (plateau at +1.13%)",
    ],
    "confirmed": [
        "phase transition is real (derivative spike ~ -442)",
        "flat-band condensation (Var*phi -> 1.780)",
        "Perron gap universality (Delta_P -> 0.8146)",
        "1/sqrt(2) commutator ratio",
        "alpha_c converges by phi=480",
    ],
    "open_questions": [
        "closed form of alpha_c ~ 1.2386 (sqrt(135/88)?)",
        "is the plateau real or a resolution artifact?",
        "anti-Perron eigenvalue lambda_- ~ -0.591",
        "flat-band constant V_inf ~ 1.780",
    ],
}

json_path = os.path.join(SCRIPT_DIR, "active_transport_cache.json")
with open(json_path, "w") as f:
    json.dump(summary, f, indent=2)

# NPZ binary cache (for fast numpy loading)
moduli = np.array([30, 210, 2310, 30030])
phis = np.array([8, 48, 480, 5760])

crash_alphas = np.array([
    np.nan,          # m=30 (no crash)
    1.2702512563,    # m=210
    1.2385702095,    # m=2310
    1.2385680356,    # m=30030
])

steepest_alphas = np.array([
    0.5128947368,    # m=30
    1.2702261307,    # m=210
    1.2385678392,    # m=2310
    1.2385678392,    # m=30030
])

steepest_derivs = np.array([0.00, -432.61, -441.30, -442.15])

ipr_ratios = np.array([1.45, 1.88, 2.28, 2.63])

flatband_vars = np.array([0.203022, 0.036839, 0.003702, 0.000309])
flatband_products = np.array([1.624, 1.768, 1.777, 1.780])

perron_gaps = np.array([0.856826, 0.813760, 0.814621, 0.814642])

npz_path = os.path.join(SCRIPT_DIR, "active_transport_cache.npz")
np.savez_compressed(npz_path,
    moduli=moduli,
    phis=phis,
    crash_alphas=crash_alphas,
    steepest_alphas=steepest_alphas,
    steepest_derivs=steepest_derivs,
    ipr_ratios=ipr_ratios,
    flatband_vars=flatband_vars,
    flatband_products=flatband_products,
    perron_gaps=perron_gaps,
    sqrt_3_2=SQRT_3_2,
    sqrt_2_3=SQRT_2_3,
    sqrt_135_88=SQRT_135_88,
)

print(f"Saved: {json_path}")
print(f"  Size: {os.path.getsize(json_path):,} bytes")
print(f"Saved: {npz_path}")
print(f"  Size: {os.path.getsize(npz_path):,} bytes")
print()
print("Load in Python:")
print("  import numpy as np, json")
print("  cache = np.load('active_transport_cache.npz')")
print("  summary = json.load(open('active_transport_cache.json'))")
print()

# Print verification
print("=" * 60)
print("  VERIFICATION — Crash Data")
print("=" * 60)
for i, m in enumerate(moduli):
    ca = crash_alphas[i]
    if np.isnan(ca):
        print(f"  m={m:>5d}  phi={phis[i]:>5d}  crash=N/A")
    else:
        err = 100 * (ca - SQRT_3_2) / SQRT_3_2
        err135 = 100 * (ca - SQRT_135_88) / SQRT_135_88
        print(f"  m={m:>5d}  phi={phis[i]:>5d}  crash={ca:.10f}  "
              f"err(√(3/2))={err:+.4f}%  err(√(135/88))={err135:+.4f}%")

print()
print(f"  √(3/2)    = {SQRT_3_2:.12f}")
print(f"  √(135/88) = {SQRT_135_88:.12f}")
print(f"  Difference: {SQRT_135_88 - SQRT_3_2:.12f} ({100*(SQRT_135_88 - SQRT_3_2)/SQRT_3_2:+.4f}%)")
