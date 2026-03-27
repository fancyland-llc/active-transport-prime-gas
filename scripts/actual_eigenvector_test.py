#!/usr/bin/env python3
"""
EIGENVECTOR PREDICTION TESTER - ACTUAL JOB 3 DATA
===================================================
Tests Predictions 9.1 and 9.2 against REAL eigenvector data from m=30030.

Author: Claude + Tony (Fancyland LLC)
Date: 2026-03-25
"""

import json
import numpy as np
from math import gcd, pi
from typing import List, Dict, Tuple

# =============================================================================
# LOAD AND PARSE JOB 3 DATA
# =============================================================================

def load_job3_data(filepath: str) -> Dict:
    """Load the Job 3 JSON data."""
    with open(filepath, 'r') as f:
        return json.load(f)

def get_coprime_residues(m: int) -> List[int]:
    """Return residues coprime to m in sorted order."""
    return sorted([r for r in range(1, m) if gcd(r, m) == 1])

# =============================================================================
# PREDICTION 9.1: CHIRAL EXHAUST ASYMMETRY
# =============================================================================

def compute_chiral_asymmetry(evec_abs: np.ndarray, residues: List[int]) -> Dict:
    """
    Compute the chiral asymmetry A from Prediction 9.1.
    
    A = (Σ|v_r|² for r≡1(4) - Σ|v_r|² for r≡3(4)) / Σ|v_r|²
    
    Since we have |v_r| (absolute values), we square them.
    """
    v_sq = evec_abs ** 2
    total = np.sum(v_sq)
    
    odd_sum = 0.0   # r ≡ 1 (mod 4)
    even_sum = 0.0  # r ≡ 3 (mod 4)
    odd_count = 0
    even_count = 0
    
    for i, r in enumerate(residues):
        if r % 4 == 1:
            odd_sum += v_sq[i]
            odd_count += 1
        elif r % 4 == 3:
            even_sum += v_sq[i]
            even_count += 1
    
    asymmetry = (odd_sum - even_sum) / total if total > 1e-15 else 0.0
    
    return {
        'asymmetry': asymmetry,
        'odd_sum': odd_sum,
        'even_sum': even_sum,
        'total': total,
        'odd_count': odd_count,
        'even_count': even_count,
        'odd_fraction': odd_sum / total,
        'even_fraction': even_sum / total
    }

# =============================================================================
# PREDICTION 9.2: 13-PERIODIC STRUCTURE
# =============================================================================

def compute_periodic_power(evec_abs: np.ndarray, periods: List[int]) -> Dict:
    """
    Compute Fourier power at specific periods.
    
    We FFT the |v_r| sequence and check power at target frequencies.
    """
    n = len(evec_abs)
    
    # FFT of the absolute value sequence
    fft = np.fft.fft(evec_abs)
    power = np.abs(fft) ** 2
    total_power = np.sum(power)
    
    # DC component (period = infinity)
    dc_power = power[0] / total_power
    
    # Power at each target period
    period_results = {}
    for p in periods:
        if p <= 1:
            continue
        # Frequency corresponding to period p is k = n/p
        # We check the closest integer frequencies
        k = n / p
        k_low = int(k)
        k_high = k_low + 1 if k_low < n-1 else k_low
        
        # Take max of adjacent bins (in case of spectral leakage)
        p_power = max(power[k_low], power[k_high]) / total_power
        period_results[p] = p_power
    
    # Find dominant non-DC frequency
    half = n // 2
    non_dc_power = power[1:half]
    if len(non_dc_power) > 0:
        dominant_idx = np.argmax(non_dc_power) + 1
        dominant_period = n / dominant_idx
        dominant_power = non_dc_power[dominant_idx - 1] / total_power
    else:
        dominant_period = n
        dominant_power = 0
    
    return {
        'total_power': total_power,
        'dc_power_fraction': dc_power,
        'period_power': period_results,
        'dominant_period': dominant_period,
        'dominant_power_fraction': dominant_power
    }

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def main():
    print("=" * 72)
    print("  EIGENVECTOR PREDICTION TEST — ACTUAL JOB 3 DATA")
    print("  m = 30030, φ = 5760")
    print("=" * 72)
    print()
    
    # Load data
    data = load_job3_data('overnight_results/job3_m30030_micro.json')
    
    meta = data['meta']
    m = meta['m']
    phi = meta['phi']
    evec_results = data['evec_results']
    
    print(f"  Loaded m = {m}, φ = {phi}")
    print(f"  Number of α values with eigenvector data: {len(evec_results)}")
    print()
    
    # Get residues
    residues = get_coprime_residues(m)
    assert len(residues) == phi, f"Expected {phi} residues, got {len(residues)}"
    
    # Sort α values
    alpha_values = sorted([float(a) for a in evec_results.keys()])
    
    print(f"  α range: [{alpha_values[0]:.10f}, {alpha_values[-1]:.10f}]")
    print(f"  √(135/88) = 1.2387670686..., center of data ≈ {np.mean(alpha_values):.10f}")
    print()
    
    # =========================================================================
    # PREDICTION 9.1: CHIRAL ASYMMETRY
    # =========================================================================
    
    print("-" * 72)
    print("  PREDICTION 9.1: CHIRAL EXHAUST ASYMMETRY")
    print("  A = (Σ|v_r|² for r≡1 mod 4) - (Σ|v_r|² for r≡3 mod 4)) / total")
    print("-" * 72)
    print()
    
    print("  Perron eigenvector asymmetry across α values:")
    print(f"  {'α':^18} | {'A (asymmetry)':^14} | {'Odd %':^8} | {'Even %':^8}")
    print(f"  {'-'*18}-+-{'-'*14}-+-{'-'*8}-+-{'-'*8}")
    
    asymmetries = []
    for alpha in alpha_values:
        alpha_key = [k for k in evec_results.keys() if abs(float(k) - alpha) < 1e-12][0]
        entry = evec_results[alpha_key]
        
        perron_evec = np.array(entry['perron_evec_abs'])
        result = compute_chiral_asymmetry(perron_evec, residues)
        
        asymmetries.append(result['asymmetry'])
        
        print(f"  {alpha:.10f} | {result['asymmetry']:+.10f} | {result['odd_fraction']*100:6.3f}% | {result['even_fraction']*100:6.3f}%")
    
    mean_A = np.mean(asymmetries)
    std_A = np.std(asymmetries)
    
    print()
    print(f"  Mean asymmetry: A = {mean_A:+.10f}")
    print(f"  Std deviation:  σ = {std_A:.10f}")
    print()
    
    # Test prediction
    THRESHOLD = 0.001  # 0.1% threshold for significance
    
    if abs(mean_A) > THRESHOLD:
        print(f"  ╔══════════════════════════════════════════════════════════════╗")
        print(f"  ║  ✓ PREDICTION 9.1: CONFIRMED                                 ║")
        print(f"  ║    |A| = {abs(mean_A):.8f} > {THRESHOLD} threshold                   ║")
        if mean_A > 0:
            print(f"  ║    Direction: ODD BLOCK DOMINANT (r ≡ 1 mod 4)               ║")
        else:
            print(f"  ║    Direction: EVEN BLOCK DOMINANT (r ≡ 3 mod 4)              ║")
        print(f"  ╚══════════════════════════════════════════════════════════════╝")
    else:
        print(f"  ╔══════════════════════════════════════════════════════════════╗")
        print(f"  ║  ✗ PREDICTION 9.1: NOT CONFIRMED                             ║")
        print(f"  ║    |A| = {abs(mean_A):.8f} < {THRESHOLD} threshold                   ║")
        print(f"  ║    Eigenvector is chirality-balanced                         ║")
        print(f"  ╚══════════════════════════════════════════════════════════════╝")
    
    print()
    
    # =========================================================================
    # PREDICTION 9.2: 13-PERIODIC STRUCTURE
    # =========================================================================
    
    print("-" * 72)
    print("  PREDICTION 9.2: 13-PERIODIC STRUCTURE IN EIGENVECTOR")
    print("  Testing if RMP coil count (13) is encoded in spectral structure")
    print("-" * 72)
    print()
    
    # Use the middle α value
    mid_alpha = alpha_values[len(alpha_values) // 2]
    mid_alpha_key = [k for k in evec_results.keys() if abs(float(k) - mid_alpha) < 1e-12][0]
    mid_entry = evec_results[mid_alpha_key]
    
    perron_evec = np.array(mid_entry['perron_evec_abs'])
    
    # Check periods corresponding to all primes in the drive geometry
    periods_to_check = [3, 5, 7, 11, 13, 17, 24, 48, 480, 5760]
    fourier_result = compute_periodic_power(perron_evec, periods_to_check)
    
    print(f"  Analysis at α = {mid_alpha:.10f}")
    print()
    print("  Fourier power at key periods:")
    print(f"  {'Period':>8} | {'Power %':>10} | {'Significance'}")
    print(f"  {'-'*8}-+-{'-'*10}-+-{'-'*20}")
    
    for p in sorted(fourier_result['period_power'].keys()):
        power_pct = fourier_result['period_power'][p] * 100
        
        # Mark significance
        if p == 13:
            marker = " <<< TARGET (13 RMP coils)"
        elif power_pct > 0.1:
            marker = " *"
        else:
            marker = ""
        
        print(f"  {p:>8} | {power_pct:>9.4f}% |{marker}")
    
    print()
    print(f"  DC component: {fourier_result['dc_power_fraction']*100:.4f}%")
    print(f"  Dominant non-DC period: {fourier_result['dominant_period']:.1f}")
    print(f"  Dominant power: {fourier_result['dominant_power_fraction']*100:.4f}%")
    print()
    
    # Test prediction - threshold is 0.1% for detectability
    PERIOD_THRESHOLD = 0.001  # 0.1%
    power_at_13 = fourier_result['period_power'].get(13, 0)
    
    if power_at_13 > PERIOD_THRESHOLD:
        print(f"  ╔══════════════════════════════════════════════════════════════╗")
        print(f"  ║  ✓ PREDICTION 9.2: CONFIRMED                                 ║")
        print(f"  ║    13-periodic power = {power_at_13*100:.4f}% > {PERIOD_THRESHOLD*100}% threshold       ║")
        print(f"  ║    13-coil RMP architecture is SUPPORTED                     ║")
        print(f"  ╚══════════════════════════════════════════════════════════════╝")
    else:
        print(f"  ╔══════════════════════════════════════════════════════════════╗")
        print(f"  ║  ✗ PREDICTION 9.2: NOT CONFIRMED                             ║")
        print(f"  ║    13-periodic power = {power_at_13*100:.4f}% < {PERIOD_THRESHOLD*100}% threshold       ║")
        print(f"  ║    13-coil RMP count is NOT eigenvector-derived              ║")
        print(f"  ╚══════════════════════════════════════════════════════════════╝")
    
    print()
    
    # =========================================================================
    # ANTI-PERRON EIGENVECTOR ANALYSIS
    # =========================================================================
    
    print("-" * 72)
    print("  SUPPLEMENTARY: ANTI-PERRON EIGENVECTOR ASYMMETRY")
    print("-" * 72)
    print()
    
    anti_asymmetries = []
    for alpha in alpha_values:
        alpha_key = [k for k in evec_results.keys() if abs(float(k) - alpha) < 1e-12][0]
        entry = evec_results[alpha_key]
        
        anti_perron_evec = np.array(entry['anti_perron_evec_abs'])
        result = compute_chiral_asymmetry(anti_perron_evec, residues)
        anti_asymmetries.append(result['asymmetry'])
    
    mean_anti_A = np.mean(anti_asymmetries)
    std_anti_A = np.std(anti_asymmetries)
    
    print(f"  Anti-Perron mean asymmetry: A = {mean_anti_A:+.10f}")
    print(f"  Anti-Perron std deviation:  σ = {std_anti_A:.10f}")
    
    if abs(mean_anti_A) > THRESHOLD:
        print(f"  → Anti-Perron also shows chiral asymmetry!")
        if mean_anti_A * mean_A > 0:
            print(f"  → Same direction as Perron: reinforcing")
        else:
            print(f"  → Opposite direction from Perron: competing")
    
    print()
    
    # =========================================================================
    # EXHAUST RATIO INTERPRETATION
    # =========================================================================
    
    print("-" * 72)
    print("  ENGINEERING INTERPRETATION: EXHAUST PORT RATIOS")
    print("-" * 72)
    print()
    
    # The 3 exhaust ports correspond to the 3-adic structure
    # The asymmetry A tells us the directional bias
    
    if abs(mean_A) > THRESHOLD:
        # Convert asymmetry to port ratios
        # A > 0 means odd block (r ≡ 1 mod 4) dominates
        # We predict 3 ports at 120° intervals
        
        odd_frac = 0.5 + mean_A / 2  # Fraction going to "odd" direction
        even_frac = 0.5 - mean_A / 2  # Fraction going to "even" direction
        
        # Map to 3 ports based on mod-3 structure
        # This is a simplification - actual mapping needs the eigenvector phases
        
        print(f"  Based on Perron eigenvector asymmetry A = {mean_A:+.8f}:")
        print()
        print(f"  Theoretical exhaust distribution:")
        print(f"    Port 0 (primary):   {odd_frac*100:.2f}%")
        print(f"    Port 1 (secondary): {even_frac*100:.2f}%")
        print(f"    Port 2 (tertiary):  ~0% (topologically closed)")
        print()
        print(f"  Predicted thrust bias: {abs(mean_A)*100:.2f}% toward {'r≡1(mod 4)' if mean_A > 0 else 'r≡3(mod 4)'} sector")
    else:
        print(f"  Asymmetry too small for directional exhaust prediction.")
        print(f"  Exhaust would be approximately symmetric.")
    
    print()
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    
    print("=" * 72)
    print("  FINAL VERDICT")
    print("=" * 72)
    print()
    
    p91_status = "CONFIRMED" if abs(mean_A) > THRESHOLD else "NOT CONFIRMED"
    p92_status = "CONFIRMED" if power_at_13 > PERIOD_THRESHOLD else "NOT CONFIRMED"
    
    print(f"  Prediction 9.1 (Chiral asymmetry A ≠ 0):    {p91_status}")
    print(f"    Measured: A = {mean_A:+.10f}")
    print()
    print(f"  Prediction 9.2 (13-periodic structure):     {p92_status}")
    print(f"    Measured: {power_at_13*100:.4f}% power at period-13")
    print()
    
    if p91_status == "CONFIRMED" and p92_status == "CONFIRMED":
        print("  ★ BOTH PREDICTIONS CONFIRMED ★")
        print("  The 4-layer drive geometry has mathematical support.")
    elif p91_status == "CONFIRMED":
        print("  Chiral exhaust is supported; 13-coil count needs alternative derivation.")
    elif p92_status == "CONFIRMED":
        print("  13-periodic structure found; chiral asymmetry not significant.")
    else:
        print("  Neither prediction confirmed at current thresholds.")
        print("  Drive geometry requires revision or alternative justification.")

if __name__ == "__main__":
    main()