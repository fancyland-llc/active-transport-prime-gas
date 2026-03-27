#!/usr/bin/env python3
"""
M=15015 ANTI-PERRON EXTRACTION — THE NOZZLE TEST
==================================================
Tests Gemini's hypothesis: Does m=15015 (odd modulus, no factor of 2)
have chiral asymmetry A ≈ 0?

If A(15015) = 0 but A(30030) = -0.003:
  → The factor of 2 IS the nozzle
  → Primorial geometry is MANDATORY for directional thrust
  → m=15015 is a thermal bomb, not a rocket

If A(15015) ≠ 0:
  → The asymmetry comes from something else
  → The factor of 2 only affects sharpness, not direction

Author: Claude + Gemini + Tony (Fancyland LLC)
Date: 2026-03-25
"""

import json
import numpy as np
from math import gcd
from typing import List, Dict

def get_coprime_residues(m: int) -> List[int]:
    """Return residues coprime to m in sorted order."""
    return sorted([r for r in range(1, m) if gcd(r, m) == 1])

def compute_chiral_asymmetry(evec_abs: np.ndarray, residues: List[int]) -> Dict:
    """Compute chiral asymmetry A = (odd - even) / total."""
    v_sq = evec_abs ** 2
    total = np.sum(v_sq)
    
    odd_sum = 0.0   # r ≡ 1 (mod 4)
    even_sum = 0.0  # r ≡ 3 (mod 4)
    
    for i, r in enumerate(residues):
        if r % 4 == 1:
            odd_sum += v_sq[i]
        elif r % 4 == 3:
            even_sum += v_sq[i]
    
    asymmetry = (odd_sum - even_sum) / total if total > 1e-15 else 0.0
    
    return {
        'asymmetry': asymmetry,
        'odd_fraction': odd_sum / total,
        'even_fraction': even_sum / total
    }

def analyze_m15015(filepath: str):
    """Full analysis of m=15015 eigenvector data."""
    
    print("=" * 72)
    print("  M=15015 ANTI-PERRON EXTRACTION — THE NOZZLE TEST")
    print("  Is m=15015 a rocket or a bomb?")
    print("=" * 72)
    print()
    
    # Load data
    with open(filepath, 'r') as f:
        data = json.load(f)
    
    meta = data['meta']
    m = meta['m']
    phi = meta['phi']
    
    print(f"  m = {m}")
    print(f"  φ = {phi}")
    print(f"  Primes: {meta.get('primes', 'N/A')}")
    print(f"  NOTE: m is {'ODD' if m % 2 == 1 else 'EVEN'} — factor of 2 is {'ABSENT' if m % 2 == 1 else 'PRESENT'}")
    print()
    
    # Get residues
    residues = get_coprime_residues(m)
    
    # Check for eigenvector data
    if 'evec_results' in data:
        evec_results = data['evec_results']
    elif 'phase_c' in data:
        evec_results = data['phase_c']
    else:
        print("  ERROR: No eigenvector data found!")
        print(f"  Available keys: {list(data.keys())}")
        return
    
    print(f"  Eigenvector data at {len(evec_results)} α values")
    print()
    
    # =========================================================================
    # EXTRACTION 1: Eigenvalue Gap
    # =========================================================================
    
    print("-" * 72)
    print("  EXTRACTION 1: EIGENVALUE GAP Δλ")
    print("-" * 72)
    print()
    
    gaps_15015 = []
    lambda_P_values = []
    lambda_minus_values = []
    
    for alpha_key in sorted(evec_results.keys(), key=lambda x: float(x)):
        alpha = float(alpha_key)
        entry = evec_results[alpha_key]
        
        if 'eigenvalues' in entry:
            eigenvalues = np.array(entry['eigenvalues'])
            lambda_P = np.max(eigenvalues)
            lambda_minus = np.min(eigenvalues)
            gap = lambda_P - lambda_minus
            
            gaps_15015.append(gap)
            lambda_P_values.append(lambda_P)
            lambda_minus_values.append(lambda_minus)
            
            print(f"  α = {alpha:.8f}: λ_P = {lambda_P:+.8f}, λ_- = {lambda_minus:+.8f}, Δλ = {gap:.8f}")
    
    if gaps_15015:
        mean_gap = np.mean(gaps_15015)
        print()
        print(f"  Mean gap (m=15015): Δλ = {mean_gap:.10f}")
        print(f"  Gap (m=30030):      Δλ = 1.7441353394")
        print(f"  Difference:         {mean_gap - 1.7441353394:+.10f} ({100*(mean_gap - 1.7441353394)/1.7441353394:+.4f}%)")
        
        if abs(mean_gap - 1.7441353394) < 0.001:
            print()
            print("  → Gap is IDENTICAL (within 0.06%)")
            print("  → Resonance frequency is dimension-determined, not 2-torsion")
        else:
            print()
            print("  → Gap SHIFTED!")
            print("  → The factor of 2 tunes the resonance frequency")
    
    # =========================================================================
    # EXTRACTION 2: Anti-Perron Asymmetry
    # =========================================================================
    
    print()
    print("-" * 72)
    print("  EXTRACTION 2: ANTI-PERRON CHIRAL ASYMMETRY")
    print("-" * 72)
    print()
    
    anti_perron_asymmetries = []
    perron_asymmetries = []
    
    for alpha_key in sorted(evec_results.keys(), key=lambda x: float(x)):
        alpha = float(alpha_key)
        entry = evec_results[alpha_key]
        
        # Check for eigenvector data
        if 'anti_perron_evec_abs' in entry:
            anti_perron = np.array(entry['anti_perron_evec_abs'])
            result = compute_chiral_asymmetry(anti_perron, residues)
            anti_perron_asymmetries.append(result['asymmetry'])
            
        if 'perron_evec_abs' in entry:
            perron = np.array(entry['perron_evec_abs'])
            result = compute_chiral_asymmetry(perron, residues)
            perron_asymmetries.append(result['asymmetry'])
    
    if anti_perron_asymmetries:
        mean_anti_A = np.mean(anti_perron_asymmetries)
        print(f"  Anti-Perron asymmetry (m=15015): A = {mean_anti_A:+.10f}")
        print(f"  Anti-Perron asymmetry (m=30030): A = -0.0029603701")
        print()
        
        if abs(mean_anti_A) < 0.0001:
            print("  ╔══════════════════════════════════════════════════════════════╗")
            print("  ║  ★ GEMINI'S HYPOTHESIS CONFIRMED ★                           ║")
            print("  ║                                                              ║")
            print("  ║  m=15015 asymmetry A ≈ 0                                     ║")
            print("  ║  m=30030 asymmetry A = -0.003                                ║")
            print("  ║                                                              ║")
            print("  ║  THE FACTOR OF 2 IS THE NOZZLE.                              ║")
            print("  ║  Without it, exhaust is SYMMETRIC → thermal bomb.            ║")
            print("  ║  With it, exhaust is DIRECTIONAL → rocket engine.            ║")
            print("  ║                                                              ║")
            print("  ║  PRIMORIAL GEOMETRY IS MANDATORY FOR THRUST.                 ║")
            print("  ╚══════════════════════════════════════════════════════════════╝")
        else:
            print("  ╔══════════════════════════════════════════════════════════════╗")
            print("  ║  GEMINI'S HYPOTHESIS FALSIFIED                               ║")
            print("  ║                                                              ║")
            print(f"  ║  m=15015 asymmetry A = {mean_anti_A:+.8f}                      ║")
            print("  ║  The asymmetry survives without the factor of 2.             ║")
            print("  ║  The nozzle mechanism is NOT 2-dependent.                    ║")
            print("  ╚══════════════════════════════════════════════════════════════╝")
    else:
        print("  No Anti-Perron eigenvector data found!")
        print("  Checking for alternative formats...")
        
        # Try to find any eigenvector data
        sample_key = list(evec_results.keys())[0]
        sample_entry = evec_results[sample_key]
        print(f"  Available fields in eigenvector entry: {list(sample_entry.keys())}")
    
    if perron_asymmetries:
        mean_perron_A = np.mean(perron_asymmetries)
        print()
        print(f"  Perron asymmetry (m=15015):    A = {mean_perron_A:+.10f}")
        print(f"  Perron asymmetry (m=30030):    A = -0.0000024600")
    
    # =========================================================================
    # EXTRACTION 3: Mod-4 Partition Analysis
    # =========================================================================
    
    print()
    print("-" * 72)
    print("  EXTRACTION 3: MOD-4 PARTITION STRUCTURE")
    print("-" * 72)
    print()
    
    # For odd m, ALL residues are odd (coprime to m means coprime to all odd primes)
    # So all residues are either 1 or 3 mod 4
    
    odd_count = sum(1 for r in residues if r % 4 == 1)
    even_count = sum(1 for r in residues if r % 4 == 3)
    other_count = sum(1 for r in residues if r % 4 not in [1, 3])
    
    print(f"  m = {m} ({'ODD' if m % 2 == 1 else 'EVEN'})")
    print(f"  Total coprime residues: {len(residues)}")
    print(f"  r ≡ 1 (mod 4): {odd_count} ({100*odd_count/len(residues):.1f}%)")
    print(f"  r ≡ 3 (mod 4): {even_count} ({100*even_count/len(residues):.1f}%)")
    print(f"  r ≡ 0,2 (mod 4): {other_count}")
    
    # The KEY observation: for ODD m, we still have the mod-4 partition
    # But the multiplicative structure is different
    
    print()
    print("  For m=30030 (EVEN):")
    print("    All coprime residues are ODD (can't be divisible by 2)")
    print("    Partition by r mod 4 gives 1 vs 3")
    print()
    print("  For m=15015 (ODD):")
    print("    Residues can be EVEN (coprime to odd m)")
    print("    r ≡ 0 or 2 (mod 4) are now allowed!")
    
    even_residues = [r for r in residues if r % 2 == 0]
    print(f"    Even residues in m=15015: {len(even_residues)} ({100*len(even_residues)/len(residues):.1f}%)")
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    
    print()
    print("=" * 72)
    print("  VERDICT: ROCKET OR BOMB?")
    print("=" * 72)
    print()
    
    if anti_perron_asymmetries:
        A_15015 = np.mean(anti_perron_asymmetries)
        A_30030 = -0.0029603701
        
        if abs(A_15015) < 0.0001:
            print("  m=15015: A ≈ 0 → THERMAL BOMB (symmetric exhaust)")
            print("  m=30030: A = -0.003 → DIRECTED ROCKET (asymmetric exhaust)")
            print()
            print("  ★ THE FACTOR OF 2 IS THE ENGINE NOZZLE ★")
            print()
            print("  Engineering implication:")
            print("    TF coils MUST use primorial geometry")
            print("    Odd modulus = explosion, not propulsion")
        elif abs(A_15015 - A_30030) < 0.001:
            print(f"  m=15015: A = {A_15015:+.6f}")
            print(f"  m=30030: A = {A_30030:+.6f}")
            print("  → Asymmetries are IDENTICAL")
            print()
            print("  The nozzle mechanism is NOT 2-dependent.")
            print("  The factor of 2 only affects transition sharpness.")
        else:
            print(f"  m=15015: A = {A_15015:+.6f}")
            print(f"  m=30030: A = {A_30030:+.6f}")
            print(f"  → Asymmetry SHIFTED by {A_15015 - A_30030:+.6f}")
            print()
            print("  The factor of 2 MODULATES the nozzle direction.")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        filepath = sys.argv[1]
    else:
        # Try default locations
        paths = [
            'overnight_results/job4_m15015.json',
        ]
        filepath = None
        for p in paths:
            try:
                with open(p, 'r') as f:
                    filepath = p
                    break
            except:
                continue
        
        if filepath is None:
            print("ERROR: job4_m15015.json not found!")
            print("Usage: python nozzle_test.py <path_to_job4_m15015.json>")
            sys.exit(1)
    
    analyze_m15015(filepath)