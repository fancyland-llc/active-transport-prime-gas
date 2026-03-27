#!/usr/bin/env python3
"""
4-CLASS PARTITION ANALYSIS — RESOLVING THE SIGN FLIP
=====================================================
The nozzle_test.py showed anti-Perron asymmetry flips sign between
m=30030 (A=-0.003) and m=15015 (A=+0.002). But the mod-4 partition
covers 100% of residues at even m and only 50% at odd m.

This script does the honest comparison:
  1. Full 4-class partition {0,1,2,3 mod 4} at both moduli
  2. Parity partition {even, odd} at m=15015 (impossible at m=30030)
  3. Spatial autocorrelation of |v_r|^2 — partition-independent
  4. Entropy of eigenvector mass across mod-4 classes

If the sign flip survives the proper analysis → real topological switch.
If it's a partition artifact → we document it and move on.

Author: Claude (Anthropic) — adversarial review mode
Date: 2026-03-25
"""

import json
import numpy as np
from math import gcd
from typing import List, Dict, Tuple


def get_coprime_residues(m: int) -> List[int]:
    return sorted([r for r in range(1, m) if gcd(r, m) == 1])


def load_eigenvectors(filepath: str) -> Tuple[Dict, int, int]:
    """Load eigenvector data, return (evec_results, m, phi)."""
    with open(filepath, 'r') as f:
        data = json.load(f)
    meta = data['meta']
    m, phi = meta['m'], meta['phi']
    evec = data.get('evec_results', data.get('phase_c', {}))
    return evec, m, phi


def four_class_analysis(evec_abs: np.ndarray, residues: List[int]) -> Dict:
    """Partition eigenvector mass into all four mod-4 classes."""
    v_sq = evec_abs ** 2
    total = np.sum(v_sq)
    
    classes = {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0}
    counts = {0: 0, 1: 0, 2: 0, 3: 0}
    
    for i, r in enumerate(residues):
        c = r % 4
        classes[c] += v_sq[i]
        counts[c] += 1
    
    fractions = {c: classes[c] / total for c in range(4)}
    expected = {c: counts[c] / len(residues) for c in range(4)}
    
    # Enrichment: fraction of mass / fraction of residues
    enrichment = {}
    for c in range(4):
        if expected[c] > 0:
            enrichment[c] = fractions[c] / expected[c]
        else:
            enrichment[c] = float('nan')
    
    return {
        'counts': counts,
        'mass': {c: classes[c] for c in range(4)},
        'fractions': fractions,
        'expected': expected,
        'enrichment': enrichment
    }


def spatial_autocorrelation(evec_abs: np.ndarray, residues: List[int], m: int) -> Dict:
    """Compute spatial autocorrelation of |v_r|^2 along the residue ring.
    
    This is partition-independent: it measures how much nearby residues
    on Z/mZ have correlated eigenvector weight.
    """
    v_sq = evec_abs ** 2
    v_centered = v_sq - np.mean(v_sq)
    var = np.var(v_sq)
    
    if var < 1e-30:
        return {'lag1': 0.0, 'lag2': 0.0, 'interpretation': 'flat eigenvector'}
    
    # Sort by residue position on the circle
    order = np.argsort(residues)
    v_ordered = v_centered[order]
    n = len(v_ordered)
    
    # Autocorrelation at lag 1, 2, 5
    results = {}
    for lag in [1, 2, 5, 10]:
        if lag < n:
            corr = np.sum(v_ordered[:-lag] * v_ordered[lag:]) / ((n - lag) * var)
            results[f'lag{lag}'] = corr
    
    return results


def parity_analysis(evec_abs: np.ndarray, residues: List[int]) -> Dict:
    """For odd m: partition into even vs odd residues."""
    v_sq = evec_abs ** 2
    total = np.sum(v_sq)
    
    even_mass = sum(v_sq[i] for i, r in enumerate(residues) if r % 2 == 0)
    odd_mass = sum(v_sq[i] for i, r in enumerate(residues) if r % 2 == 1)
    even_count = sum(1 for r in residues if r % 2 == 0)
    odd_count = sum(1 for r in residues if r % 2 == 1)
    
    return {
        'even_mass_frac': even_mass / total,
        'odd_mass_frac': odd_mass / total,
        'even_count_frac': even_count / len(residues),
        'odd_count_frac': odd_count / len(residues),
        'even_enrichment': (even_mass / total) / (even_count / len(residues)) if even_count > 0 else float('nan'),
        'odd_enrichment': (odd_mass / total) / (odd_count / len(residues)) if odd_count > 0 else float('nan'),
        'parity_asymmetry': (odd_mass - even_mass) / total
    }


def run_analysis():
    print("=" * 72)
    print("  4-CLASS PARTITION ANALYSIS — RESOLVING THE SIGN FLIP")
    print("=" * 72)
    print()
    
    # Load both datasets
    datasets = {}
    for label, path in [
        ('m=30030', 'overnight_results/job3_m30030_micro.json'),
        ('m=15015', 'overnight_results/job4_m15015.json'),
    ]:
        try:
            evec, m, phi = load_eigenvectors(path)
            residues = get_coprime_residues(m)
            datasets[label] = {'evec': evec, 'm': m, 'phi': phi, 'residues': residues}
            print(f"  Loaded {label}: m={m}, φ={phi}, {len(evec)} α values")
        except FileNotFoundError:
            print(f"  WARNING: {path} not found, skipping {label}")
    
    print()
    
    # =========================================================================
    # ANALYSIS 1: Full 4-class partition
    # =========================================================================
    
    print("=" * 72)
    print("  ANALYSIS 1: FULL 4-CLASS PARTITION OF ANTI-PERRON EIGENVECTOR")
    print("=" * 72)
    print()
    
    for label, ds in datasets.items():
        m = ds['m']
        residues = ds['residues']
        evec_data = ds['evec']
        
        print(f"  --- {label} (m={m}, {'EVEN' if m % 2 == 0 else 'ODD'}) ---")
        
        # Aggregate across all α values
        all_results = []
        for alpha_key in sorted(evec_data.keys(), key=lambda x: float(x)):
            entry = evec_data[alpha_key]
            if 'anti_perron_evec_abs' in entry:
                ap = np.array(entry['anti_perron_evec_abs'])
                result = four_class_analysis(ap, residues)
                all_results.append(result)
        
        if not all_results:
            print("    No anti-Perron eigenvector data!\n")
            continue
        
        # Average across α values
        print(f"    Averaged over {len(all_results)} α values:")
        print()
        print(f"    {'Class':>8} {'Count':>8} {'Expected%':>10} {'Mass%':>10} {'Enrichment':>12}")
        print(f"    {'-'*8:>8} {'-'*8:>8} {'-'*10:>10} {'-'*10:>10} {'-'*12:>12}")
        
        for c in range(4):
            count = all_results[0]['counts'][c]
            expected = all_results[0]['expected'][c]
            mass_frac = np.mean([r['fractions'][c] for r in all_results])
            enrichment = np.mean([r['enrichment'][c] for r in all_results])
            
            marker = ""
            if count == 0:
                marker = " (empty)"
            elif abs(enrichment - 1.0) > 0.001:
                marker = f" {'↑' if enrichment > 1.0 else '↓'}"
            
            print(f"    r≡{c}(4) {count:>8} {100*expected:>9.1f}% {100*mass_frac:>9.4f}% {enrichment:>11.6f}{marker}")
        
        # The CORRECT asymmetry: within the {1,3} sector only
        A_13 = np.mean([r['fractions'][1] - r['fractions'][3] for r in all_results])
        # Old-style asymmetry (normalized by total, as in nozzle_test.py)
        A_old = A_13  # same formula, just being explicit
        # Sector-normalized: within {1,3} only
        sector_13_mass = np.mean([r['fractions'][1] + r['fractions'][3] for r in all_results])
        A_sector = A_13 / sector_13_mass if sector_13_mass > 1e-15 else 0.0
        
        print()
        print(f"    A (total-normalized, nozzle_test.py):  {A_old:+.10f}")
        print(f"    A (sector-normalized, within {{1,3}}):  {A_sector:+.10f}")
        print(f"    Sector {{1,3}} carries {100*sector_13_mass:.2f}% of total mass")
        print()
        
        # For odd m: also show the {0,2} sector
        if m % 2 == 1:
            A_02 = np.mean([r['fractions'][0] - r['fractions'][2] for r in all_results])
            sector_02_mass = np.mean([r['fractions'][0] + r['fractions'][2] for r in all_results])
            A_02_sector = A_02 / sector_02_mass if sector_02_mass > 1e-15 else 0.0
            print(f"    A in {{0,2}} sector (total-normalized):   {A_02:+.10f}")
            print(f"    A in {{0,2}} sector (sector-normalized):  {A_02_sector:+.10f}")
            print(f"    Sector {{0,2}} carries {100*sector_02_mass:.2f}% of total mass")
            print()
    
    # =========================================================================
    # ANALYSIS 2: Parity partition (even vs odd residues) — odd m only
    # =========================================================================
    
    print("=" * 72)
    print("  ANALYSIS 2: PARITY PARTITION (ODD m ONLY)")
    print("=" * 72)
    print()
    
    if 'm=15015' in datasets:
        ds = datasets['m=15015']
        m = ds['m']
        residues = ds['residues']
        evec_data = ds['evec']
        
        print(f"  m={m} (ODD): even vs odd residue eigenvector mass")
        print()
        
        # Anti-Perron
        for mode_name, mode_key in [('Anti-Perron', 'anti_perron_evec_abs'), ('Perron', 'perron_evec_abs')]:
            parity_results = []
            for alpha_key in sorted(evec_data.keys(), key=lambda x: float(x)):
                entry = evec_data[alpha_key]
                if mode_key in entry:
                    evec = np.array(entry[mode_key])
                    parity_results.append(parity_analysis(evec, residues))
            
            if parity_results:
                print(f"  {mode_name} eigenvector:")
                mean_even = np.mean([r['even_mass_frac'] for r in parity_results])
                mean_odd = np.mean([r['odd_mass_frac'] for r in parity_results])
                mean_even_enr = np.mean([r['even_enrichment'] for r in parity_results])
                mean_odd_enr = np.mean([r['odd_enrichment'] for r in parity_results])
                asym = np.mean([r['parity_asymmetry'] for r in parity_results])
                
                print(f"    Even residues: {100*mean_even:.4f}% mass (expected 50.0%), enrichment = {mean_even_enr:.6f}")
                print(f"    Odd residues:  {100*mean_odd:.4f}% mass (expected 50.0%), enrichment = {mean_odd_enr:.6f}")
                print(f"    Parity asymmetry (odd - even): {asym:+.10f}")
                print()
    
    # =========================================================================
    # ANALYSIS 3: Spatial autocorrelation — partition-independent
    # =========================================================================
    
    print("=" * 72)
    print("  ANALYSIS 3: SPATIAL AUTOCORRELATION (PARTITION-INDEPENDENT)")
    print("=" * 72)
    print()
    
    for label, ds in datasets.items():
        m = ds['m']
        residues = ds['residues']
        evec_data = ds['evec']
        
        print(f"  --- {label} Anti-Perron eigenvector ---")
        
        all_corrs = []
        for alpha_key in sorted(evec_data.keys(), key=lambda x: float(x)):
            entry = evec_data[alpha_key]
            if 'anti_perron_evec_abs' in entry:
                ap = np.array(entry['anti_perron_evec_abs'])
                corr = spatial_autocorrelation(ap, residues, m)
                all_corrs.append(corr)
        
        if all_corrs:
            for lag_key in ['lag1', 'lag2', 'lag5', 'lag10']:
                if lag_key in all_corrs[0]:
                    mean_c = np.mean([c[lag_key] for c in all_corrs])
                    print(f"    Autocorrelation at {lag_key}: {mean_c:+.6f}")
        print()
    
    # =========================================================================
    # ANALYSIS 4: Perron eigenvector (should be flat at both)
    # =========================================================================
    
    print("=" * 72)
    print("  ANALYSIS 4: PERRON EIGENVECTOR 4-CLASS (CONTROL)")
    print("=" * 72)
    print()
    
    for label, ds in datasets.items():
        m = ds['m']
        residues = ds['residues']
        evec_data = ds['evec']
        
        print(f"  --- {label} Perron eigenvector ---")
        
        all_results = []
        for alpha_key in sorted(evec_data.keys(), key=lambda x: float(x)):
            entry = evec_data[alpha_key]
            if 'perron_evec_abs' in entry:
                p = np.array(entry['perron_evec_abs'])
                result = four_class_analysis(p, residues)
                all_results.append(result)
        
        if all_results:
            for c in range(4):
                count = all_results[0]['counts'][c]
                if count == 0:
                    continue
                enrichment = np.mean([r['enrichment'][c] for r in all_results])
                print(f"    r≡{c}(4): enrichment = {enrichment:.8f} (expect 1.0)")
        print()
    
    # =========================================================================
    # VERDICT
    # =========================================================================
    
    print("=" * 72)
    print("  VERDICT")
    print("=" * 72)
    print()
    
    if 'm=30030' in datasets and 'm=15015' in datasets:
        # Get sector-normalized asymmetries
        for label in ['m=30030', 'm=15015']:
            ds = datasets[label]
            residues = ds['residues']
            evec_data = ds['evec']
            
            results = []
            for alpha_key in sorted(evec_data.keys(), key=lambda x: float(x)):
                entry = evec_data[alpha_key]
                if 'anti_perron_evec_abs' in entry:
                    ap = np.array(entry['anti_perron_evec_abs'])
                    r = four_class_analysis(ap, residues)
                    results.append(r)
            
            if results:
                A_total = np.mean([r['fractions'][1] - r['fractions'][3] for r in results])
                sector = np.mean([r['fractions'][1] + r['fractions'][3] for r in results])
                A_sect = A_total / sector if sector > 0 else 0
                print(f"  {label}:")
                print(f"    A (total-normalized): {A_total:+.8f}")
                print(f"    A (sector-normalized): {A_sect:+.8f}")
                print(f"    {1,3}-sector coverage: {100*sector:.1f}%")
                print()
        
        print("  If sector-normalized A has SAME SIGN → sign flip is partition artifact")
        print("  If sector-normalized A has DIFFERENT SIGN → genuine topological switch")
        print()


if __name__ == "__main__":
    run_analysis()
