#!/usr/bin/env python3
"""
ADDITIVE GHOST EXTRACTOR
=========================
Isolates the pure prime-gap defects from the D(m) transition matrix.

The D(m) matrix encodes TWO structures:
  1. MULTIPLICATIVE: Which residue classes have q-torsion, Sylow structure
  2. ADDITIVE: How consecutive prime gaps actually land on residue classes

Gemini's washing theorem applies to (1) but NOT to (2).

This script:
  - Constructs D_ideal: what D would look like if gaps respected only multiplicative symmetry
  - Computes D_ghost = D_actual - D_ideal: the pure additive defect
  - Analyzes the spectral structure of D_ghost

The Ghost is where the prime gas LIVES — the irreconcilability made manifest.

Author: Claude + Gemini + Tony (Fancyland LLC)
Date: 2026-03-25
"""

import numpy as np
from math import gcd
from typing import List, Dict, Tuple
from collections import defaultdict
import time

# =============================================================================
# PRIME UTILITIES
# =============================================================================

def sieve_primes(limit: int) -> List[int]:
    """Sieve of Eratosthenes."""
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return [i for i in range(limit + 1) if is_prime[i]]

def get_consecutive_primes(m: int, count: int = 100000) -> List[int]:
    """Get consecutive primes starting above m for gap statistics."""
    # Estimate upper bound using prime number theorem
    upper = max(m + count * 20, m * 2)
    primes = sieve_primes(upper)
    # Filter to primes > m
    return [p for p in primes if p > m][:count]

# =============================================================================
# RESIDUE CLASS ANALYSIS
# =============================================================================

def get_coprime_residues(m: int) -> List[int]:
    """Return all residue classes coprime to m."""
    return [r for r in range(1, m) if gcd(r, m) == 1]

def mod4_partition(residues: List[int]) -> Tuple[List[int], List[int]]:
    """
    Partition residues by mod-4 class.
    
    NOTE: This is NOT a multiplicative character when 4∤m!
    It's a "phantom character" that sits at the additive-multiplicative interface.
    """
    odd_block = [r for r in residues if r % 4 == 1]   # r ≡ 1 (mod 4)
    even_block = [r for r in residues if r % 4 == 3]  # r ≡ 3 (mod 4)
    return odd_block, even_block

def demonstrate_non_multiplicativity(m: int):
    """
    Prove that χ₋₁ (mod-4 partition) is not multiplicative for this m.
    """
    print(f"\n  Demonstrating non-multiplicativity of χ₋₁ at m={m}:")
    print(f"  m ≡ {m % 4} (mod 4)")
    
    if m % 4 == 0:
        print("  4 | m, so χ₋₁ IS multiplicative here.")
        return
    
    # Find counterexample: a, b coprime to m with χ(a)χ(b) ≠ χ(ab mod m)
    residues = get_coprime_residues(m)
    
    for a in residues[:20]:  # Check first 20
        for b in residues[:20]:
            ab = (a * b) % m
            if ab == 0 or gcd(ab, m) != 1:
                continue
            
            chi_a = 1 if a % 4 == 1 else -1
            chi_b = 1 if b % 4 == 1 else -1
            chi_ab = 1 if ab % 4 == 1 else -1
            
            if chi_a * chi_b != chi_ab:
                print(f"  COUNTEREXAMPLE: {a} × {b} ≡ {ab} (mod {m})")
                print(f"    χ({a}) = {chi_a}, χ({b}) = {chi_b}, χ({a})χ({b}) = {chi_a * chi_b}")
                print(f"    χ({ab}) = {chi_ab} ≠ {chi_a * chi_b}")
                print(f"  ∴ The mod-4 partition is NOT a group homomorphism!")
                return
    
    print("  No counterexample found in first 400 pairs (unexpected).")

# =============================================================================
# TRANSITION MATRIX CONSTRUCTION
# =============================================================================

def build_actual_D(m: int, primes: List[int]) -> np.ndarray:
    """
    Build the actual transition matrix D(m) from consecutive prime data.
    
    D[i,j] = count of consecutive prime pairs (p, p') where
             p ≡ residues[i] (mod m) and p' ≡ residues[j] (mod m)
    """
    residues = get_coprime_residues(m)
    res_to_idx = {r: i for i, r in enumerate(residues)}
    n = len(residues)
    
    D = np.zeros((n, n), dtype=np.float64)
    
    for i in range(len(primes) - 1):
        p, p_next = primes[i], primes[i + 1]
        r1 = p % m
        r2 = p_next % m
        
        if r1 in res_to_idx and r2 in res_to_idx:
            D[res_to_idx[r1], res_to_idx[r2]] += 1
    
    return D, residues

def build_ideal_D(m: int, primes: List[int]) -> np.ndarray:
    """
    Build the IDEAL transition matrix D_ideal(m).
    
    This is what D would look like if consecutive prime gaps respected
    ONLY the multiplicative structure — i.e., if the probability of
    transitioning from class i to class j depended only on the
    expected density, not on actual gap correlations.
    
    D_ideal[i,j] = (count of primes ≡ i) × (count of primes ≡ j) / total_transitions
    
    This represents the "null hypothesis" of multiplicative washing.
    """
    residues = get_coprime_residues(m)
    res_to_idx = {r: i for i, r in enumerate(residues)}
    n = len(residues)
    
    # Count primes in each residue class
    counts = np.zeros(n, dtype=np.float64)
    for p in primes:
        r = p % m
        if r in res_to_idx:
            counts[res_to_idx[r]] += 1
    
    total = np.sum(counts)
    
    # Under independence, D_ideal[i,j] ∝ counts[i] * counts[j]
    # Normalize to match total transition count
    D_ideal = np.outer(counts, counts)
    D_ideal *= (len(primes) - 1) / np.sum(D_ideal)
    
    return D_ideal, residues

def build_symmetric_ideal_D(m: int, primes: List[int]) -> np.ndarray:
    """
    Build a SYMMETRIC ideal matrix based on Dirichlet's theorem.
    
    Under GRH, primes equidistribute across coprime residue classes.
    The ideal symmetric matrix has D_ideal[i,j] = total_transitions / φ(m)²
    for all coprime pairs (i,j).
    """
    residues = get_coprime_residues(m)
    n = len(residues)
    
    total_transitions = len(primes) - 1
    
    # Perfect equidistribution
    D_ideal = np.full((n, n), total_transitions / (n * n), dtype=np.float64)
    
    return D_ideal, residues

# =============================================================================
# GHOST EXTRACTION AND ANALYSIS
# =============================================================================

def extract_ghost(D_actual: np.ndarray, D_ideal: np.ndarray) -> np.ndarray:
    """
    Extract the Additive Ghost: D_ghost = D_actual - D_ideal
    
    This matrix contains ONLY the gap-correlation defects.
    """
    return D_actual - D_ideal

def analyze_ghost_spectrum(D_ghost: np.ndarray, residues: List[int], 
                           block_name: str = "full"):
    """
    Analyze the spectral structure of the Ghost matrix.
    """
    print(f"\n  Ghost spectrum ({block_name}):")
    
    # Eigenvalue decomposition
    eigenvalues = np.linalg.eigvalsh(D_ghost + D_ghost.T)  # Symmetrize for real eigenvalues
    eigenvalues = np.sort(eigenvalues)[::-1]
    
    print(f"    Shape: {D_ghost.shape}")
    print(f"    Frobenius norm: {np.linalg.norm(D_ghost, 'fro'):.4f}")
    print(f"    Max |entry|: {np.max(np.abs(D_ghost)):.4f}")
    print(f"    Top 5 eigenvalues: {eigenvalues[:5]}")
    print(f"    Bottom 5 eigenvalues: {eigenvalues[-5:]}")
    
    # Rank (how many significant singular values)
    U, S, Vt = np.linalg.svd(D_ghost)
    rank_99 = np.sum(S > 0.01 * S[0])  # Count singular values > 1% of max
    print(f"    Effective rank (sv > 1% of max): {rank_99}")
    
    return eigenvalues, S

def analyze_ghost_blocks(D_ghost: np.ndarray, residues: List[int]):
    """
    Partition the Ghost by mod-4 blocks and analyze each.
    """
    odd_block, even_block = mod4_partition(residues)
    
    res_to_idx = {r: i for i, r in enumerate(residues)}
    odd_idx = [res_to_idx[r] for r in odd_block]
    even_idx = [res_to_idx[r] for r in even_block]
    
    # Extract blocks
    D_ghost_odd = D_ghost[np.ix_(odd_idx, odd_idx)]
    D_ghost_even = D_ghost[np.ix_(even_idx, even_idx)]
    D_ghost_cross = D_ghost[np.ix_(odd_idx, even_idx)]
    
    print("\n  Block analysis of Ghost matrix:")
    
    print(f"\n  ODD-ODD block ({len(odd_idx)}×{len(odd_idx)}):")
    print(f"    Frobenius norm: {np.linalg.norm(D_ghost_odd, 'fro'):.4f}")
    print(f"    Mean: {np.mean(D_ghost_odd):.6f}")
    print(f"    Trace: {np.trace(D_ghost_odd):.4f}")
    
    print(f"\n  EVEN-EVEN block ({len(even_idx)}×{len(even_idx)}):")
    print(f"    Frobenius norm: {np.linalg.norm(D_ghost_even, 'fro'):.4f}")
    print(f"    Mean: {np.mean(D_ghost_even):.6f}")
    print(f"    Trace: {np.trace(D_ghost_even):.4f}")
    
    print(f"\n  ODD-EVEN cross block ({len(odd_idx)}×{len(even_idx)}):")
    print(f"    Frobenius norm: {np.linalg.norm(D_ghost_cross, 'fro'):.4f}")
    print(f"    Mean: {np.mean(D_ghost_cross):.6f}")
    
    # The KEY diagnostic: does the Ghost have block asymmetry?
    odd_total = np.sum(D_ghost_odd)
    even_total = np.sum(D_ghost_even)
    cross_total = np.sum(D_ghost_cross)
    
    print(f"\n  Block totals (sum of entries):")
    print(f"    Odd-Odd:   {odd_total:+.4f}")
    print(f"    Even-Even: {even_total:+.4f}")
    print(f"    Cross:     {cross_total:+.4f}")
    print(f"    Asymmetry (Odd - Even): {odd_total - even_total:+.4f}")
    
    return D_ghost_odd, D_ghost_even, D_ghost_cross

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def analyze_primorial(m: int, prime_count: int = 50000):
    """
    Full Ghost analysis for a primorial m.
    """
    print("=" * 72)
    print(f"  ADDITIVE GHOST EXTRACTION: m = {m}")
    print("=" * 72)
    
    # Step 1: Demonstrate non-multiplicativity
    demonstrate_non_multiplicativity(m)
    
    # Step 2: Get prime data
    print(f"\n  Generating {prime_count} consecutive primes above {m}...")
    t0 = time.time()
    primes = get_consecutive_primes(m, prime_count)
    print(f"  Range: [{primes[0]}, {primes[-1]}]")
    print(f"  Time: {time.time() - t0:.2f}s")
    
    # Step 3: Build actual D matrix
    print(f"\n  Building actual transition matrix D({m})...")
    D_actual, residues = build_actual_D(m, primes)
    print(f"  Shape: {D_actual.shape}")
    print(f"  Total transitions: {np.sum(D_actual):.0f}")
    
    # Step 4: Build ideal D matrix (independence model)
    print(f"\n  Building ideal matrix (independence model)...")
    D_ideal, _ = build_ideal_D(m, primes)
    
    # Step 5: Build symmetric ideal (equidistribution model)
    print(f"  Building symmetric ideal (equidistribution model)...")
    D_sym_ideal, _ = build_symmetric_ideal_D(m, primes)
    
    # Step 6: Extract Ghosts
    print(f"\n  Extracting Additive Ghosts...")
    D_ghost_indep = extract_ghost(D_actual, D_ideal)
    D_ghost_equi = extract_ghost(D_actual, D_sym_ideal)
    
    print("\n" + "-" * 72)
    print("  GHOST FROM INDEPENDENCE MODEL (removes density correlations)")
    print("-" * 72)
    analyze_ghost_spectrum(D_ghost_indep, residues, "independence")
    analyze_ghost_blocks(D_ghost_indep, residues)
    
    print("\n" + "-" * 72)
    print("  GHOST FROM EQUIDISTRIBUTION MODEL (pure deviation from uniform)")
    print("-" * 72)
    analyze_ghost_spectrum(D_ghost_equi, residues, "equidistribution")
    analyze_ghost_blocks(D_ghost_equi, residues)
    
    # Step 7: Diagonal analysis (Lemke Oliver-Soundararajan)
    print("\n" + "-" * 72)
    print("  DIAGONAL SUPPRESSION (Lemke Oliver-Soundararajan)")
    print("-" * 72)
    
    diag_actual = np.diag(D_actual)
    diag_ideal = np.diag(D_ideal)
    diag_ghost = np.diag(D_ghost_indep)
    
    print(f"  Mean diagonal (actual): {np.mean(diag_actual):.4f}")
    print(f"  Mean diagonal (ideal):  {np.mean(diag_ideal):.4f}")
    print(f"  Mean diagonal (ghost):  {np.mean(diag_ghost):.4f}")
    print(f"  Diagonal suppression ratio: {np.mean(diag_actual) / np.mean(diag_ideal):.4f}")
    
    return D_actual, D_ideal, D_ghost_indep, D_ghost_equi, residues

# =============================================================================
# COMPARATIVE ANALYSIS
# =============================================================================

def compare_primorials():
    """
    Compare Ghost structure across primorials to find patterns.
    """
    primorials = [
        (30, [2, 3, 5], "balanced, no phase transition"),
        (210, [2, 3, 5, 7], "7≡1(mod 3) injects 3-torsion"),
        (2310, [2, 3, 5, 7, 11], "11≡1(mod 5) injects 5-torsion"),
    ]
    
    print("\n" + "=" * 72)
    print("  COMPARATIVE GHOST ANALYSIS")
    print("=" * 72)
    
    results = {}
    
    for m, primes_in_m, description in primorials:
        print(f"\n{'─' * 72}")
        print(f"  m = {m}: {description}")
        print(f"  Terminal prime: {primes_in_m[-1]}")
        print(f"{'─' * 72}")
        
        D_actual, D_ideal, D_ghost, _, residues = analyze_primorial(m, prime_count=30000)
        
        # Store key metrics
        odd_block, even_block = mod4_partition(residues)
        res_to_idx = {r: i for i, r in enumerate(residues)}
        odd_idx = [res_to_idx[r] for r in odd_block]
        even_idx = [res_to_idx[r] for r in even_block]
        
        D_ghost_odd = D_ghost[np.ix_(odd_idx, odd_idx)]
        D_ghost_even = D_ghost[np.ix_(even_idx, even_idx)]
        
        results[m] = {
            'ghost_asymmetry': np.sum(D_ghost_odd) - np.sum(D_ghost_even),
            'ghost_norm': np.linalg.norm(D_ghost, 'fro'),
            'diag_suppression': np.mean(np.diag(D_actual)) / np.mean(np.diag(D_ideal)),
        }
    
    print("\n" + "=" * 72)
    print("  SUMMARY TABLE")
    print("=" * 72)
    print(f"\n  {'m':>8} | {'Ghost Asymm':>12} | {'Ghost Norm':>12} | {'Diag Supp':>10}")
    print(f"  {'-'*8}-+-{'-'*12}-+-{'-'*12}-+-{'-'*10}")
    
    for m in sorted(results.keys()):
        r = results[m]
        print(f"  {m:>8} | {r['ghost_asymmetry']:>+12.4f} | {r['ghost_norm']:>12.4f} | {r['diag_suppression']:>10.4f}")

# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    # Run comparative analysis across small primorials
    compare_primorials()
    
    print("\n" + "=" * 72)
    print("  The Ghost is the scar left by the sieve.")
    print("  The defect is where addition and multiplication cannot reconcile.")
    print("=" * 72)