#!/usr/bin/env python3
"""
M=510510 ZERO-DEFECT VERIFICATION
=================================
Tests Gemini's prediction that m=510510 is perfectly balanced.

Terminal prime 17 injects no q-torsion (17 ≢ 1 mod q for q ∈ {3,5,7,11,13})
The parity cross-product washes out all historical defects.

Prediction: Δvq = 0 for all odd primes q.

Author: Claude + Gemini + Tony (Fancyland LLC)
Date: 2026-03-25
"""

import sys
from math import gcd
from typing import List, Dict, Tuple
from collections import defaultdict
import time

# =============================================================================
# CORE FUNCTIONS
# =============================================================================

def euler_phi(n: int) -> int:
    """Compute Euler's totient function φ(n)."""
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result

def v_p(n: int, p: int) -> int:
    """Compute p-adic valuation of n."""
    if n == 0:
        return float('inf')
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

def multiplicative_order(a: int, n: int) -> int:
    """Compute multiplicative order of a modulo n."""
    if gcd(a, n) != 1:
        return 0
    
    phi = euler_phi(n)
    
    # Find divisors of phi
    divisors = []
    for i in range(1, int(phi**0.5) + 1):
        if phi % i == 0:
            divisors.append(i)
            if i != phi // i:
                divisors.append(phi // i)
    divisors.sort()
    
    for d in divisors:
        if pow(a, d, n) == 1:
            return d
    return phi

def multiplicative_order_fast(a: int, n: int, phi_n: int, phi_factors: List[Tuple[int,int]]) -> int:
    """
    Fast multiplicative order using known factorization of φ(n).
    
    phi_factors is list of (prime, exponent) pairs.
    """
    if gcd(a, n) != 1:
        return 0
    
    order = phi_n
    for p, e in phi_factors:
        while order % p == 0:
            if pow(a, order // p, n) == 1:
                order //= p
            else:
                break
    return order

def factor_phi(m: int, primes_in_m: List[int]) -> List[Tuple[int, int]]:
    """
    Factor φ(m) for a primorial m = p1·p2·...·pk.
    φ(m) = (p1-1)(p2-1)...(pk-1)
    """
    from collections import Counter
    
    phi_factors = Counter()
    for p in primes_in_m:
        # Factor (p-1)
        n = p - 1
        for q in range(2, n + 1):
            if q * q > n:
                break
            while n % q == 0:
                phi_factors[q] += 1
                n //= q
        if n > 1:
            phi_factors[n] += 1
    
    return list(phi_factors.items())

# =============================================================================
# MAIN VERIFICATION
# =============================================================================

def verify_510510():
    """
    Compute order statistics for m=510510 to verify zero-defect prediction.
    """
    print("=" * 72)
    print("  M=510510 ZERO-DEFECT VERIFICATION")
    print("  Testing Gemini's Markovian Quasicrystal Hypothesis")
    print("=" * 72)
    
    # m = 2 × 3 × 5 × 7 × 11 × 13 × 17 = 510510
    primes_in_m = [2, 3, 5, 7, 11, 13, 17]
    m = 1
    for p in primes_in_m:
        m *= p
    
    phi_m = euler_phi(m)
    
    print(f"\n  m = {m}")
    print(f"  φ(m) = {phi_m}")
    print(f"  Primes in m: {primes_in_m}")
    print(f"  Terminal prime: 17")
    print()
    print("  Terminal prime analysis:")
    for q in [3, 5, 7, 11, 13]:
        print(f"    17 mod {q} = {17 % q} {'→ INJECTS' if 17 % q == 1 else '→ clean'}")
    print()
    
    # Factor φ(m) for fast order computation
    print("  Factoring φ(m) for fast order computation...")
    phi_factors = factor_phi(m, primes_in_m)
    print(f"  φ(m) = {' × '.join(f'{p}^{e}' for p,e in sorted(phi_factors))}")
    print()
    
    # Enumerate residue classes
    print(f"  Enumerating {phi_m} residue classes...")
    t0 = time.time()
    
    residues = [r for r in range(1, m) if gcd(r, m) == 1]
    
    # Partition by χ₋₁ (quadratic character)
    # For odd r: χ₋₁(r) = +1 if r ≡ 1 (mod 4), -1 if r ≡ 3 (mod 4)
    odd_block = []
    even_block = []
    
    for r in residues:
        if r % 2 == 1:  # r is odd (must be, since 2|m)
            if r % 4 == 1:
                odd_block.append(r)
            else:
                even_block.append(r)
    
    print(f"  |odd block| = {len(odd_block)}")
    print(f"  |even block| = {len(even_block)}")
    print(f"  Enumeration time: {time.time() - t0:.2f}s")
    print()
    
    # Compute order statistics
    print("  Computing multiplicative orders (this may take a few minutes)...")
    print("  Progress: ", end="", flush=True)
    
    # Track v_q(ord) distributions
    q_values = [3, 5, 7, 11, 13]
    
    odd_vq_counts = {q: defaultdict(int) for q in q_values}
    even_vq_counts = {q: defaultdict(int) for q in q_values}
    
    total = len(residues)
    checkpoint = total // 20
    
    t0 = time.time()
    for i, r in enumerate(residues):
        if checkpoint > 0 and i % checkpoint == 0:
            print(f"{100*i//total}%.", end="", flush=True)
        
        ord_r = multiplicative_order_fast(r, m, phi_m, phi_factors)
        
        is_odd = (r % 4 == 1)
        
        for q in q_values:
            vq = v_p(ord_r, q)
            if is_odd:
                odd_vq_counts[q][vq] += 1
            else:
                even_vq_counts[q][vq] += 1
    
    print(f" Done! ({time.time() - t0:.1f}s)")
    print()
    
    # Print results
    print("-" * 72)
    print("  ORDER STATISTICS RESULTS")
    print("-" * 72)
    
    defects_found = []
    
    for q in q_values:
        print(f"\n  v_{q}(ord(r,m)) distribution:")
        print(f"    v{q:<3} {'odd':>8} {'even':>8} {'diff':>8}")
        
        all_vq = set(odd_vq_counts[q].keys()) | set(even_vq_counts[q].keys())
        total_delta = 0
        
        for vq in sorted(all_vq):
            o = odd_vq_counts[q][vq]
            e = even_vq_counts[q][vq]
            d = o - e
            total_delta += vq * d  # Weighted by valuation level
            
            marker = " ***" if d != 0 else ""
            sign = "+" if d >= 0 else ""
            print(f"    {vq:<4} {o:>8} {e:>8}   {sign}{d}{marker}")
        
        print(f"    {'─'*36}")
        
        # Compute net Δvq as sum of (count_odd - count_even) weighted by vq
        # Actually, the defect is simpler: just sum of differences at each level
        simple_delta = sum(odd_vq_counts[q][v] - even_vq_counts[q][v] 
                          for v in all_vq if v > 0)
        
        print(f"    Σ(diff for v>0) = {simple_delta}")
        
        if simple_delta != 0:
            defects_found.append((q, simple_delta))
    
    # Summary
    print()
    print("=" * 72)
    print("  VERDICT")
    print("=" * 72)
    
    if not defects_found:
        print("""
  ╔══════════════════════════════════════════════════════════════════╗
  ║                                                                  ║
  ║   🎯 GEMINI'S PREDICTION CONFIRMED 🎯                            ║
  ║                                                                  ║
  ║   m = 510510 is PERFECTLY BALANCED                               ║
  ║   Δvq = 0 for all odd primes q ∈ {3, 5, 7, 11, 13}               ║
  ║                                                                  ║
  ║   The Markovian Quasicrystal hypothesis holds:                   ║
  ║   - Terminal prime 17 injects no torsion                         ║
  ║   - Parity cross-product washed out all historical defects       ║
  ║   - This is a TRANSPARENT primorial (lattice node)               ║
  ║                                                                  ║
  ╚══════════════════════════════════════════════════════════════════╝
        """)
    else:
        print(f"\n  ⚠️  DEFECTS FOUND: {defects_found}")
        print("  Gemini's prediction FALSIFIED!")
        print("  The parity cross-product does NOT fully wash out defects.")
    
    print()
    return len(defects_found) == 0

# =============================================================================
# BONUS: Compare with m=30030 (should have defect)
# =============================================================================

def verify_30030():
    """
    Quick verification of m=30030 defect for comparison.
    """
    print("=" * 72)
    print("  M=30030 DEFECT CHECK (for comparison)")
    print("=" * 72)
    
    primes_in_m = [2, 3, 5, 7, 11, 13]
    m = 30030
    phi_m = 5760
    
    print(f"\n  m = {m}, φ(m) = {phi_m}")
    print(f"  Terminal prime: 13")
    print(f"  13 mod 3 = {13 % 3} → INJECTS 3-torsion")
    print()
    
    phi_factors = factor_phi(m, primes_in_m)
    
    residues = [r for r in range(1, m) if gcd(r, m) == 1]
    
    odd_v3 = defaultdict(int)
    even_v3 = defaultdict(int)
    
    print("  Computing orders...", end="", flush=True)
    t0 = time.time()
    
    for r in residues:
        ord_r = multiplicative_order_fast(r, m, phi_m, phi_factors)
        v3 = v_p(ord_r, 3)
        
        if r % 4 == 1:
            odd_v3[v3] += 1
        else:
            even_v3[v3] += 1
    
    print(f" Done ({time.time() - t0:.1f}s)")
    print()
    
    print("  v_3(ord) distribution:")
    print(f"    v3   {'odd':>8} {'even':>8} {'diff':>8}")
    
    all_v3 = set(odd_v3.keys()) | set(even_v3.keys())
    for v3 in sorted(all_v3):
        o = odd_v3[v3]
        e = even_v3[v3]
        d = o - e
        marker = " ***" if d != 0 else ""
        sign = "+" if d >= 0 else ""
        print(f"    {v3:<4} {o:>8} {e:>8}   {sign}{d}{marker}")
    
    delta = sum(odd_v3[v] - even_v3[v] for v in all_v3 if v > 0)
    print(f"\n  Δv₃ = {delta}")
    
    if delta != 0:
        print("  ✓ Defect confirmed (as expected from 13 ≡ 1 mod 3)")
    print()

# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    # First verify 30030 has defect (sanity check)
    verify_30030()
    
    # Then test 510510
    result = verify_510510()
    
    sys.exit(0 if result else 1)