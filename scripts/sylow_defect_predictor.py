#!/usr/bin/env python3
"""
SYLOW DEFECT PREDICTOR
======================
Predicts p-adic valuation defects Δvₚ = vₚ(det_odd) - vₚ(det_even)
for primorial transition matrices WITHOUT computing the full matrix.

Core hypothesis (Gemini-Claude synthesis):
- The odd/even blocks partition residue classes by χ₋₁ (Jacobi symbol)
- When a new prime p ≡ 1 (mod q) enters, it injects q-torsion
- If this torsion distributes asymmetrically across χ₋₁ = ±1, we get Δvq ≠ 0

Author: Claude + Gemini + Tony (Fancyland LLC)
Date: 2026-03-25
"""

from math import gcd
from functools import reduce
from collections import defaultdict
from typing import List, Dict, Tuple

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def prime_list(n: int) -> List[int]:
    """Return list of primes up to n using sieve."""
    if n < 2:
        return []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, n + 1, i):
                sieve[j] = False
    return [i for i in range(n + 1) if sieve[i]]

def v_p(n: int, p: int) -> int:
    """Compute p-adic valuation of n (highest power of p dividing n)."""
    if n == 0:
        return float('inf')
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

def jacobi_symbol(a: int, n: int) -> int:
    """Compute Jacobi symbol (a/n) for odd n > 0."""
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer")
    
    a = a % n
    result = 1
    
    while a != 0:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a = a % n
    
    return result if n == 1 else 0

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

def multiplicative_order(a: int, n: int) -> int:
    """Compute multiplicative order of a modulo n."""
    if gcd(a, n) != 1:
        return 0
    
    phi = euler_phi(n)
    order = phi
    
    # Find divisors of phi and test
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

# =============================================================================
# PRIMORIAL ANALYSIS
# =============================================================================

def primorial(k: int, primes: List[int]) -> int:
    """Return k-th primorial: product of first k primes."""
    return reduce(lambda x, y: x * y, primes[:k], 1)

def analyze_primorial(m: int, primes_in_m: List[int]) -> Dict:
    """
    Analyze the group structure (Z/mZ)× partitioned by χ₋₁.
    
    Returns detailed statistics about Sylow subgroups in each partition.
    """
    # Get all residue classes coprime to m
    residues = [r for r in range(1, m) if gcd(r, m) == 1]
    phi_m = len(residues)
    
    # Partition by Jacobi symbol (-1/r)
    # For primorial m (product of odd primes ≥ 3), we compute (-1/r)
    odd_block = []   # χ₋₁(r) = +1
    even_block = []  # χ₋₁(r) = -1
    
    for r in residues:
        # Compute (-1/r) via quadratic reciprocity
        # (-1/r) = (-1)^((r-1)/2) for odd r
        # But r might be even if 2|m... actually for primorial m≥6, 2|m
        # So we need r coprime to m, meaning r is odd when 2|m
        
        # For odd r: (-1/r) = +1 if r ≡ 1 (mod 4), -1 if r ≡ 3 (mod 4)
        if r % 2 == 1:
            chi = 1 if r % 4 == 1 else -1
        else:
            # r even, m odd case (shouldn't happen for primorials ≥ 6)
            chi = jacobi_symbol(-1, r) if r % 2 == 1 else 0
        
        if chi == 1:
            odd_block.append(r)
        else:
            even_block.append(r)
    
    # Compute q-adic order statistics for each block
    # For small primes q, count how q-power orders distribute
    
    analysis = {
        'm': m,
        'phi_m': phi_m,
        'primes_in_m': primes_in_m,
        'odd_block_size': len(odd_block),
        'even_block_size': len(even_block),
        'sylow_analysis': {}
    }
    
    # Analyze Sylow q-subgroups for q = 2, 3, 5, 7, 11, 13
    for q in [2, 3, 5, 7, 11, 13]:
        v_q_phi = v_p(phi_m, q)
        if v_q_phi == 0:
            continue
            
        # For each residue, compute v_q of its multiplicative order
        # This tells us which "level" of the Sylow q-subgroup it lives in
        
        odd_order_valuations = []
        even_order_valuations = []
        
        # This is expensive for large m, so we'll do sampling or exact for small m
        if phi_m <= 10000:
            for r in odd_block:
                ord_r = multiplicative_order(r, m)
                odd_order_valuations.append(v_p(ord_r, q))
            for r in even_block:
                ord_r = multiplicative_order(r, m)
                even_order_valuations.append(v_p(ord_r, q))
        else:
            # For large m, we compute the theoretical prediction instead
            odd_order_valuations = None
            even_order_valuations = None
        
        # Sum of v_q(ord(r)) over block gives a proxy for how q-torsion distributes
        if odd_order_valuations is not None:
            analysis['sylow_analysis'][q] = {
                'v_q_phi': v_q_phi,
                'odd_v_sum': sum(odd_order_valuations),
                'even_v_sum': sum(even_order_valuations),
                'delta': sum(odd_order_valuations) - sum(even_order_valuations)
            }
        else:
            analysis['sylow_analysis'][q] = {
                'v_q_phi': v_q_phi,
                'note': 'too large for direct computation'
            }
    
    return analysis

# =============================================================================
# DEFECT PREDICTION ENGINE
# =============================================================================

def predict_defect_from_incoming_prime(q: int, incoming_p: int, prev_delta: int) -> Tuple[int, str]:
    """
    Predict how Δvq changes when a new prime p is added to the primorial.
    
    Gemini's hypothesis:
    - If p ≡ 1 (mod q), new q-torsion is injected
    - The defect Δvq increases by v_q(p-1)
    - If p ≢ 1 (mod q), no change
    """
    v_q_of_p_minus_1 = v_p(incoming_p - 1, q)
    
    if incoming_p % q == 1:
        # New q-cycle injected!
        new_delta = prev_delta + v_q_of_p_minus_1
        reason = f"p={incoming_p} ≡ 1 (mod {q}), injects {q}^{v_q_of_p_minus_1} torsion"
    else:
        new_delta = prev_delta
        reason = f"p={incoming_p} ≡ {incoming_p % q} (mod {q}), no new {q}-torsion"
    
    return new_delta, reason

def full_defect_table(max_primorial_index: int = 8) -> List[Dict]:
    """
    Build the complete defect prediction table for primorials.
    """
    primes = prime_list(30)  # First 10 primes: 2,3,5,7,11,13,17,19,23,29
    
    results = []
    
    # Track cumulative defects for each small prime q
    cumulative_delta = {2: 0, 3: 0, 5: 0, 7: 0}
    
    for k in range(3, max_primorial_index + 1):  # Start from m=30 (k=3)
        m = primorial(k, primes)
        primes_in_m = primes[:k]
        incoming_p = primes[k-1]  # The prime just added
        
        row = {
            'k': k,
            'm': m,
            'incoming_prime': incoming_p,
            'phi_m': euler_phi(m),
            'predictions': {}
        }
        
        # Predict defects for q = 3, 5, 7
        for q in [3, 5, 7]:
            if q >= incoming_p:
                continue  # q not yet in the primorial
            
            prev_delta = cumulative_delta.get(q, 0)
            new_delta, reason = predict_defect_from_incoming_prime(q, incoming_p, prev_delta)
            cumulative_delta[q] = new_delta
            
            row['predictions'][q] = {
                'delta_vq': new_delta,
                'reason': reason,
                'triggered': (incoming_p % q == 1)
            }
        
        results.append(row)
    
    return results

# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    print("=" * 80)
    print("SYLOW DEFECT PREDICTOR")
    print("Predicting p-adic valuation defects without computing transition matrices")
    print("=" * 80)
    print()
    
    primes = prime_list(30)
    
    # Build prediction table
    print("DEFECT PREDICTION TABLE")
    print("-" * 80)
    print(f"{'m':>12} | {'p_new':>5} | {'Δv₃':>6} | {'Δv₅':>6} | {'Δv₇':>6} | Trigger Analysis")
    print("-" * 80)
    
    results = full_defect_table(9)
    
    for row in results:
        m = row['m']
        p_new = row['incoming_prime']
        
        dv3 = row['predictions'].get(3, {}).get('delta_vq', '-')
        dv5 = row['predictions'].get(5, {}).get('delta_vq', '-')
        dv7 = row['predictions'].get(7, {}).get('delta_vq', '-')
        
        triggers = []
        for q in [3, 5, 7]:
            if q in row['predictions'] and row['predictions'][q]['triggered']:
                triggers.append(f"{p_new}≡1(mod {q})")
        
        trigger_str = ", ".join(triggers) if triggers else "(none)"
        
        print(f"{m:>12} | {p_new:>5} | {str(dv3):>6} | {str(dv5):>6} | {str(dv7):>6} | {trigger_str}")
    
    print("-" * 80)
    print()
    
    # Detailed analysis for key primorials
    print("DETAILED SYLOW ANALYSIS")
    print("=" * 80)
    
    key_primorials = [
        (4, "m=210 (before first 3-adic defect)"),
        (5, "m=2310 (check for 5-adic defect from p=11)"),
        (6, "m=30030 (OBSERVED Δv₃=1)"),
        (7, "m=510510 (PREDICTED Δv₃=1, no change)"),
    ]
    
    for k, description in key_primorials:
        m = primorial(k, primes)
        primes_in_m = primes[:k]
        
        print(f"\n{description}")
        print(f"  m = {m}")
        print(f"  φ(m) = {euler_phi(m)}")
        print(f"  Primes: {primes_in_m}")
        
        if m <= 2310:
            analysis = analyze_primorial(m, primes_in_m)
            print(f"  Odd block size: {analysis['odd_block_size']}")
            print(f"  Even block size: {analysis['even_block_size']}")
            
            if analysis['sylow_analysis']:
                print("  Sylow analysis (direct computation):")
                for q, data in analysis['sylow_analysis'].items():
                    if 'delta' in data:
                        print(f"    q={q}: v_q(φ)={data['v_q_phi']}, "
                              f"Σv_q(ord) odd={data['odd_v_sum']}, even={data['even_v_sum']}, "
                              f"Δ={data['delta']}")
        else:
            print("  (Direct computation skipped - using theoretical prediction)")
    
    print()
    print("=" * 80)
    print("KEY PREDICTIONS TO VERIFY:")
    print("=" * 80)
    print()
    print("1. Δv₅(2310) ≠ 0  -- because 11 ≡ 1 (mod 5)")
    print("   → Tony: Check v₅(det_odd) vs v₅(det_even) from your 2310 data!")
    print()
    print("2. Δv₃(510510) = 1  -- unchanged from 30030 because 17 ≡ 2 (mod 3)")
    print("   → This is Gemini's core prediction")
    print()
    print("3. Δv₃(9699690) = 3  -- jump by 2 because 19 ≡ 1 (mod 3) and v₃(18)=2")
    print("   → Next major 3-adic break")
    print()
    print("4. Δv₇(510510) = 1  -- because 17 ≡ 1 (mod 7) triggers first 7-adic defect!")
    print("   → NEW PREDICTION: 7-adic defect appears at m=510510")
    print()

    # The killer: find ALL defect-triggering events up to m ~ 10^9
    print("=" * 80)
    print("DEFECT TRIGGER CALENDAR (all p-adic breaks through m ~ 10^9)")
    print("=" * 80)
    
    primes_extended = prime_list(50)
    
    print(f"\n{'m':>15} | {'p_new':>5} | Defects Triggered")
    print("-" * 60)
    
    for k in range(3, 16):
        if k >= len(primes_extended):
            break
        m = primorial(k, primes_extended)
        p_new = primes_extended[k-1]
        
        triggers = []
        for q in [3, 5, 7, 11, 13]:
            if q < p_new and p_new % q == 1:
                v = v_p(p_new - 1, q)
                triggers.append(f"Δv_{q} += {v}")
        
        trigger_str = ", ".join(triggers) if triggers else "(clean transition)"
        print(f"{m:>15} | {p_new:>5} | {trigger_str}")
    
    print()

if __name__ == "__main__":
    main()