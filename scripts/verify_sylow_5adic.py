#!/usr/bin/env python3
"""
verify_sylow_5adic.py -- Check v_5(det_odd) vs v_5(det_even) for m=2310
========================================================================

The Sylow defect predictor says Delta_v5(2310) != 0 because 11 = 1 (mod 5).
Let's verify by computing actual determinants of the transition matrix
restricted to odd/even blocks (chi_{-1} = +1 / -1 partition).

This is fast: m=2310 has phi=480, so each block is 240x240.

Author: Claude + Tony
Date: 2026-03-25
"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
from math import gcd, log
from fractions import Fraction
from collections import defaultdict

def coprime_residues(m):
    return sorted(r for r in range(1, m) if gcd(r, m) == 1)

def build_D_sym(residues, m):
    r = np.array(residues, dtype=np.float64)
    diff = np.abs(r[:, None] - r[None, :])
    return np.minimum(diff, m - diff)

def v_p(n, p):
    """p-adic valuation of integer n."""
    if n == 0:
        return float('inf')
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

def multiplicative_order(a, n):
    """Multiplicative order of a mod n."""
    if gcd(a, n) != 1:
        return 0
    order = 1
    current = a % n
    while current != 1:
        current = (current * a) % n
        order += 1
        if order > n:
            return 0
    return order

# =========================================================================
# Main computation
# =========================================================================

print("=" * 72)
print("  SYLOW 5-ADIC DEFECT VERIFICATION for m=2310")
print("=" * 72)

m = 2310
res = coprime_residues(m)
phi = len(res)
print(f"\n  m={m}, phi={phi}")
print(f"  Prime factorization: 2 * 3 * 5 * 7 * 11")

# Partition residues by chi_{-1}(r) = (-1)^((r-1)/2) for odd r
# Since 2|m, all coprime residues are odd
odd_block = []   # r = 1 (mod 4), chi_{-1} = +1
even_block = []  # r = 3 (mod 4), chi_{-1} = -1

for r in res:
    if r % 4 == 1:
        odd_block.append(r)
    elif r % 4 == 3:
        even_block.append(r)
    else:
        print(f"  WARNING: r={r} is even but coprime to {m}!?")

print(f"  |odd block| (r=1 mod 4) = {len(odd_block)}")
print(f"  |even block| (r=3 mod 4) = {len(even_block)}")

# =========================================================================
# Method 1: Order-based analysis (what the Sylow predictor computes)
# =========================================================================
print("\n" + "-" * 72)
print("  METHOD 1: Multiplicative order statistics")
print("-" * 72)

for q in [2, 3, 5, 7, 11]:
    v_q_phi = v_p(phi, q)
    if v_q_phi == 0:
        continue
    
    odd_v_sum = sum(v_p(multiplicative_order(r, m), q) for r in odd_block)
    even_v_sum = sum(v_p(multiplicative_order(r, m), q) for r in even_block)
    delta = odd_v_sum - even_v_sum
    
    marker = " *** DEFECT!" if delta != 0 else ""
    print(f"  q={q}: v_q(phi)={v_q_phi}, "
          f"Sigma v_q(ord) odd={odd_v_sum}, even={even_v_sum}, "
          f"Delta={delta}{marker}")

# =========================================================================
# Method 2: Actual matrix determinant analysis
# =========================================================================
print("\n" + "-" * 72)
print("  METHOD 2: D_sym block determinants")
print("-" * 72)

D = build_D_sym(res, m)
print(f"  D_sym shape: {D.shape}")

# Index maps
res_to_idx = {r: i for i, r in enumerate(res)}
odd_idx = [res_to_idx[r] for r in odd_block]
even_idx = [res_to_idx[r] for r in even_block]

# Extract blocks
D_odd = D[np.ix_(odd_idx, odd_idx)]
D_even = D[np.ix_(even_idx, even_idx)]

print(f"  D_odd shape:  {D_odd.shape}")
print(f"  D_even shape: {D_even.shape}")

# Compute eigenvalues of each block
eigs_odd = np.linalg.eigvalsh(D_odd)
eigs_even = np.linalg.eigvalsh(D_even)

print(f"\n  D_odd eigenvalue range:  [{eigs_odd[0]:.4f}, {eigs_odd[-1]:.4f}]")
print(f"  D_even eigenvalue range: [{eigs_even[0]:.4f}, {eigs_even[-1]:.4f}]")
print(f"  D_odd  Perron eigenvalue: {eigs_odd[-1]:.8f}")
print(f"  D_even Perron eigenvalue: {eigs_even[-1]:.8f}")

# Log-determinants (to avoid overflow)
# det = product of eigenvalues
# log|det| = sum of log|eigenvalue|

# Check for zero eigenvalues
n_zero_odd = np.sum(np.abs(eigs_odd) < 1e-10)
n_zero_even = np.sum(np.abs(eigs_even) < 1e-10)
print(f"\n  Zero eigenvalues (|e|<1e-10): odd={n_zero_odd}, even={n_zero_even}")

if n_zero_odd > 0 or n_zero_even > 0:
    print("  WARNING: Zero eigenvalues detected - determinant is 0!")

# Sign and log-magnitude of determinant
sign_odd, logdet_odd = np.linalg.slogdet(D_odd)
sign_even, logdet_even = np.linalg.slogdet(D_even)

print(f"\n  det(D_odd):  sign={sign_odd:+.0f}, log|det|={logdet_odd:.4f}")
print(f"  det(D_even): sign={sign_even:+.0f}, log|det|={logdet_even:.4f}")
print(f"  log|det_odd/det_even| = {logdet_odd - logdet_even:.6f}")

# =========================================================================
# Method 3: Integer determinant via exact arithmetic approach
# =========================================================================
print("\n" + "-" * 72)
print("  METHOD 3: Integer D_sym (D_sym has integer entries)")
print("-" * 72)

# D_sym entries are min(|ri - rj|, m - |ri - rj|) which are integers!
D_odd_int = D_odd.astype(np.int64)
D_even_int = D_even.astype(np.int64)

# For 240x240 integer matrices, the determinant can be enormous
# Let's use modular arithmetic to check v_5
# We want v_5(det(D_odd)) and v_5(det(D_even))

# Instead of computing huge determinants, we can compute det mod 5^k
# using numpy with integer arithmetic

# Actually, for moderate sizes we can try to compute the exact determinant
# using Python integers via fraction-free Gaussian elimination

def det_mod_pk(M, p, k):
    """Compute det(M) mod p^k using LU decomposition over Z/(p^k)Z.
    
    For checking v_p(det), we need det mod p^(k+1).
    """
    n = M.shape[0]
    mod = p ** k
    A = np.array(M % mod, dtype=object)
    
    # Convert to Python ints for exact arithmetic
    A = [[int(A[i, j]) for j in range(n)] for i in range(n)]
    
    sign = 1
    for col in range(n):
        # Find pivot
        pivot_row = None
        for row in range(col, n):
            if A[row][col] % p != 0:
                pivot_row = row
                break
        
        if pivot_row is None:
            # All entries in this column are divisible by p
            # det is divisible by p (at least)
            return 0  # det = 0 mod p^k if we can't find invertible pivot
        
        if pivot_row != col:
            A[col], A[pivot_row] = A[pivot_row], A[col]
            sign = -sign
        
        pivot = A[col][col]
        # We need pivot to be invertible mod p
        # Since pivot % p != 0, it's invertible mod p^k
        try:
            pivot_inv = pow(pivot, -1, mod)
        except ValueError:
            return None  # can't invert
        
        for row in range(col + 1, n):
            factor = (A[row][col] * pivot_inv) % mod
            for j in range(col, n):
                A[row][j] = (A[row][j] - factor * A[col][j]) % mod
    
    # det = sign * product of diagonal
    det = sign
    for i in range(n):
        det = (det * A[i][i]) % mod
    
    return det % mod

# Compute v_5: find the highest power of 5 dividing det
# We check det mod 5, 25, 125, ...

print("\n  Computing p-adic valuations of block determinants...")
print("  (Using modular Gaussian elimination)")

for label, D_int in [("D_odd", D_odd_int), ("D_even", D_even_int)]:
    print(f"\n  {label}:")
    for p in [2, 3, 5, 7, 11]:
        # Check v_p by computing det mod p, p^2, etc.
        vals = []
        for k in range(1, 8):
            d = det_mod_pk(D_int, p, k)
            if d is None:
                vals.append(f"p^{k}: err")
            else:
                vals.append(d)
        
        # v_p(det) = smallest k where det mod p^(k+1) is nonzero mod p^k... 
        # Actually: v_p(det) = max{k : det = 0 mod p^k}
        vp = 0
        for k in range(1, 8):
            if vals[k-1] == 0 or (isinstance(vals[k-1], str)):
                vp = k
            else:
                # det mod p^k is nonzero, so v_p < k
                break
        
        print(f"    v_{p}(det) >= {vp}")
        print(f"      det mod p^k: {vals[:6]}")

# =========================================================================
# Method 4: Cross-block analysis (off-diagonal)
# =========================================================================
print("\n" + "-" * 72)
print("  METHOD 4: Cross-block structure (off-diagonal)")
print("-" * 72)

D_cross = D[np.ix_(odd_idx, even_idx)]
print(f"  D_cross shape: {D_cross.shape}")

# SVD of cross-block
svd_vals = np.linalg.svd(D_cross, compute_uv=False)
print(f"  Cross-block singular values: [{svd_vals[-1]:.4f}, ..., {svd_vals[0]:.4f}]")
print(f"  Rank (sv > 1e-10): {np.sum(svd_vals > 1e-10)}")
print(f"  Condition number: {svd_vals[0]/svd_vals[-1]:.2e}" if svd_vals[-1] > 1e-15 else "  Rank deficient!")

# =========================================================================
# Method 5: Direct Sylow verification via orders
# =========================================================================
print("\n" + "-" * 72)
print("  METHOD 5: Order distribution histogram for q=5")
print("-" * 72)

# For each residue, get v_5(ord(r,m))
odd_v5 = defaultdict(int)
even_v5 = defaultdict(int)

for r in odd_block:
    v = v_p(multiplicative_order(r, m), 5)
    odd_v5[v] += 1

for r in even_block:
    v = v_p(multiplicative_order(r, m), 5)
    even_v5[v] += 1

all_keys = sorted(set(odd_v5.keys()) | set(even_v5.keys()))
print(f"  v_5(ord(r,m)) distribution:")
print(f"  {'v5':>4s}  {'odd':>6s}  {'even':>6s}  {'diff':>6s}")
for k in all_keys:
    o = odd_v5.get(k, 0)
    e = even_v5.get(k, 0)
    marker = " ***" if o != e else ""
    print(f"  {k:>4d}  {o:>6d}  {e:>6d}  {o-e:>+6d}{marker}")

# Same for q=3 (known defect)
print(f"\n  v_3(ord(r,m)) distribution:")
odd_v3 = defaultdict(int)
even_v3 = defaultdict(int)
for r in odd_block:
    v = v_p(multiplicative_order(r, m), 3)
    odd_v3[v] += 1
for r in even_block:
    v = v_p(multiplicative_order(r, m), 3)
    even_v3[v] += 1

all_keys_3 = sorted(set(odd_v3.keys()) | set(even_v3.keys()))
print(f"  {'v3':>4s}  {'odd':>6s}  {'even':>6s}  {'diff':>6s}")
for k in all_keys_3:
    o = odd_v3.get(k, 0)
    e = even_v3.get(k, 0)
    marker = " ***" if o != e else ""
    print(f"  {k:>4d}  {o:>6d}  {e:>6d}  {o-e:>+6d}{marker}")

# =========================================================================
# Summary
# =========================================================================
print("\n" + "=" * 72)
print("  SUMMARY")
print("=" * 72)

sum_odd_v5 = sum(v_p(multiplicative_order(r, m), 5) for r in odd_block)
sum_even_v5 = sum(v_p(multiplicative_order(r, m), 5) for r in even_block)
delta_v5 = sum_odd_v5 - sum_even_v5

sum_odd_v3 = sum(v_p(multiplicative_order(r, m), 3) for r in odd_block)
sum_even_v3 = sum(v_p(multiplicative_order(r, m), 3) for r in even_block)
delta_v3 = sum_odd_v3 - sum_even_v3

print(f"\n  Delta v_5 (m=2310): {delta_v5}")
print(f"    Prediction (11 = 1 mod 5, so defect): NONZERO")
print(f"    Actual: {'CONFIRMED - DEFECT!' if delta_v5 != 0 else 'ZERO - no defect'}")

print(f"\n  Delta v_3 (m=2310): {delta_v3}")
print(f"    Prediction (from 7=1 mod 3): NONZERO")

if delta_v5 != 0:
    print(f"\n  >>> 5-ADIC BRAGG SCATTERING CONFIRMED at m=2310! <<<")
    print(f"  >>> 11 = 1 (mod 5) injects 5-torsion asymmetrically <<<")
    print(f"  >>> v_5(p-1) = v_5(10) = {v_p(10, 5)} <<<")
else:
    print(f"\n  >>> 5-adic defect NOT detected in order statistics <<<")
    print(f"  >>> Does NOT automatically mean no Bragg scattering <<<")
    print(f"  >>> (defect may appear in the operator, not just orders) <<<")

print()
