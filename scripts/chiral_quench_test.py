#!/usr/bin/env python3
"""
chiral_quench_test.py — Spectral Phase Transition in the Prime Gas
===================================================================

Constructs the Active Transport Operator:

    H(alpha) = D_sym / lambda_P  +  i*alpha * [D_sym, P_tau] / lambda_P

where D_sym is the palindromic distance matrix and P_tau is the
multiplicative inversion permutation for coprime residues mod m.

PREDICTION (Google Deep Think, March 25, 2026):
  Spectral topological phase transition at alpha = sqrt(2/3) ~ 0.8165.
    * Balanced primorials (chirality I = 0): clean block-diagonalization
    * Unbalanced primorials (I != 0): noisy, leaking spectrum

CHIRALITY TABLE:
  m=30     2*3*5           I=0  BALANCED    phi=8
  m=210    2*3*5*7         I=1  UNBALANCED  phi=48
  m=2310   2*3*5*7*11      I=2  UNBALANCED  phi=480
  m=30030  2*3*5*7*11*13   I=1  UNBALANCED  phi=5760  <-- hierarchy break

Usage:
  python chiral_quench_test.py           # full run (m up to 30030)
  python chiral_quench_test.py --quick   # m=30, 210 only (seconds)
  python chiral_quench_test.py --medium  # m up to 2310 (< 1 min)

Author: Antonio P. Matos / Fancyland LLC
Date: March 25, 2026
"""

import numpy as np
from math import gcd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import sys
import os

# ================================================================
ALPHA_CRIT = np.sqrt(2.0 / 3.0)   # 0.816496580927726...
# ================================================================


def coprime_residues(m):
    """Sorted coprime residues r with 1 <= r < m, gcd(r,m)=1."""
    return sorted(r for r in range(1, m) if gcd(r, m) == 1)


def build_D_sym(residues, m):
    """Palindromic distance matrix: D[i,j] = min(|ri-rj|, m-|ri-rj|)."""
    r = np.array(residues, dtype=np.float64)
    diff = np.abs(r[:, None] - r[None, :])
    return np.minimum(diff, m - diff)


def inversion_perm(residues, m):
    """Permutation sigma where residues[sigma[i]] = residues[i]^{-1} mod m."""
    idx = {r: i for i, r in enumerate(residues)}
    return np.array([idx[pow(r, -1, m)] for r in residues])


def commutator_fast(D, perm):
    """[D, P_tau] = D*P_tau - P_tau*D via column/row permutation (no matmul)."""
    return D[:, perm] - D[perm, :]


def perron_eigenvalue(D, tol=1e-14, max_iter=500):
    """Perron eigenvalue of non-negative symmetric D via power iteration."""
    n = D.shape[0]
    v = np.ones(n, dtype=np.float64) / np.sqrt(n)
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


def chirality(m):
    """(I, n_plus, n_minus, odd_primes) for chirality imbalance."""
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
    return abs(np_ - nm_), np_, nm_, odd


# ================================================================
# Per-modulus analysis
# ================================================================

def analyze_modulus(m, alphas):
    """Full spectral analysis of H(alpha) for modulus m."""
    I, np_, nm_, odd = chirality(m)
    res = coprime_residues(m)
    n = len(res)
    tag = "BALANCED" if I == 0 else "UNBALANCED"

    print()
    print("-" * 72)
    print("  m = {:>6d}  |  phi = {:>5d}  |  I = {}  |  {}".format(m, n, I, tag))
    print("  p = 1 (4): {} (n+={})   p = 3 (4): {} (n-={})".format(
        [p for p in odd if p % 4 == 1], np_,
        [p for p in odd if p % 4 == 3], nm_))
    print("-" * 72)

    # D_sym
    t0 = time.time()
    sys.stdout.write("  D_sym ({0}x{0})... ".format(n))
    sys.stdout.flush()
    D = build_D_sym(res, m)
    print("{:.1f}s".format(time.time() - t0))

    # Perron eigenvalue
    t0 = time.time()
    sys.stdout.write("  Perron eigenvalue... ")
    sys.stdout.flush()
    lam_P = perron_eigenvalue(D)
    print("lam_P = {:.6f}  ({:.1f}s)".format(lam_P, time.time() - t0))

    # Commutator
    perm = inversion_perm(res, m)
    C = commutator_fast(D, perm)
    C_fro = np.linalg.norm(C, 'fro')
    D_fro = np.linalg.norm(D, 'fro')
    cd_ratio = C_fro / D_fro
    skew_err = np.linalg.norm(C + C.T) / max(C_fro, 1e-15)
    print("  ||[D,P_tau]|| / ||D|| = {:.6f}  (sqrt(2/3) = {:.6f})".format(
        cd_ratio, ALPHA_CRIT))
    print("  Skew-symmetry residual: {:.2e}".format(skew_err))

    # Normalized operators
    A = D / lam_P
    B = C / lam_P
    del D, C, perm   # free memory for large m

    # --- Spectral sweep ---
    n_alpha = len(alphas)
    gaps         = np.empty(n_alpha)
    widths       = np.empty(n_alpha)
    int_gaps     = np.empty(n_alpha)
    gap_pos      = np.empty(n_alpha)
    eig_variance = np.empty(n_alpha)
    eig_store    = []

    t0 = time.time()
    print("  Sweeping {} alpha values...".format(n_alpha))

    for i, alpha in enumerate(alphas):
        H = A + (1j * alpha) * B
        eigs = np.linalg.eigvalsh(H)          # sorted ascending, real

        gaps[i]     = eigs[-1] - eigs[-2]      # Perron gap
        widths[i]   = eigs[-1] - eigs[0]       # spectral width
        eig_variance[i] = np.var(eigs)         # flat-band indicator

        consec = np.diff(eigs)
        j_max = int(np.argmax(consec))
        int_gaps[i] = consec[j_max]            # max internal gap
        gap_pos[i]  = j_max / max(len(consec) - 1, 1)

        eig_store.append(eigs)

        # progress for large matrices
        if n > 200 and ((i + 1) % 10 == 0 or i == 0):
            elapsed = time.time() - t0
            rate = elapsed / (i + 1)
            eta = rate * (n_alpha - i - 1)
            print("    [{:>3d}/{}] {:.1f}s elapsed, ~{:.0f}s ETA".format(
                i + 1, n_alpha, elapsed, eta))

    elapsed = time.time() - t0
    print("  Done: {:.1f}s ({:.2f}s/alpha)".format(elapsed, elapsed / n_alpha))

    ratios = np.where(widths > 1e-15, gaps / widths, 0.0)

    # --- IPR at critical alpha (eigenvector localization) ---
    ic = int(np.argmin(np.abs(alphas - ALPHA_CRIT)))
    H_crit = A + (1j * ALPHA_CRIT) * B
    sys.stdout.write("  Computing eigenvectors at alpha_crit for IPR... ")
    sys.stdout.flush()
    t0 = time.time()
    evals_c, evecs_c = np.linalg.eigh(H_crit)
    ipr = np.sum(np.abs(evecs_c)**4, axis=0)  # IPR per eigenvector
    mean_ipr = float(np.mean(ipr))
    max_ipr  = float(np.max(ipr))
    delocalized_ipr = 1.0 / n                  # reference: perfectly delocalized
    print("{:.1f}s".format(time.time() - t0))

    # Report
    print()
    print("  * At alpha = sqrt(2/3) ~ {:.4f}:".format(alphas[ic]))
    print("    Perron gap:       {:.6f}".format(gaps[ic]))
    print("    Max internal gap: {:.6f}  (block separation)".format(int_gaps[ic]))
    print("    Gap position:     {:.4f}  (0=bottom, 0.5=center, 1=top)".format(
        gap_pos[ic]))
    print("    Spectral width:   {:.4f}".format(widths[ic]))
    print("    Eigenvalue var:   {:.6f}  (0 = flat band)".format(eig_variance[ic]))
    print("    Mean IPR:         {:.6f}  (delocalized ref: {:.6f})".format(
        mean_ipr, delocalized_ipr))
    print("    Max IPR:          {:.6f}".format(max_ipr))
    print("    IPR ratio:        {:.2f}x delocalized".format(
        mean_ipr / delocalized_ipr if delocalized_ipr > 0 else 0))

    ipk = int(np.argmax(ratios))
    print("  * Peak gap ratio at alpha = {:.4f} (ratio = {:.6f})".format(
        alphas[ipk], ratios[ipk]))

    return dict(
        m=m, phi=n, I=I, cd_ratio=cd_ratio,
        alphas=alphas, gaps=gaps, widths=widths, ratios=ratios,
        int_gaps=int_gaps, gap_pos=gap_pos, eig_variance=eig_variance,
        eigs=eig_store, ic=ic,
        mean_ipr=mean_ipr, max_ipr=max_ipr, delocalized_ipr=delocalized_ipr,
    )


# ================================================================
# Main
# ================================================================

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Mode selection
    if '--quick' in sys.argv:
        primorials = [30, 210]
        mode = 'QUICK'
    elif '--medium' in sys.argv:
        primorials = [30, 210, 2310]
        mode = 'MEDIUM'
    else:
        primorials = [30, 210, 2310, 30030]
        mode = 'FULL'

    alphas = np.sort(np.unique(np.append(np.linspace(0.0, 1.5, 61), ALPHA_CRIT)))

    print("=" * 72)
    print("  CHIRAL QUENCH TEST  ({})".format(mode))
    print("  Double Helix Spectral Phase Transition in the Prime Gas")
    print("  Critical shear: alpha = sqrt(2/3) ~ {:.6f}".format(ALPHA_CRIT))
    print("  Moduli: {}".format(primorials))
    print("  alpha sweep: {} values in [0, 1.5]".format(len(alphas)))
    print("=" * 72)

    all_results = {}
    t_total = time.time()

    for m in primorials:
        try:
            all_results[m] = analyze_modulus(m, alphas)
        except MemoryError:
            print("  *** MemoryError at m={} -- skipping ***".format(m))
            continue

    total_time = time.time() - t_total
    print()
    print("=" * 72)
    print("  Total computation: {:.1f}s ({:.1f} min)".format(
        total_time, total_time / 60))
    print("=" * 72)

    # ============================================================
    # VERDICT
    # ============================================================
    print()
    print("  VERDICT -- Balanced vs Unbalanced at alpha = sqrt(2/3)")
    print("  " + "-" * 68)
    header = "  {:>7s} {:>3s} {:>5s} {:>9s} {:>9s} {:>7s} {:>8s} {:>7s} {:>7s}"
    print(header.format(
        "m", "I", "phi", "P.gap", "int.gap", "gapPos", "ratio", "C/D", "IPR_r"))
    print("  " + "-" * 68)
    for m in primorials:
        if m not in all_results:
            continue
        r = all_results[m]
        ic = r['ic']
        bal = "*" if r['I'] == 0 else " "
        ipr_r = r['mean_ipr'] / r['delocalized_ipr'] if r['delocalized_ipr'] > 0 else 0
        print("  {}{:>6d} {:>3d} {:>5d} {:>9.5f} {:>9.5f} {:>7.3f} {:>8.5f} {:>7.4f} {:>7.1f}x".format(
            bal, m, r['I'], r['phi'],
            r['gaps'][ic], r['int_gaps'][ic], r['gap_pos'][ic],
            r['ratios'][ic], r['cd_ratio'], ipr_r))
    print("  * = chirally balanced (I=0)")
    print("  IPR_r = mean IPR / (1/phi), higher = more localized")

    # ============================================================
    # SAVE TEXT SUMMARY
    # ============================================================
    summary_path = os.path.join(script_dir, 'chiral_quench_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("CHIRAL QUENCH TEST RESULTS\n")
        f.write("Date: {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S")))
        f.write("Critical alpha = sqrt(2/3) = {:.10f}\n\n".format(ALPHA_CRIT))
        for m in primorials:
            if m not in all_results:
                continue
            r = all_results[m]
            ic = r['ic']
            ipr_r = r['mean_ipr'] / r['delocalized_ipr'] if r['delocalized_ipr'] > 0 else 0
            tag = "BALANCED" if r['I'] == 0 else "UNBALANCED"
            f.write("m={}, phi={}, I={}, {}\n".format(m, r['phi'], r['I'], tag))
            f.write("  ||[D,P_tau]||/||D|| = {:.8f}\n".format(r['cd_ratio']))
            f.write("  At alpha_crit:\n")
            f.write("    Perron gap     = {:.8f}\n".format(r['gaps'][ic]))
            f.write("    Internal gap   = {:.8f}\n".format(r['int_gaps'][ic]))
            f.write("    Gap position   = {:.6f}\n".format(r['gap_pos'][ic]))
            f.write("    Spectral width = {:.8f}\n".format(r['widths'][ic]))
            f.write("    Eigenvalue var = {:.8f}\n".format(r['eig_variance'][ic]))
            f.write("    Mean IPR       = {:.8f} ({:.1f}x delocalized)\n".format(
                r['mean_ipr'], ipr_r))
            f.write("\n")
    print("\n  Summary: {}".format(summary_path))

    # ============================================================
    # PLOTS
    # ============================================================
    print("  Generating plots...")
    colors = {30: '#27ae60', 210: '#e74c3c', 2310: '#2980b9', 30030: '#8e44ad'}

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(
        'Chiral Quench Test\n'
        r'Predicted spectral quench at $\alpha = \sqrt{{2/3}}$'
        ' = {:.4f}'.format(ALPHA_CRIT),
        fontsize=14, fontweight='bold', y=0.98)

    avail = [m for m in primorials if m in all_results]

    # Panel A: Max internal gap (block separation detector)
    ax = axes[0, 0]
    for m in avail:
        r = all_results[m]
        sym = "* " if r['I'] == 0 else "  "
        ax.plot(r['alphas'], r['int_gaps'], color=colors[m], lw=1.8,
                label='{}m={} (I={}, phi={})'.format(sym, m, r['I'], r['phi']))
    ax.axvline(ALPHA_CRIT, color='red', ls='--', alpha=0.6, lw=1.2,
               label=r'$\alpha=\sqrt{2/3}$')
    ax.set_xlabel(r'$\alpha$ (shear amplitude)')
    ax.set_ylabel('Max Consecutive Eigenvalue Gap')
    ax.set_title('A.  Block Separation Indicator')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.25)

    # Panel B: Gap ratio
    ax = axes[0, 1]
    for m in avail:
        r = all_results[m]
        ax.plot(r['alphas'], r['ratios'], color=colors[m], lw=1.8,
                label='m={} (I={})'.format(m, r['I']))
    ax.axvline(ALPHA_CRIT, color='red', ls='--', alpha=0.6, lw=1.2)
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel('Perron Gap / Spectral Width')
    ax.set_title('B.  Normalized Perron Gap Ratio')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.25)

    # Panel C: Eigenvalue density at alpha_crit
    ax = axes[1, 0]
    for m in avail:
        r = all_results[m]
        eigs = r['eigs'][r['ic']]
        span = eigs[-1] - eigs[0]
        e_norm = (eigs - eigs[0]) / span if span > 1e-15 else eigs - eigs[0]
        nbins = min(80, max(8, r['phi'] // 5))
        ax.hist(e_norm, bins=nbins, alpha=0.45, color=colors[m], density=True,
                label='m={} (I={})'.format(m, r['I']))
    ax.set_xlabel('Normalized eigenvalue')
    ax.set_ylabel('Density')
    ax.set_title(r'C.  Eigenvalue Density at $\alpha = \sqrt{2/3}$')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.25)

    # Panel D: Eigenvalue variance (flat-band detector)
    ax = axes[1, 1]
    for m in avail:
        r = all_results[m]
        # normalize variance by its value at alpha=0
        var0 = r['eig_variance'][0] if r['eig_variance'][0] > 1e-15 else 1.0
        ax.plot(r['alphas'], r['eig_variance'] / var0, color=colors[m], lw=1.8,
                label='m={} (I={})'.format(m, r['I']))
    ax.axvline(ALPHA_CRIT, color='red', ls='--', alpha=0.6, lw=1.2)
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'Var($\lambda$) / Var($\lambda$)$_{\alpha=0}$')
    ax.set_title('D.  Eigenvalue Variance (flat-band detector)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.25)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    p1 = os.path.join(script_dir, 'chiral_quench_results.png')
    plt.savefig(p1, dpi=150, bbox_inches='tight')
    print("  Saved: {}".format(p1))

    # Bonus: eigenvalue scatter at alpha_crit
    ncols = len(avail)
    fig2, axes2 = plt.subplots(1, ncols, figsize=(4.5 * ncols, 4.5), squeeze=False)
    fig2.suptitle(
        r'Eigenvalue Spectrum at $\alpha = \sqrt{2/3}$ -- Block Structure',
        fontsize=13, fontweight='bold')
    for j, m in enumerate(avail):
        ax = axes2[0, j]
        r = all_results[m]
        eigs = r['eigs'][r['ic']]
        ax.plot(range(len(eigs)), eigs, '.', color=colors[m],
                ms=max(1.0, 6 - j * 1.2))
        tag = "BAL" if r['I'] == 0 else "I={}".format(r['I'])
        ax.set_title('m={} ({}, phi={})'.format(m, tag, r['phi']), fontsize=10)
        ax.set_xlabel('Index')
        ax.set_ylabel('Eigenvalue')
        ax.grid(True, alpha=0.25)
    plt.tight_layout()
    p2 = os.path.join(script_dir, 'chiral_quench_spectrum.png')
    plt.savefig(p2, dpi=150, bbox_inches='tight')
    print("  Saved: {}".format(p2))

    # Bonus 2: Gap position heatmap
    fig3, ax3 = plt.subplots(figsize=(12, 5))
    for m in avail:
        r = all_results[m]
        sym = "* " if r['I'] == 0 else "  "
        ax3.plot(r['alphas'], r['gap_pos'], color=colors[m], lw=1.8,
                 label='{}m={} (I={})'.format(sym, m, r['I']))
    ax3.axvline(ALPHA_CRIT, color='red', ls='--', alpha=0.6, lw=1.2,
                label=r'$\alpha=\sqrt{2/3}$')
    ax3.axhline(0.5, color='gray', ls=':', alpha=0.4, label='center split')
    ax3.set_xlabel(r'$\alpha$')
    ax3.set_ylabel('Position of max gap (0=bottom, 0.5=center, 1=top)')
    ax3.set_title('Gap Position: Where the Block Separation Occurs')
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.25)
    ax3.set_ylim(-0.05, 1.05)
    plt.tight_layout()
    p3 = os.path.join(script_dir, 'chiral_quench_gap_position.png')
    plt.savefig(p3, dpi=150, bbox_inches='tight')
    print("  Saved: {}".format(p3))

    print()
    print("  === ANALYSIS COMPLETE ===")
    print("  Look for:")
    print("    1. Sharp feature at alpha ~ {:.4f} in balanced (*) curves".format(
        ALPHA_CRIT))
    print("    2. Block separation (panel A) stronger for balanced moduli")
    print("    3. Clean bimodal density (panel C) for balanced vs noisy otherwise")
    print("    4. Variance drop (panel D) at critical alpha for balanced")
    print("    5. Gap position near 0.5 (center split) for balanced")


if __name__ == '__main__':
    main()
