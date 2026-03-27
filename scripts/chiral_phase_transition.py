#!/usr/bin/env python3
"""
chiral_phase_transition.py — Fine-Resolution Phase Transition Hunt
===================================================================

High-resolution alpha sweep around the gap position crash at alpha ~ 1.2.

FALSIFIABLE PREDICTION (Claude Opus 4.6, March 25, 2026):
  The gap position inflection point occurs at alpha_c = sqrt(3/2) = 1.224744...
  across all primorials >= 210.

  If true: the irreconcilability coupling constant sqrt(2/3) and the
  critical shear sqrt(3/2) are multiplicative inverses.
  kappa * alpha_c = 1.

STRATEGY:
  1. Fine sweep: 200 points in [1.0, 1.45] (include sqrt(3/2) as mandatory)
  2. Context: 20 points in [0.5, 1.0] to see the approach
  3. At each alpha: eigenvalues only (fast)
  4. At selected alpha values near transition: full eigenvectors for IPR
  5. Numerical derivative of gap position to find the crash point

Usage:
  python chiral_phase_transition.py           # full (m up to 30030)
  python chiral_phase_transition.py --quick   # m=30, 210 only
  python chiral_phase_transition.py --medium  # m up to 2310

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
ALPHA_CRIT = np.sqrt(3.0 / 2.0)   # 1.224744871391589...
KAPPA      = np.sqrt(2.0 / 3.0)   # 0.816496580927726...
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
    """[D, P_tau] = D*P_tau - P_tau*D via column/row permutation."""
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


def gap_position(eigs):
    """Position of largest consecutive gap: 0=bottom, 0.5=center, 1=top."""
    consec = np.diff(eigs)
    j_max = int(np.argmax(consec))
    return j_max / max(len(consec) - 1, 1)


def find_crash_alpha(alphas, positions):
    """Find the alpha where gap position drops through 0.5 (midpoint crash).
    Uses linear interpolation between the two bracketing alpha values."""
    for i in range(len(positions) - 1):
        if positions[i] > 0.5 and positions[i + 1] <= 0.5:
            # Linear interpolation
            frac = (0.5 - positions[i]) / (positions[i + 1] - positions[i])
            return alphas[i] + frac * (alphas[i + 1] - alphas[i])
    return None


def find_steepest_drop(alphas, positions):
    """Find alpha of steepest descent in gap position (most negative derivative)."""
    derivs = np.diff(positions) / np.diff(alphas)
    j = int(np.argmin(derivs))  # most negative
    return 0.5 * (alphas[j] + alphas[j + 1]), derivs[j]


# ================================================================
# Per-modulus fine sweep
# ================================================================

def analyze_fine(m, alphas, ipr_alphas=None):
    """Fine-resolution spectral sweep for the phase transition."""
    I, np_, nm_, odd = chirality(m)
    res = coprime_residues(m)
    n = len(res)
    tag = "BALANCED" if I == 0 else "UNBALANCED"

    print()
    print("=" * 72)
    print("  m = {:>6d}  |  phi = {:>5d}  |  I = {}  |  {}".format(m, n, I, tag))
    print("=" * 72)

    # Build matrices
    t0 = time.time()
    sys.stdout.write("  Building D_sym ({0}x{0})... ".format(n))
    sys.stdout.flush()
    D = build_D_sym(res, m)
    print("{:.1f}s".format(time.time() - t0))

    sys.stdout.write("  Perron eigenvalue... ")
    sys.stdout.flush()
    t0 = time.time()
    lam_P = perron_eigenvalue(D)
    print("lam_P = {:.6f}  ({:.1f}s)".format(lam_P, time.time() - t0))

    perm = inversion_perm(res, m)
    C = commutator_fast(D, perm)

    # Normalized operators
    A = D / lam_P
    B = C / lam_P
    del D, C, perm

    # --- Main sweep: eigenvalues only ---
    n_alpha = len(alphas)
    gap_pos_arr   = np.empty(n_alpha)
    perron_gap    = np.empty(n_alpha)
    max_int_gap   = np.empty(n_alpha)
    eig_var       = np.empty(n_alpha)

    t0 = time.time()
    print("  Sweeping {} alpha values in [{:.3f}, {:.3f}]...".format(
        n_alpha, alphas[0], alphas[-1]))

    for i, alpha in enumerate(alphas):
        H = A + (1j * alpha) * B
        eigs = np.linalg.eigvalsh(H)

        gap_pos_arr[i] = gap_position(eigs)
        perron_gap[i]  = eigs[-1] - eigs[-2]
        consec = np.diff(eigs)
        max_int_gap[i] = float(np.max(consec))
        eig_var[i]     = np.var(eigs)

        # Progress for large matrices
        if n > 200 and ((i + 1) % 20 == 0 or i == 0):
            elapsed = time.time() - t0
            rate = elapsed / (i + 1)
            eta = rate * (n_alpha - i - 1)
            print("    [{:>3d}/{}] {:.1f}s elapsed, ~{:.0f}s ETA".format(
                i + 1, n_alpha, elapsed, eta))

    elapsed = time.time() - t0
    print("  Done: {:.1f}s ({:.2f}s/alpha)".format(elapsed, elapsed / max(n_alpha, 1)))

    # --- Find the crash point ---
    crash_alpha = find_crash_alpha(alphas, gap_pos_arr)
    steepest_alpha, steepest_deriv = find_steepest_drop(alphas, gap_pos_arr)

    print()
    print("  TRANSITION ANALYSIS:")
    print("  Predicted alpha_c = sqrt(3/2) = {:.10f}".format(ALPHA_CRIT))
    if crash_alpha is not None:
        error_pct = 100 * (crash_alpha - ALPHA_CRIT) / ALPHA_CRIT
        print("  Measured crash (gap_pos = 0.5): alpha = {:.10f}".format(crash_alpha))
        print("  Error from sqrt(3/2): {:.6f}  ({:+.4f}%)".format(
            crash_alpha - ALPHA_CRIT, error_pct))
    else:
        print("  Crash (gap_pos = 0.5): NOT REACHED in sweep range")
    print("  Steepest descent at alpha = {:.10f}  (deriv = {:.2f})".format(
        steepest_alpha, steepest_deriv))
    error_steep = steepest_alpha - ALPHA_CRIT
    print("  Steepest vs sqrt(3/2): {:+.10f}  ({:+.4f}%)".format(
        error_steep, 100 * error_steep / ALPHA_CRIT))

    # --- IPR at selected alpha values ---
    ipr_results = {}
    if ipr_alphas is not None and len(ipr_alphas) > 0:
        print()
        print("  Computing IPR at {} selected alpha values...".format(len(ipr_alphas)))
        for ja, alpha_ipr in enumerate(ipr_alphas):
            sys.stdout.write("    alpha = {:.6f} ... ".format(alpha_ipr))
            sys.stdout.flush()
            t0 = time.time()
            H = A + (1j * alpha_ipr) * B
            evals, evecs = np.linalg.eigh(H)
            ipr = np.sum(np.abs(evecs)**4, axis=0)
            mean_ipr = float(np.mean(ipr))
            max_ipr  = float(np.max(ipr))
            delocalized = 1.0 / n
            ipr_results[alpha_ipr] = {
                'mean': mean_ipr, 'max': max_ipr,
                'ratio': mean_ipr / delocalized if delocalized > 0 else 0,
            }
            print("{:.1f}s  (IPR ratio: {:.2f}x)".format(
                time.time() - t0, mean_ipr / delocalized))

    # Value at the exact predicted critical point
    ic = int(np.argmin(np.abs(alphas - ALPHA_CRIT)))
    print()
    print("  At alpha = sqrt(3/2) ~ {:.6f}:".format(alphas[ic]))
    print("    Gap position:  {:.6f}".format(gap_pos_arr[ic]))
    print("    Perron gap:    {:.6f}".format(perron_gap[ic]))
    print("    Max int gap:   {:.6f}".format(max_int_gap[ic]))
    print("    Eigenvalue var: {:.8f}".format(eig_var[ic]))

    return dict(
        m=m, phi=n, I=I,
        alphas=alphas, gap_pos=gap_pos_arr,
        perron_gap=perron_gap, max_int_gap=max_int_gap,
        eig_var=eig_var,
        crash_alpha=crash_alpha,
        steepest_alpha=steepest_alpha, steepest_deriv=steepest_deriv,
        ipr_results=ipr_results,
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

    # Build alpha grid: fine in [1.0, 1.45], coarser context in [0.5, 1.0]
    context = np.linspace(0.5, 0.99, 20)       # context approach
    fine    = np.linspace(1.0, 1.45, 200)       # high-resolution transition zone
    mandatory = np.array([ALPHA_CRIT, KAPPA])   # exact predicted values
    alphas = np.sort(np.unique(np.concatenate([context, fine, mandatory])))

    # IPR probe values: below, at, and above the predicted transition
    ipr_probes = [1.0, 1.15, ALPHA_CRIT, 1.30, 1.40]

    print("=" * 72)
    print("  CHIRAL PHASE TRANSITION HUNT  ({})".format(mode))
    print("  Fine-resolution sweep around gap position crash")
    print("  Predicted: alpha_c = sqrt(3/2) = {:.10f}".format(ALPHA_CRIT))
    print("  Reciprocal: kappa  = sqrt(2/3) = {:.10f}".format(KAPPA))
    print("  Product: alpha_c * kappa = {:.15f}  (should be 1.0)".format(
        ALPHA_CRIT * KAPPA))
    print("  Moduli: {}".format(primorials))
    print("  Alpha sweep: {} values in [{:.2f}, {:.2f}]".format(
        len(alphas), alphas[0], alphas[-1]))
    print("  Fine zone: 200 values in [1.00, 1.45]")
    print("  IPR probes: {}".format([round(a, 6) for a in ipr_probes]))
    print("=" * 72)

    all_results = {}
    t_total = time.time()

    for m in primorials:
        # Only compute IPR for smaller moduli (eigenvector solve is expensive)
        if m <= 2310:
            ipr = ipr_probes
        else:
            # For m=30030, only IPR at the predicted critical point
            ipr = [ALPHA_CRIT]
        try:
            all_results[m] = analyze_fine(m, alphas, ipr_alphas=ipr)
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
    avail = [m for m in primorials if m in all_results]

    print()
    print("  VERDICT — Transition Location vs sqrt(3/2)")
    print("  " + "-" * 68)
    header = "  {:>7s} {:>3s} {:>5s} {:>14s} {:>14s} {:>10s} {:>10s}"
    print(header.format(
        "m", "I", "phi", "crash_alpha", "steepest", "err_crash%", "err_steep%"))
    print("  " + "-" * 68)

    for m in avail:
        r = all_results[m]
        crash_str = "{:.10f}".format(r['crash_alpha']) if r['crash_alpha'] else "   N/A"
        steep_str = "{:.10f}".format(r['steepest_alpha'])
        err_c = "N/A"
        if r['crash_alpha'] is not None:
            err_c = "{:+.4f}%".format(
                100 * (r['crash_alpha'] - ALPHA_CRIT) / ALPHA_CRIT)
        err_s = "{:+.4f}%".format(
            100 * (r['steepest_alpha'] - ALPHA_CRIT) / ALPHA_CRIT)
        print("  {:>7d} {:>3d} {:>5d} {:>14s} {:>14s} {:>10s} {:>10s}".format(
            m, r['I'], r['phi'], crash_str, steep_str, err_c, err_s))

    print("  " + "-" * 68)
    print("  Prediction: alpha_c = sqrt(3/2) = {:.10f}".format(ALPHA_CRIT))

    # ============================================================
    # SAVE SUMMARY
    # ============================================================
    summary_path = os.path.join(script_dir, 'phase_transition_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("CHIRAL PHASE TRANSITION HUNT\n")
        f.write("Date: {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S")))
        f.write("Predicted alpha_c = sqrt(3/2) = {:.15f}\n".format(ALPHA_CRIT))
        f.write("Reciprocal kappa  = sqrt(2/3) = {:.15f}\n".format(KAPPA))
        f.write("\n")
        for m in avail:
            r = all_results[m]
            tag = "BALANCED" if r['I'] == 0 else "UNBALANCED"
            f.write("m={}, phi={}, I={}, {}\n".format(m, r['phi'], r['I'], tag))
            if r['crash_alpha']:
                f.write("  Crash (gap=0.5): alpha = {:.15f}\n".format(r['crash_alpha']))
                f.write("  Error from sqrt(3/2): {:+.2e}\n".format(
                    r['crash_alpha'] - ALPHA_CRIT))
            else:
                f.write("  Crash: NOT REACHED\n")
            f.write("  Steepest descent: alpha = {:.15f}  (deriv = {:.4f})\n".format(
                r['steepest_alpha'], r['steepest_deriv']))
            f.write("  Steepest error from sqrt(3/2): {:+.2e}\n".format(
                r['steepest_alpha'] - ALPHA_CRIT))
            if r['ipr_results']:
                f.write("  IPR at selected alphas:\n")
                for alpha_ipr in sorted(r['ipr_results'].keys()):
                    ir = r['ipr_results'][alpha_ipr]
                    f.write("    alpha={:.6f}: mean={:.8f}, ratio={:.2f}x\n".format(
                        alpha_ipr, ir['mean'], ir['ratio']))
            f.write("\n")
    print("\n  Summary: {}".format(summary_path))

    # ============================================================
    # PLOTS
    # ============================================================
    print("  Generating plots...")
    colors = {30: '#27ae60', 210: '#e74c3c', 2310: '#2980b9', 30030: '#8e44ad'}

    # --- Figure 1: Gap Position (the main event) ---
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(
        'Phase Transition Hunt\n'
        r'Predicted $\alpha_c = \sqrt{{3/2}}$'
        ' = {:.6f}'.format(ALPHA_CRIT),
        fontsize=14, fontweight='bold', y=0.98)

    # Panel A: Gap position vs alpha (full view)
    ax = axes[0, 0]
    for m in avail:
        r = all_results[m]
        sym = "* " if r['I'] == 0 else ""
        ax.plot(r['alphas'], r['gap_pos'], color=colors[m], lw=1.5,
                label='{}m={} (I={}, phi={})'.format(sym, m, r['I'], r['phi']))
        if r['crash_alpha'] is not None:
            ax.axvline(r['crash_alpha'], color=colors[m], ls=':', alpha=0.4, lw=0.8)
    ax.axvline(ALPHA_CRIT, color='red', ls='--', alpha=0.8, lw=1.5,
               label=r'$\alpha_c=\sqrt{3/2}$')
    ax.axhline(0.5, color='gray', ls=':', alpha=0.4)
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel('Gap position (0=bottom, 1=top)')
    ax.set_title('A.  Gap Position — Full Sweep')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.25)
    ax.set_ylim(-0.05, 1.05)

    # Panel B: Gap position zoomed into transition
    ax = axes[0, 1]
    for m in avail:
        r = all_results[m]
        mask = (r['alphas'] >= 1.1) & (r['alphas'] <= 1.35)
        ax.plot(r['alphas'][mask], r['gap_pos'][mask], 'o-', color=colors[m],
                ms=2, lw=1.2, label='m={}'.format(m))
        if r['crash_alpha'] is not None:
            ax.axvline(r['crash_alpha'], color=colors[m], ls=':', alpha=0.5, lw=1)
    ax.axvline(ALPHA_CRIT, color='red', ls='--', alpha=0.8, lw=1.5,
               label=r'$\sqrt{3/2}$')
    ax.axhline(0.5, color='gray', ls=':', alpha=0.4)
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel('Gap position')
    ax.set_title('B.  Transition Zone [1.10, 1.35]')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.25)
    ax.set_ylim(-0.05, 1.05)

    # Panel C: Perron gap vs max internal gap (crossover detector)
    ax = axes[1, 0]
    for m in avail:
        r = all_results[m]
        mask = (r['alphas'] >= 1.0) & (r['alphas'] <= 1.4)
        ax.plot(r['alphas'][mask], r['perron_gap'][mask], '-', color=colors[m],
                lw=1.5, label='Perron gap (m={})'.format(m))
        ax.plot(r['alphas'][mask], r['max_int_gap'][mask], '--', color=colors[m],
                lw=1.0, alpha=0.7, label='Max bulk gap (m={})'.format(m))
    ax.axvline(ALPHA_CRIT, color='red', ls='--', alpha=0.8, lw=1.5)
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel('Gap size')
    ax.set_title('C.  Perron Gap vs Max Internal Gap')
    ax.legend(fontsize=6, ncol=2)
    ax.grid(True, alpha=0.25)

    # Panel D: Numerical derivative of gap position
    ax = axes[1, 1]
    for m in avail:
        r = all_results[m]
        mask = (r['alphas'] >= 1.05) & (r['alphas'] <= 1.40)
        a_masked = r['alphas'][mask]
        g_masked = r['gap_pos'][mask]
        if len(a_masked) > 1:
            derivs = np.diff(g_masked) / np.diff(a_masked)
            a_mid = 0.5 * (a_masked[:-1] + a_masked[1:])
            ax.plot(a_mid, derivs, color=colors[m], lw=1.2,
                    label='m={}'.format(m))
    ax.axvline(ALPHA_CRIT, color='red', ls='--', alpha=0.8, lw=1.5,
               label=r'$\sqrt{3/2}$')
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r"$d(\mathrm{gap\_pos})/d\alpha$")
    ax.set_title('D.  Derivative of Gap Position (crash detector)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.25)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    p1 = os.path.join(script_dir, 'phase_transition_hunt.png')
    plt.savefig(p1, dpi=150, bbox_inches='tight')
    print("  Saved: {}".format(p1))

    # --- Figure 2: IPR comparison ---
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    fig2.suptitle(r'IPR Ratio at Selected $\alpha$ Values', fontsize=13, fontweight='bold')
    for m in avail:
        r = all_results[m]
        if not r['ipr_results']:
            continue
        a_vals = sorted(r['ipr_results'].keys())
        ratios = [r['ipr_results'][a]['ratio'] for a in a_vals]
        ax2.plot(a_vals, ratios, 'o-', color=colors[m], ms=6, lw=1.5,
                 label='m={} (phi={})'.format(m, r['phi']))
    ax2.axvline(ALPHA_CRIT, color='red', ls='--', alpha=0.8, lw=1.5,
                label=r'$\sqrt{3/2}$')
    ax2.axvline(KAPPA, color='blue', ls='--', alpha=0.5, lw=1,
                label=r'$\sqrt{2/3}$')
    ax2.set_xlabel(r'$\alpha$')
    ax2.set_ylabel('Mean IPR / (1/phi)')
    ax2.set_title('Localization Across the Transition')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.25)
    plt.tight_layout()
    p2 = os.path.join(script_dir, 'phase_transition_ipr.png')
    plt.savefig(p2, dpi=150, bbox_inches='tight')
    print("  Saved: {}".format(p2))

    # --- Figure 3: Crash alpha convergence ---
    fig3, ax3 = plt.subplots(figsize=(10, 6))
    crash_ms = []
    crash_as = []
    steep_ms = []
    steep_as = []
    for m in avail:
        r = all_results[m]
        if r['crash_alpha'] is not None:
            crash_ms.append(m)
            crash_as.append(r['crash_alpha'])
        steep_ms.append(m)
        steep_as.append(r['steepest_alpha'])

    if crash_ms:
        ax3.semilogx(crash_ms, crash_as, 's-', color='#e74c3c', ms=8, lw=1.5,
                      label='Crash point (gap_pos=0.5)')
    ax3.semilogx(steep_ms, steep_as, 'D-', color='#2980b9', ms=8, lw=1.5,
                  label='Steepest descent')
    ax3.axhline(ALPHA_CRIT, color='red', ls='--', alpha=0.8, lw=1.5,
                label=r'$\sqrt{3/2}$ = ' + '{:.6f}'.format(ALPHA_CRIT))
    ax3.set_xlabel('Primorial m')
    ax3.set_ylabel(r'$\alpha$ at transition')
    ax3.set_title(r'Convergence of Transition Point to $\sqrt{3/2}$')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.25)
    plt.tight_layout()
    p3 = os.path.join(script_dir, 'phase_transition_convergence.png')
    plt.savefig(p3, dpi=150, bbox_inches='tight')
    print("  Saved: {}".format(p3))

    print()
    print("  === ANALYSIS COMPLETE ===")
    print("  Key question: Does crash_alpha converge to sqrt(3/2)?")
    print("  If errors shrink by ~10x per primorial level: YES")
    print("  If errors plateau or diverge: NO, look elsewhere")


if __name__ == '__main__':
    main()
