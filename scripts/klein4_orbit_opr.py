"""
Klein Four-Group Orbit Participation Ratio (OPR)
=================================================
Tests Gemini's prediction: the extreme-IPR flat-band modes at odd m
localize onto single Klein-4 orbits O(r) = {r, m-r, r^{-1}, (m-r)^{-1}}
under the joint action of σ (additive palindrome) and τ (multiplicative inverse).

Key predictions:
  - At m=1155 (odd): top-IPR mode has >90% mass in O(2) = {2, 577, 578, 1153}
  - At m=2310 (even): O(2) doesn't exist (gcd(2,2310)=2), so no single orbit dominates
  - The kurtosis drop 114 → 4 is the deletion of O(2) by the coprime sieve

Runtime: < 30s
"""

import numpy as np
from math import gcd, sqrt


def get_coprime_residues(m):
    return sorted([r for r in range(1, m) if gcd(r, m) == 1])


def build_D_sym(residues, m):
    r = np.array(residues)
    diff = np.abs(r[:, None] - r[None, :])
    return np.minimum(diff, m - diff).astype(float)


def build_P_tau(residues, m):
    n = len(residues)
    res_to_idx = {r: i for i, r in enumerate(residues)}
    perm = np.array([res_to_idx[pow(r, -1, m)] for r in residues])
    P = np.zeros((n, n))
    P[np.arange(n), perm] = 1.0
    return P


def build_klein4_orbits(residues, m):
    """Partition coprime residues into Klein-4 orbits under ⟨σ, τ⟩."""
    res_set = set(residues)
    res_to_idx = {r: i for i, r in enumerate(residues)}
    visited = set()
    orbits = []

    for r in residues:
        if r in visited:
            continue
        r_inv = pow(r, -1, m)
        mr = m - r
        mr_inv = m - r_inv
        # Intersect with coprime residues (all should be coprime, but be safe)
        orbit_residues = sorted({x for x in [r, mr, r_inv, mr_inv] if x in res_set})
        orbit_indices = [res_to_idx[x] for x in orbit_residues]
        orbits.append({
            'residues': orbit_residues,
            'indices': orbit_indices,
            'size': len(orbit_residues),
            'seed': min(orbit_residues),
        })
        visited.update(orbit_residues)

    return orbits


def ipr(v):
    p = np.abs(v) ** 2
    return np.sum(p ** 2)


def run_experiment(m, label):
    print(f"\n{'='*70}")
    print(f"  {label}: m = {m}")
    print(f"{'='*70}")

    residues = get_coprime_residues(m)
    phi = len(residues)
    res_to_idx = {r: i for i, r in enumerate(residues)}
    print(f"  φ = {phi}")

    # Build operators
    D_sym = build_D_sym(residues, m)
    P_tau = build_P_tau(residues, m)

    # Build Klein-4 orbits
    orbits = build_klein4_orbits(residues, m)
    orbit_sizes = [o['size'] for o in orbits]
    n_orbits = len(orbits)
    size_counts = {}
    for s in orbit_sizes:
        size_counts[s] = size_counts.get(s, 0) + 1
    print(f"\n  Klein-4 orbit decomposition:")
    print(f"    Total orbits: {n_orbits}")
    print(f"    Size distribution: {dict(sorted(size_counts.items()))}")
    print(f"    Total residues covered: {sum(orbit_sizes)}")

    # Show first few orbits (smallest seeds)
    print(f"\n    First 10 orbits (by seed):")
    for i, orb in enumerate(orbits[:10]):
        print(f"      O({orb['seed']:4d}) = {orb['residues']} (size {orb['size']})")

    # Check if O(2) exists
    o2_exists = any(o['seed'] == 2 for o in orbits)
    print(f"\n    O(2) orbit exists: {o2_exists}")
    if o2_exists:
        o2 = next(o for o in orbits if o['seed'] == 2)
        print(f"    O(2) = {o2['residues']}")
        # Verify the Klein-4 structure
        if len(o2['residues']) == 4:
            r = 2
            r_inv = pow(r, -1, m)
            print(f"    Verification: τ(2) = 2⁻¹ mod {m} = {r_inv}")
            print(f"                  σ(2) = {m} - 2 = {m-2}")
            print(f"                  σ(τ(2)) = {m} - {r_inv} = {m - r_inv}")

    # ==== EIGENDECOMPOSITION ====
    print(f"\n  ======= H(α_c) eigendecomposition =======")
    alpha_c = sqrt(135 / 88)
    lambda_P = np.max(np.linalg.eigvalsh(D_sym))
    C = D_sym @ P_tau - P_tau @ D_sym
    H = D_sym / lambda_P + 1j * alpha_c * C / lambda_P
    evals, evecs = np.linalg.eigh(H)

    # Flat band
    flat_mask = np.abs(evals) < 0.05
    n_flat = np.sum(flat_mask)
    print(f"  Flat band: {n_flat}/{phi}")

    # IPR for all eigenvectors
    ipr_vals = np.array([ipr(evecs[:, j]) * phi for j in range(phi)])

    # ==== ORBIT PARTICIPATION RATIO ====
    print(f"\n  ======= Orbit Participation Ratio (OPR) =======")

    # Get flat-band eigenvectors sorted by IPR
    flat_indices = np.where(flat_mask)[0]
    flat_ipr = [(ipr_vals[j], j) for j in flat_indices]
    flat_ipr.sort(reverse=True)

    print(f"\n  Top 20 most localized flat-band states:")
    print(f"  {'Rank':>4} {'IPR×φ':>8} {'MaxOrb%':>8} {'Seed':>6} {'OrbSz':>5} "
          f"{'Top2Orb%':>9} {'Residues in top orbit':>30}")

    results_top = []
    for rank in range(min(20, len(flat_ipr))):
        ipr_val, j = flat_ipr[rank]
        v = evecs[:, j]
        probs = np.abs(v) ** 2

        # Compute mass in each orbit
        orbit_masses = []
        for orb in orbits:
            mass = np.sum(probs[orb['indices']])
            orbit_masses.append(mass)
        orbit_masses = np.array(orbit_masses)

        # Sort orbits by mass
        mass_order = np.argsort(orbit_masses)[::-1]
        top_orb_idx = mass_order[0]
        top_orb = orbits[top_orb_idx]
        top_mass = orbit_masses[top_orb_idx]
        second_mass = orbit_masses[mass_order[1]] if len(mass_order) > 1 else 0
        top2_mass = top_mass + second_mass

        res_str = str(top_orb['residues'])
        if len(res_str) > 30:
            res_str = res_str[:27] + "..."

        print(f"  {rank+1:4d} {ipr_val:8.4f} {top_mass*100:7.1f}% "
              f"O({top_orb['seed']:4d}) {top_orb['size']:5d} "
              f"{top2_mass*100:8.1f}% {res_str:>30}")

        results_top.append({
            'rank': rank + 1,
            'ipr': ipr_val,
            'top_orbit_mass': top_mass,
            'top_orbit_seed': top_orb['seed'],
            'top_orbit_size': top_orb['size'],
            'top2_orbit_mass': top2_mass,
        })

    # ==== STATISTICAL SUMMARY ====
    print(f"\n  ======= OPR Statistics =======")

    # Top-5 average
    top5_masses = [r['top_orbit_mass'] for r in results_top[:5]]
    print(f"  Top-5 localized flat-band modes:")
    print(f"    Avg max-orbit mass: {np.mean(top5_masses)*100:.1f}%")
    print(f"    Min max-orbit mass: {np.min(top5_masses)*100:.1f}%")
    print(f"    Max max-orbit mass: {np.max(top5_masses)*100:.1f}%")

    # Is O(2) the dominant orbit?
    if o2_exists:
        top_seeds = [r['top_orbit_seed'] for r in results_top[:10]]
        o2_count = top_seeds.count(2)
        print(f"\n  O(2) dominance among top-10 localized modes: {o2_count}/10")

    # ==== ORBIT INVERSE PARTICIPATION RATIO ====
    # How many orbits does each mode span?
    print(f"\n  ======= Orbit-space IPR (how many orbits per mode) =======")
    print(f"  {'Rank':>4} {'IPR×φ':>8} {'OPR':>6} {'N_orb':>6} {'OPR/N':>6}")
    for rank in range(min(10, len(flat_ipr))):
        ipr_val, j = flat_ipr[rank]
        v = evecs[:, j]
        probs = np.abs(v) ** 2

        orbit_masses = np.array([np.sum(probs[orb['indices']]) for orb in orbits])
        # OPR = 1 / sum(w_k^2) where w_k = orbit mass fraction
        orbit_masses_norm = orbit_masses / orbit_masses.sum()
        opr = 1.0 / np.sum(orbit_masses_norm ** 2)

        print(f"  {rank+1:4d} {ipr_val:8.4f} {opr:6.2f} {n_orbits:6d} {opr/n_orbits:6.4f}")

    # ==== COMPARISON: ALL FLAT-BAND MODES ====
    print(f"\n  ======= All flat-band modes: orbit mass distribution =======")
    all_top_masses = []
    for ipr_val, j in flat_ipr:
        v = evecs[:, j]
        probs = np.abs(v) ** 2
        orbit_masses = np.array([np.sum(probs[orb['indices']]) for orb in orbits])
        all_top_masses.append(orbit_masses.max())

    all_top_masses = np.array(all_top_masses)
    from scipy.stats import kurtosis as scipy_kurtosis
    print(f"  Max-orbit mass across all {len(flat_ipr)} flat-band modes:")
    print(f"    Mean:   {all_top_masses.mean()*100:.2f}%")
    print(f"    Median: {np.median(all_top_masses)*100:.2f}%")
    print(f"    Max:    {all_top_masses.max()*100:.2f}%")
    print(f"    Std:    {all_top_masses.std()*100:.2f}%")
    # How many modes have >50% in one orbit?
    n_trapped = np.sum(all_top_masses > 0.5)
    n_strongly_trapped = np.sum(all_top_masses > 0.9)
    print(f"    Modes with >50% in single orbit: {n_trapped}/{len(flat_ipr)}")
    print(f"    Modes with >90% in single orbit: {n_strongly_trapped}/{len(flat_ipr)}")

    # ==== PERRON AND ANTI-PERRON IN ORBIT SPACE ====
    print(f"\n  ======= Perron/Anti-Perron orbit decomposition =======")
    idx_sorted = np.argsort(evals)
    for name, idx in [("Perron", idx_sorted[-1]), ("Anti-Perron", idx_sorted[0])]:
        v = evecs[:, idx]
        probs = np.abs(v) ** 2
        orbit_masses = np.array([np.sum(probs[orb['indices']]) for orb in orbits])
        mass_order = np.argsort(orbit_masses)[::-1]
        top3 = mass_order[:3]
        print(f"  {name} (λ = {evals[idx]:.6f}):")
        for k, oi in enumerate(top3):
            print(f"    #{k+1}: O({orbits[oi]['seed']:4d}) mass = {orbit_masses[oi]*100:.2f}% "
                  f"(size {orbits[oi]['size']})")
        opr = 1.0 / np.sum((orbit_masses / orbit_masses.sum()) ** 2)
        print(f"    OPR = {opr:.2f} / {n_orbits}")

    # ==== τ-PARITY CHECK (the test Gemini warned about) ====
    print(f"\n  ======= τ-parity check (Gemini's warning) =======")
    print(f"  Testing v @ P_tau @ v for top-10 flat-band modes:")
    print(f"  {'Rank':>4} {'IPR×φ':>8} {'v·Pτ·v':>10} {'|v·Pτ·v|':>10}")
    for rank in range(min(10, len(flat_ipr))):
        ipr_val, j = flat_ipr[rank]
        v = evecs[:, j]
        tau_parity = np.real(v.conj() @ P_tau @ v)
        print(f"  {rank+1:4d} {ipr_val:8.4f} {tau_parity:10.6f} {abs(tau_parity):10.6f}")

    return {
        'phi': phi,
        'n_orbits': n_orbits,
        'o2_exists': o2_exists,
        'top1_orbit_mass': results_top[0]['top_orbit_mass'] if results_top else 0,
        'top1_orbit_seed': results_top[0]['top_orbit_seed'] if results_top else None,
        'top5_avg_mass': np.mean(top5_masses),
        'n_trapped_50': n_trapped,
        'n_trapped_90': n_strongly_trapped,
        'kurtosis_all': scipy_kurtosis(ipr_vals, fisher=True),
    }


# ==============================================================
if __name__ == '__main__':
    results = {}
    results['odd'] = run_experiment(1155, "ODD: m = 1155")
    results['even'] = run_experiment(2310, "EVEN: m = 2310")

    # ==== SCORECARD ====
    print(f"\n{'='*70}")
    print(f"  SCORECARD: Klein Four-Group Orbit Participation")
    print(f"{'='*70}")

    r_odd = results['odd']
    r_even = results['even']

    print(f"\n  Prediction 1: O(2) orbit traps the top mode at odd m")
    if r_odd['o2_exists']:
        print(f"    O(2) exists at m=1155: YES")
        print(f"    Top mode mass in dominant orbit: {r_odd['top1_orbit_mass']*100:.1f}%")
        print(f"    Dominant orbit seed: O({r_odd['top1_orbit_seed']})")
        if r_odd['top1_orbit_seed'] == 2 and r_odd['top1_orbit_mass'] > 0.9:
            print(f"    ✓ CONFIRMED: >90% mass in O(2) — wave function collapses onto Klein-4 orbit")
        elif r_odd['top1_orbit_seed'] == 2 and r_odd['top1_orbit_mass'] > 0.5:
            print(f"    ~ PARTIAL: >50% in O(2) but not >90%")
        elif r_odd['top1_orbit_seed'] == 2:
            print(f"    ~ O(2) is dominant but weak ({r_odd['top1_orbit_mass']*100:.1f}%)")
        else:
            print(f"    ✗ O(2) is NOT the dominant orbit — O({r_odd['top1_orbit_seed']}) dominates")
    else:
        print(f"    O(2) does not exist at m=1155 (unexpected!)")

    print(f"\n  Prediction 2: O(2) deleted at even m")
    print(f"    O(2) exists at m=2310: {r_even['o2_exists']}")
    if not r_even['o2_exists']:
        print(f"    ✓ CONFIRMED: coprime sieve deletes O(2)")
        print(f"    Top mode mass at even m: {r_even['top1_orbit_mass']*100:.1f}% in O({r_even['top1_orbit_seed']})")
    else:
        print(f"    ✗ O(2) exists at even m (unexpected!)")

    print(f"\n  Prediction 3: Kurtosis drop = O(2) deletion")
    print(f"    Kurtosis: odd = {r_odd['kurtosis_all']:.1f}, even = {r_even['kurtosis_all']:.1f}")
    print(f"    Modes >90% in single orbit: odd = {r_odd['n_trapped_90']}, even = {r_even['n_trapped_90']}")
    print(f"    Modes >50% in single orbit: odd = {r_odd['n_trapped_50']}, even = {r_even['n_trapped_50']}")

    print(f"\n  Key question: Is the flat band a collection of decoupled Klein-4 rings?")
    if r_odd['n_trapped_90'] > 10:
        print(f"    → YES: {r_odd['n_trapped_90']} modes trapped in individual orbits")
    elif r_odd['n_trapped_50'] > 20:
        print(f"    → PARTIALLY: {r_odd['n_trapped_50']} modes with >50% orbit concentration")
    else:
        print(f"    → NO: modes are orbit-delocalized despite position-space hotspots")
    print()
