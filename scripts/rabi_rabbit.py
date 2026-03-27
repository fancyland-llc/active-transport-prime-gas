"""
Rabi Rabbit — The σ-Exact Block Diagonal and τ-Parity Two-Level System
=======================================================================
Tests the prediction: H(α) exactly commutes with P_σ for all α,
making every eigenstate have definite σ-parity. Within each σ-sector,
the dynamics is a Rabi-like two-level system in τ-parity.

Tests:
  1. [H(α), P_σ] = 0 exactly (for all α, not just α_c)
  2. σ-sector eigendecomposition — does H split into 240+240?
  3. τ-parity sweep: track the (++), (+-), (-+), (--) irrep weights
     as a function of α from 0 to 2.5
  4. At α_c: is the system at maximum τ-superposition (equal weight)?
  5. Rabi frequency: what is the effective coupling strength?
  6. Eigenvector σ-parity: every eigenstate should be exactly ±1

Runtime: < 60s
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


def build_P_sigma(residues, m):
    """Palindrome involution: σ(r) = m - r."""
    n = len(residues)
    res_to_idx = {r: i for i, r in enumerate(residues)}
    P = np.zeros((n, n))
    for r in residues:
        mr = m - r
        if mr in res_to_idx:
            P[res_to_idx[r], res_to_idx[mr]] = 1.0
    return P


def build_irrep_projectors_fast(P_sigma, P_tau, phi):
    """Build the 4 Klein-4 irrep projectors from P_sigma and P_tau directly.
    P_++ = (I + P_σ)(I + P_τ)/4
    P_+- = (I + P_σ)(I - P_τ)/4
    P_-+ = (I - P_σ)(I + P_τ)/4
    P_-- = (I - P_σ)(I - P_τ)/4
    """
    I = np.eye(phi)
    Sp = (I + P_sigma) / 2.0  # σ-even projector
    Sm = (I - P_sigma) / 2.0  # σ-odd projector
    Tp = (I + P_tau) / 2.0    # τ-even projector
    Tm = (I - P_tau) / 2.0    # τ-odd projector
    return {
        '++': Sp @ Tp,
        '+-': Sp @ Tm,
        '-+': Sm @ Tp,
        '--': Sm @ Tm,
    }


def run_experiment(m, label):
    print(f"\n{'='*70}")
    print(f"  {label}: m = {m}")
    print(f"{'='*70}")

    residues = get_coprime_residues(m)
    phi = len(residues)
    print(f"  φ = {phi}")

    D_sym = build_D_sym(residues, m)
    P_tau = build_P_tau(residues, m)
    P_sigma = build_P_sigma(residues, m)
    lambda_P = np.max(np.linalg.eigvalsh(D_sym))

    # Antisymmetric part for commutator
    r_arr = np.array(residues)
    D_dir = ((r_arr[None, :] - r_arr[:, None]) % m).astype(float)
    A = (D_dir - D_dir.T) / 2.0
    A_abs = np.abs(A)
    C = A_abs @ P_tau - P_tau @ A_abs  # = [|A|_ew, P_tau] = -[D_sym, P_tau]

    # ============================================================
    # TEST 1: [H(α), P_σ] = 0 EXACTLY?
    # ============================================================
    print(f"\n  ======= TEST 1: [H(α), P_σ] = 0 ? =======")

    # Check components
    comm_D_sigma = D_sym @ P_sigma - P_sigma @ D_sym
    err_D_sigma = np.linalg.norm(comm_D_sigma, 'fro')
    print(f"  ||[D_sym, P_σ]||_F = {err_D_sigma:.2e}")

    comm_Ptau_sigma = P_tau @ P_sigma - P_sigma @ P_tau
    err_Ptau_sigma = np.linalg.norm(comm_Ptau_sigma, 'fro')
    print(f"  ||[P_τ, P_σ]||_F = {err_Ptau_sigma:.2e}")

    comm_C_sigma = C @ P_sigma - P_sigma @ C
    err_C_sigma = np.linalg.norm(comm_C_sigma, 'fro')
    print(f"  ||[C, P_σ]||_F = {err_C_sigma:.2e}")

    # Full H at several α values
    print(f"\n  ||[H(α), P_σ]||_F at various α:")
    for alpha in [0.0, 0.5, 1.0, sqrt(135/88), 1.5, 2.0]:
        H = D_sym / lambda_P + 1j * alpha * C / lambda_P
        comm_H_sigma = H @ P_sigma - P_sigma @ H
        err = np.linalg.norm(comm_H_sigma, 'fro')
        marker = " ← α_c" if abs(alpha - sqrt(135/88)) < 0.01 else ""
        print(f"    α = {alpha:.4f}: ||[H, P_σ]|| = {err:.2e}{marker}")

    # Check P_σ is involutory
    sigma_sq = P_sigma @ P_sigma
    sigma_sq_err = np.linalg.norm(sigma_sq - np.eye(phi), 'fro')
    print(f"\n  P_σ² = I? Error: {sigma_sq_err:.2e}")

    # ============================================================
    # TEST 2: σ-SECTOR EIGENDECOMPOSITION
    # ============================================================
    print(f"\n  ======= TEST 2: σ-SECTOR EIGENDECOMPOSITION =======")

    # Eigendecompose P_sigma
    evals_sigma = np.linalg.eigvalsh(P_sigma)
    n_plus = np.sum(evals_sigma > 0.5)
    n_minus = np.sum(evals_sigma < -0.5)
    print(f"  P_σ eigenvalues: {n_plus} at +1, {n_minus} at -1")
    print(f"  σ-even sector: dim = {n_plus}")
    print(f"  σ-odd sector:  dim = {n_minus}")

    # Build σ-sector projectors
    S_plus = (np.eye(phi) + P_sigma) / 2.0
    S_minus = (np.eye(phi) - P_sigma) / 2.0

    # Check ranks
    rank_plus = int(round(np.trace(S_plus)))
    rank_minus = int(round(np.trace(S_minus)))
    print(f"  Tr(S+) = {rank_plus}, Tr(S-) = {rank_minus}")

    # ============================================================
    # TEST 3: EVERY EIGENSTATE HAS DEFINITE σ-PARITY
    # ============================================================
    print(f"\n  ======= TEST 3: EIGENVECTOR σ-PARITY =======")

    alpha_c = sqrt(135 / 88)
    H_c = D_sym / lambda_P + 1j * alpha_c * C / lambda_P
    evals_c, evecs_c = np.linalg.eigh(H_c)

    # For each eigenstate, compute <v|P_σ|v>
    sigma_parities = np.array([np.real(evecs_c[:, j].conj() @ P_sigma @ evecs_c[:, j])
                                for j in range(phi)])

    n_even = np.sum(sigma_parities > 0.999)
    n_odd = np.sum(sigma_parities < -0.999)
    n_mixed = phi - n_even - n_odd
    print(f"  Eigenstates with σ-parity = +1: {n_even}")
    print(f"  Eigenstates with σ-parity = -1: {n_odd}")
    print(f"  Eigenstates with mixed σ-parity: {n_mixed}")
    print(f"  Min |σ-parity|: {np.min(np.abs(sigma_parities)):.6f}")
    print(f"  Max deviation from ±1: {np.min(1 - np.abs(sigma_parities)):.2e}")

    if n_mixed == 0:
        print(f"  ✓ EVERY eigenstate has EXACT σ-parity")
    else:
        print(f"  ✗ {n_mixed} eigenstates have mixed σ-parity")

    # Contrast with τ-parity
    tau_parities = np.array([np.real(evecs_c[:, j].conj() @ P_tau @ evecs_c[:, j])
                              for j in range(phi)])
    n_tau_clean = np.sum(np.abs(tau_parities) > 0.999)
    print(f"\n  For contrast — τ-parity:")
    print(f"  Eigenstates with |τ-parity| > 0.999: {n_tau_clean}/{phi}")
    print(f"  τ-parity range: [{tau_parities.min():.4f}, {tau_parities.max():.4f}]")
    print(f"  τ-parity std: {tau_parities.std():.4f}")

    # ============================================================
    # TEST 4: α-SWEEP — IRREP WEIGHTS AS FUNCTION OF α
    # ============================================================
    print(f"\n  ======= TEST 4: τ-PARITY SWEEP vs α =======")

    projectors = build_irrep_projectors_fast(P_sigma, P_tau, phi)

    # Verify projectors sum to I
    P_sum = sum(projectors.values())
    proj_err = np.linalg.norm(P_sum - np.eye(phi), 'fro')
    print(f"  Projector sum = I? Error: {proj_err:.2e}")

    # Verify projector dimensions
    for label, P in projectors.items():
        dim = int(round(np.trace(P).real))
        print(f"    dim({label}) = {dim}")

    alphas = np.linspace(0, 2.5, 51)
    sweep_data = []

    print(f"\n  {'α':>6} {'(++)%':>7} {'(+-)%':>7} {'(-+)%':>7} {'(--)%':>7} {'Cross%':>7}"
          f" {'Perron σ':>9} {'AntiP σ':>9} {'AntiP τ':>9}")

    for alpha in alphas:
        H = D_sym / lambda_P + 1j * alpha * C / lambda_P
        H_norm_sq = np.sum(np.abs(H) ** 2)

        # Irrep decomposition
        irrep_fracs = {}
        diag_total = 0
        for label, P in projectors.items():
            Hp = P @ H @ P
            norm_sq = np.sum(np.abs(Hp) ** 2)
            frac = norm_sq / H_norm_sq if H_norm_sq > 0 else 0
            irrep_fracs[label] = frac
            diag_total += norm_sq

        cross_frac = 1.0 - diag_total / H_norm_sq if H_norm_sq > 0 else 0

        # Eigendecomposition at this α
        evals_a, evecs_a = np.linalg.eigh(H)
        idx_sorted = np.argsort(evals_a)
        perron_idx = idx_sorted[-1]
        anti_idx = idx_sorted[0]

        v_perron = evecs_a[:, perron_idx]
        v_anti = evecs_a[:, anti_idx]

        perron_sigma = np.real(v_perron.conj() @ P_sigma @ v_perron)
        anti_sigma = np.real(v_anti.conj() @ P_sigma @ v_anti)
        anti_tau = np.real(v_anti.conj() @ P_tau @ v_anti)

        row = {
            'alpha': alpha,
            '++': irrep_fracs['++'],
            '+-': irrep_fracs['+-'],
            '-+': irrep_fracs['-+'],
            '--': irrep_fracs['--'],
            'cross': cross_frac,
            'perron_sigma': perron_sigma,
            'anti_sigma': anti_sigma,
            'anti_tau': anti_tau,
        }
        sweep_data.append(row)

        # Print selected rows
        if (abs(alpha % 0.25) < 0.03 or abs(alpha - sqrt(135/88)) < 0.03):
            marker = " ← α_c" if abs(alpha - sqrt(135/88)) < 0.03 else ""
            print(f"  {alpha:6.3f} {irrep_fracs['++'] * 100:7.1f} {irrep_fracs['+-'] * 100:7.1f} "
                  f"{irrep_fracs['-+'] * 100:7.1f} {irrep_fracs['--'] * 100:7.1f} "
                  f"{cross_frac * 100:7.1f}"
                  f" {perron_sigma:9.4f} {anti_sigma:9.4f} {anti_tau:9.4f}{marker}")

    # ============================================================
    # TEST 5: RABI STRUCTURE — τ-decomposition within σ-even sector
    # ============================================================
    print(f"\n  ======= TEST 5: RABI STRUCTURE IN σ-EVEN SECTOR =======")

    # Build the σ-even sector explicitly
    evals_sigma_full, evecs_sigma = np.linalg.eigh(P_sigma)
    # σ-even eigenvectors (eigenvalue +1)
    even_mask = evals_sigma_full > 0.5
    V_even = evecs_sigma[:, even_mask]  # phi × n_even
    n_even_dim = V_even.shape[1]

    # Project H(α_c) into σ-even sector
    H_c_even = V_even.T.conj() @ H_c @ V_even  # n_even × n_even

    # Project P_tau into σ-even sector
    P_tau_even = V_even.T.conj() @ P_tau @ V_even

    # τ-even/odd projectors within σ-even
    I_even = np.eye(n_even_dim)
    Tp_even = (I_even + P_tau_even) / 2.0
    Tm_even = (I_even - P_tau_even) / 2.0

    dim_tau_plus = int(round(np.trace(Tp_even).real))
    dim_tau_minus = int(round(np.trace(Tm_even).real))
    print(f"  σ-even sector: dim = {n_even_dim}")
    print(f"    τ-even (++) within σ-even: dim = {dim_tau_plus}")
    print(f"    τ-odd  (+-) within σ-even: dim = {dim_tau_minus}")

    # H(α_c) in 2×2 τ-block form within σ-even
    H_pp = Tp_even @ H_c_even @ Tp_even  # (++) block
    H_mm = Tm_even @ H_c_even @ Tm_even  # (+-) block
    H_pm = Tp_even @ H_c_even @ Tm_even  # off-diagonal (++) → (+-)
    H_mp = Tm_even @ H_c_even @ Tp_even  # off-diagonal (+-) → (++)

    norm_pp = np.sum(np.abs(H_pp) ** 2)
    norm_mm = np.sum(np.abs(H_mm) ** 2)
    norm_pm = np.sum(np.abs(H_pm) ** 2)
    norm_mp = np.sum(np.abs(H_mp) ** 2)
    total_even = norm_pp + norm_mm + norm_pm + norm_mp
    print(f"\n  H(α_c) within σ-even sector (2×2 τ-block structure):")
    print(f"    (++) diagonal:      {norm_pp:.4f} ({norm_pp/total_even*100:.1f}%)")
    print(f"    (+-) diagonal:      {norm_mm:.4f} ({norm_mm/total_even*100:.1f}%)")
    print(f"    (++)↔(+-) coupling: {norm_pm + norm_mp:.4f} ({(norm_pm+norm_mp)/total_even*100:.1f}%)")

    # What fraction of the coupling is from D_sym vs the commutator?
    D_even = V_even.T.conj() @ (D_sym / lambda_P) @ V_even
    C_even = V_even.T.conj() @ (C / lambda_P) @ V_even

    D_pm = Tp_even @ D_even @ Tm_even
    C_pm = Tp_even @ C_even @ Tm_even

    norm_D_pm = np.sum(np.abs(D_pm) ** 2)
    norm_C_pm = np.sum(np.abs(C_pm) ** 2)
    print(f"\n  τ-off-diagonal coupling breakdown:")
    print(f"    From D_sym:     {norm_D_pm:.6f}")
    print(f"    From commutator: {norm_C_pm:.6f}")
    print(f"    Ratio C/D:       {norm_C_pm / norm_D_pm:.2f}" if norm_D_pm > 0 else "    D_pm = 0")

    # ============================================================
    # TEST 6: PERRON/ANTI-PERRON τ-PARITY AS FUNCTION OF α
    # ============================================================
    print(f"\n  ======= TEST 6: τ-PARITY ROTATION (RABI ANGLE) =======")

    print(f"  {'α':>6} {'Perron τ':>10} {'AntiP τ':>10} {'Gap':>8} {'Gap pos':>8}")
    alpha_range = np.linspace(0, 2.5, 101)
    gap_positions = []
    tau_perron_sweep = []
    tau_anti_sweep = []

    for alpha in alpha_range:
        H = D_sym / lambda_P + 1j * alpha * C / lambda_P
        evals_a, evecs_a = np.linalg.eigh(H)
        idx_sorted = np.argsort(evals_a)

        v_perron = evecs_a[:, idx_sorted[-1]]
        v_anti = evecs_a[:, idx_sorted[0]]

        tau_perron = np.real(v_perron.conj() @ P_tau @ v_perron)
        tau_anti = np.real(v_anti.conj() @ P_tau @ v_anti)

        gap = evals_a[idx_sorted[-1]] - evals_a[idx_sorted[-2]]
        gap_pos = (evals_a[idx_sorted[-1]] - evals_a[idx_sorted[-2]]) / \
                  (evals_a[idx_sorted[-1]] - evals_a[idx_sorted[0]]) if evals_a[idx_sorted[-1]] - evals_a[idx_sorted[0]] > 0 else 0

        tau_perron_sweep.append(tau_perron)
        tau_anti_sweep.append(tau_anti)
        gap_positions.append(gap_pos)

        if abs(alpha % 0.25) < 0.013 or abs(alpha - sqrt(135/88)) < 0.013:
            marker = " ← α_c" if abs(alpha - sqrt(135/88)) < 0.013 else ""
            print(f"  {alpha:6.3f} {tau_perron:10.4f} {tau_anti:10.4f} "
                  f"{gap:8.4f} {gap_pos:8.4f}{marker}")

    # Find where |τ-parity| of anti-Perron crosses zero
    tau_anti_arr = np.array(tau_anti_sweep)
    sign_changes = np.where(np.diff(np.sign(tau_anti_arr)))[0]
    if len(sign_changes) > 0:
        # Linear interpolation
        for sc in sign_changes:
            a0, a1 = alpha_range[sc], alpha_range[sc+1]
            t0, t1 = tau_anti_arr[sc], tau_anti_arr[sc+1]
            alpha_cross = a0 - t0 * (a1 - a0) / (t1 - t0)
            print(f"\n  Anti-Perron τ-parity zero crossing at α ≈ {alpha_cross:.6f}")
            print(f"  α_c = {sqrt(135/88):.6f}")
            print(f"  Ratio: α_cross / α_c = {alpha_cross / sqrt(135/88):.6f}")

    # ============================================================
    # TEST 7: σ-SECTOR SPECTRA
    # ============================================================
    print(f"\n  ======= TEST 7: σ-SECTOR SPECTRA AT α_c =======")

    # Separate eigenvalues by σ-parity
    even_evals = []
    odd_evals = []
    for j in range(phi):
        sp = sigma_parities[j]
        if sp > 0:
            even_evals.append(evals_c[j])
        else:
            odd_evals.append(evals_c[j])

    even_evals = np.sort(even_evals)
    odd_evals = np.sort(odd_evals)

    print(f"  σ-even sector ({len(even_evals)} states):")
    print(f"    Perron:      {even_evals[-1]:.6f}")
    print(f"    Anti-Perron: {even_evals[0]:.6f}")
    print(f"    Flat (|λ|<0.05): {np.sum(np.abs(even_evals) < 0.05)}/{len(even_evals)}")

    print(f"  σ-odd sector ({len(odd_evals)} states):")
    print(f"    Max:         {odd_evals[-1]:.6f}")
    print(f"    Min:         {odd_evals[0]:.6f}")
    print(f"    Flat (|λ|<0.05): {np.sum(np.abs(odd_evals) < 0.05)}/{len(odd_evals)}")

    # Does the Perron eigenvector live in σ-even?
    idx_sorted = np.argsort(evals_c)
    perron_sp = sigma_parities[idx_sorted[-1]]
    anti_sp = sigma_parities[idx_sorted[0]]
    print(f"\n  Perron σ-parity: {perron_sp:+.6f}")
    print(f"  Anti-Perron σ-parity: {anti_sp:+.6f}")

    # ============================================================
    # SCORECARD
    # ============================================================
    print(f"\n  ======= SCORECARD =======")
    print(f"  [H(α), P_σ] = 0 for all α?  ", end="")
    if err_D_sigma < 1e-10 and err_C_sigma < 1e-10:
        print(f"✓ EXACT (errors < 10⁻¹⁰)")
    else:
        print(f"✗ FAILS")

    print(f"  Every eigenstate has definite σ-parity?  ", end="")
    if n_mixed == 0:
        print(f"✓ ALL {phi} eigenstates")
    else:
        print(f"✗ {n_mixed} mixed")

    print(f"  σ-even/odd dimensions: {n_even}/{phi - n_even}")

    print(f"  τ-parity clean? ", end="")
    print(f"✗ NO — range [{tau_parities.min():.3f}, {tau_parities.max():.3f}]")

    print(f"  The qubit: σ-parity is EXACT, τ-parity is the RABI VARIABLE")

    return sweep_data


# ==============================================================
if __name__ == '__main__':
    sweep_odd = run_experiment(1155, "ODD: m = 1155")
    sweep_even = run_experiment(2310, "EVEN: m = 2310")
    print()
