# Active Transport on the Prime Gas

**Flat-Band Condensation, Spectral Phase Transition, and the CRT Block / Qubit Architecture**

Antonio P. Matos — [ORCID 0009-0002-0722-3752](https://orcid.org/0009-0002-0722-3752)  
Fancyland LLC / Lattice OS  
March 2026 — Preprint

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19243258.svg)](https://doi.org/10.5281/zenodo.19243258)

---

## Abstract

We construct the Hermitian *active transport operator*

$$H(\alpha) = \frac{D_{\text{sym}}}{\lambda_P} + i\alpha \, \frac{[D_{\text{sym}},\, P_\tau]}{\lambda_P}$$

on the coprime residues modulo a primorial *m*, where $D_{\text{sym}}$ is the palindromic distance matrix, $P_\tau$ is the multiplicative inversion permutation, $\lambda_P$ is the Perron eigenvalue, and $\alpha \geq 0$ is a real shear parameter. Full eigendecomposition at four primorial levels ($m \in \{30, 210, 2310, 30030\}$; matrices up to $5760 \times 5760$) reveals:

1. **Flat-band condensation** — eigenvalue variance decays as $\text{Var}(\lambda) \approx 1.78 / \varphi(m)$
2. **Perron gap universality** — the Perron–anti-Perron gap locks to $\Delta_P \to 0.8146$
3. **Spectral phase transition** — a sharp gap-position crossover at critical $\alpha_c$
4. **$1/\sqrt{2}$ commutator ratio** — $\|[D, P_\tau]\|_F / \|D\|_F \to 1/\sqrt{2}$
5. **Sylow torsion defects** and the mod-4 partition
6. **CRT block architecture** — Chinese Remainder Theorem organizes the Hilbert space into 30 blocks of dimension 16
7. **Exact $\sigma$-parity conservation** — an arithmetic qubit whose states are $\tau$-parity within a $\sigma$-sector, protected by four layers of number-theoretic error resistance

The paper documents 8 falsified predictions (primorial qubit, CRT qudit, chirality balance, washing conjecture, and others), resolved by computational execution rather than model assertion.

---

## Repository Structure

```
├── paper/
│   ├── ACTIVE_TRANSPORT_PRIME_GAS.pdf    # Compiled manuscript (66 pages)
│   ├── ACTIVE_TRANSPORT_PRIME_GAS.tex    # LaTeX source
│   └── ACTIVE_TRANSPORT_PRIME_GAS.md     # Markdown source (approved draft)
│
├── scripts/                               # All 35 computational scripts cited in the paper
│   ├── quantum_prime_dynamics.py          # Core Hamiltonian constructor + α-sweep
│   ├── overnight_batch.py                 # Parallel batch driver for m = 30030
│   ├── cache_results.py                   # Result caching infrastructure
│   ├── data_run_ipr_gap.py               # IPR and gap measurements across α
│   ├── micro_sweep.py                     # Fine-resolution α sweep near α_c
│   ├── chirality_analysis.py              # Chirality identity chain (parts 1–3)
│   ├── chirality_analysis_part2.py
│   ├── chirality_analysis_part3.py
│   ├── chiral_phase_transition.py         # Phase transition characterization
│   ├── chiral_quench_test.py              # Loschmidt echo / quench dynamics
│   ├── partition_4class.py                # Mod-4 partition analysis
│   ├── verify_sylow_5adic.py              # 5-adic torsion defect verification
│   ├── verify_510510_zero_defect.py       # Zero-defect control at m = 510510
│   ├── sylow_defect_predictor.py          # Sylow defect prediction framework
│   ├── nozzle_test.py                     # Even/odd modulus nozzle test
│   ├── additive_ghost_extractor.py        # Additive ghost mode extraction
│   ├── flatband_hotspot_hunter.py         # Flat-band hotspot identification
│   ├── perturbation_robustness.py         # Perturbation robustness testing
│   ├── actual_eigenvector_test.py         # Eigenvector structure verification
│   ├── rabi_rabbit.py                     # ASDEX/ITER-style Rabi analogy
│   ├── primorial_logic_gate.py            # Primorial qubit gate attempt
│   ├── flatband_cavity_resist.py          # Flat-band cavity resistance
│   ├── corrected_rabi.py                  # Corrected Rabi analysis
│   ├── reconstructed_rabi.py              # Independent Rabi reconstruction
│   ├── true_qubit.py                      # CRT-block qubit search
│   ├── the_connection.py                  # V₂ = −I, V₃ = −½I proof
│   ├── crt_block_hamiltonian.py           # CRT block decomposition
│   ├── block9_deep_dive.py                # All-block survey
│   ├── flatband_qubit_search.py           # Screening hierarchy discovery
│   ├── find_connection.py                 # Algebraic link search
│   ├── invariant_subspace.py              # CRT block invariance testing
│   ├── crt_chirality_bridge.py            # CRT-chirality bridge
│   ├── orbit_hamiltonian.py               # Klein-4 orbit Hamiltonian
│   ├── klein4_orbit_opr.py                # Klein-4 orbit OPR analysis
│   └── error_protection_test.py           # Four-layer error protection test
│
├── LICENSE                                # MIT License
└── README.md
```

---

## Requirements

The scripts use standard scientific Python:

- **Python** ≥ 3.10
- **NumPy** — linear algebra and eigendecomposition
- **SciPy** — sparse matrices and special functions
- **Matplotlib** — plotting (optional, for visualization)
- **multiprocessing** — parallel α-sweeps (used by `overnight_batch.py`)

Install dependencies:

```bash
pip install numpy scipy matplotlib
```

---

## Running the Scripts

Each script is self-contained and can be run independently:

```bash
cd scripts
python quantum_prime_dynamics.py      # Core sweep (start here)
python the_connection.py              # V₂ = −I, V₃ = −½I proof
python error_protection_test.py       # Four-layer arithmetic error protection
```

The core script `quantum_prime_dynamics.py` constructs $H(\alpha)$ and sweeps over $\alpha$ values. Most scripts complete in under 15 minutes on a consumer PC for $m \leq 2310$. The full $m = 30030$ sweep ($5760 \times 5760$ matrices) takes approximately 60 minutes single-threaded; `overnight_batch.py` parallelizes this across available cores.

---

## Companion Papers

1. A. P. Matos, "The Universal Two-Prime Formula and Spectral Structure of Palindromic Distance Matrices," preprint (2026). DOI: [10.5281/zenodo.19210625](https://doi.org/10.5281/zenodo.19210625)
2. A. P. Matos, "Primorial Thermodynamics: Gauss Sums, the IPB98 Scaling Law, and the Hierarchy Break at m = 30030," preprint (2026). DOI: [10.5281/zenodo.19188924](https://doi.org/10.5281/zenodo.19188924)
3. A. P. Matos, "Spectral Isotropy and the Exact Temperature of the Prime Gas," preprint (2026). DOI: [10.5281/zenodo.19156532](https://doi.org/10.5281/zenodo.19156532)

---

## Citation

```bibtex
@article{matos2026activetransport,
  title   = {Active Transport on the Prime Gas: Flat-Band Condensation,
             Spectral Phase Transition, and the CRT Block Qubit Architecture},
  author  = {Matos, Antonio P.},
  year    = {2026},
  doi     = {10.5281/zenodo.19243258},
  note    = {Preprint}
}
```

---

## License

This project is licensed under the [MIT License](LICENSE).

© 2026 Fancyland LLC
