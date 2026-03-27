# Active Transport on the Prime Gas: Flat-Band Condensation, the Rabi Phase Transition, and the Arithmetic Qubit

**Author:** Antonio P. Matos  
**ORCID:** [0009-0002-0722-3752](https://orcid.org/0009-0002-0722-3752)  
**Date:** March 27, 2026  
**Affiliation:** Independent Researcher; Fancyland LLC / Lattice OS  
**Status:** Preprint  
**DOI:** [10.5281/zenodo.19243258](https://doi.org/10.5281/zenodo.19243258)  

**MSC 2020 (Primary):** 11A25 (Arithmetic functions; other than those in 11Nxx), 15A18 (Eigenvalues, singular values, and eigenvectors)  
**MSC 2020 (Secondary):** 11A07 (Congruences; primitive roots; residue systems), 81Q10 (Selfadjoint operator theory in quantum theory)  

**Keywords:** active transport operator, palindromic distance matrix, flat-band condensation, spectral phase transition, Rabi crossing, σ-parity conservation, τ-parity, arithmetic qubit, commutator ratio, prime gas, coprime residues, Perron eigenvalue, inverse participation ratio, primorial hierarchy, selective screening, Loschmidt echo, CRT block decomposition, Klein four-group, Sylow torsion defect

**Companion papers:**

- [1] "The Universal Two-Prime Formula and Spectral Structure of Palindromic Distance Matrices," Matos (2026), DOI: [10.5281/zenodo.19210625](https://doi.org/10.5281/zenodo.19210625)
- [2] "Primorial Thermodynamics: Gauss Sums, the IPB98 Scaling Law, and the Hierarchy Break at $m = 30030$," Matos (2026), DOI: [10.5281/zenodo.19188924](https://doi.org/10.5281/zenodo.19188924)
- [3] "Spectral Isotropy and the Exact Temperature of the Prime Gas," Matos (2026), DOI: [10.5281/zenodo.19156532](https://doi.org/10.5281/zenodo.19156532)

---

## Abstract

We construct the Hermitian *active transport operator*

$$H(\alpha) = \frac{D_{\text{sym}}}{\lambda_P} + i\alpha \, \frac{[D_{\text{sym}},\, P_\tau]}{\lambda_P}$$

on the coprime residues modulo a primorial $m$, where $D_{\text{sym}}$ is the palindromic distance matrix, $P_\tau$ is the multiplicative inversion permutation, $\lambda_P$ is the Perron eigenvalue, and $\alpha \geq 0$ is a real shear parameter. Because $[D_{\text{sym}}, P_\tau]$ is real and skew-symmetric, the operator $H(\alpha) = A + i\alpha B$ with $A$ real symmetric and $B$ real skew-symmetric is Hermitian with purely real spectrum for every $\alpha$.

Full eigendecomposition of $H(\alpha)$ at the first four primorial levels $m \in \{30, 210, 2310, 30030\}$ (matrices up to $5760 \times 5760$; 62 $\alpha$-values per modulus; 19 minutes total on a consumer PC using parallelized sweeps; the single-threaded $m = 30030$ sweep alone takes $\sim 60$ minutes) reveals seven new spectral phenomena:

**(1) Flat-band condensation.** The normalized eigenvalue variance of $H(\alpha)$ at any fixed $\alpha > 0$ decays as $\text{Var}(\lambda) \approx 1.78/\varphi(m)$:

| $m$ | $\varphi(m)$ | $\text{Var}(\lambda)$ at $\alpha = \sqrt{2/3}$ | $\text{Var} \cdot \varphi$ |
|---|---|---|---|
| 30 | 8 | 0.203022 | 1.624 |
| 210 | 48 | 0.036839 | 1.768 |
| 2310 | 480 | 0.003702 | 1.777 |
| 30030 | 5760 | 0.000309 | 1.780 |

In the primorial limit, the spectrum collapses to a measure supported on three isolated points — the Perron eigenvalue near $+1$, an anti-Perron mode near $-0.744$ (at $\alpha_c$), and a degenerate zero mode carrying $\sim 99.74\%$ of the spectral weight (5745/5760 states at $m = 30030$) — constituting a *flat band* in the sense of [8].

**(2) Perron gap universality.** The gap between the Perron eigenvalue and the second-largest eigenvalue of $H(\alpha)$, evaluated at $\alpha = \sqrt{2/3}$, converges rapidly to a universal constant:

$$\Delta_P = \lambda_1 - \lambda_2 \longrightarrow 0.81464 \pm 0.00002$$

with the sequence $0.8568 \to 0.8138 \to 0.8146 \to 0.8146$ across $m = 30, 210, 2310, 30030$. The spectral width simultaneously converges to $1.5908 \pm 0.0001$; no closed-form candidate for this constant has been identified (Problem 3).

**(3) Gap position phase transition.** The position of the largest consecutive eigenvalue gap within the sorted spectrum undergoes a sharp transition from $1.0$ (Perron-dominated) to $0.0$ (bulk-dominated) at a critical shear amplitude $\alpha_c$. A fine-resolution sweep (222 $\alpha$-values in $[0.5, 1.45]$) initially suggested convergence to $\alpha_c \approx 1.2386$ with an apparent plateau between $m = 2310$ and $m = 30030$. A subsequent micro-sweep (2003 $\alpha$-values, spacing $7.5 \times 10^{-6}$) revealed this plateau was a **resolution artifact** — the transition width ($\sim 4 \times 10^{-6}$) was 500 times narrower than the coarse grid spacing. The true crash at $m = 2310$ is $\alpha_c = 1.2389219$, with derivative $-133{,}000$. The $m = 30030$ micro-sweep (503 phase-B points) resolves the crash at $\alpha_c = 1.2387666$ with width $1.0 \times 10^{-7}$ and derivative $-4{,}989{,}133$ — a $38\times$ narrower transition and $37\times$ steeper crash than $m = 2310$. A non-primorial control ($m = 385$, $\varphi = 240$) confirms the transition is universal.

The candidate closed form is $\alpha_c = \sqrt{135/88} = \sqrt{3^3 \cdot 5 / (2^3 \cdot 11)} \approx 1.23858$ (Conjecture 6.1):

| $m$ | Type | $\varphi(m)$ | Crash $\alpha$ | Error from $\sqrt{135/88}$ | Width | Steepest $dg/d\alpha$ |
|---|---|---|---|---|---|---|
| 210 | primorial | 48 | 1.27025 | $+2.56\%$ | — | — |
| 385 | non-primorial | 240 | 1.24227 | $+0.30\%$ | $6.0 \times 10^{-6}$ | $-82{,}817$ |
| 2310 | primorial | 480 | 1.23892 | $+0.027\%$ | $3.8 \times 10^{-6}$ | $-133{,}000$ |
| 15015 | non-primorial | 5760 | 1.23877 | $+0.015\%$ | $8.0 \times 10^{-6}$ | $-62{,}364$ |
| 30030 | primorial | 5760 | 1.23877 | $+0.015\%$ | $1.0 \times 10^{-7}$ | $-4{,}989{,}133$ |

Convergence toward $\sqrt{135/88}$ continues monotonically but with severe deceleration: the local scaling exponent drops from $\beta \approx 2.0$ ($\varphi = 48 \to 480$) to $\beta \approx 0.24$ ($\varphi = 480 \to 5760$), inconsistent with any stable power law and suggestive of logarithmic or sub-algebraic convergence. **Torsion width effect:** The non-primorial $m = 15015$ ($\varphi = 5760$) crashes at effectively the same $\alpha_c$ as $m = 30030$ ($\delta = +0.00018\%$), but with an $80\times$ wider transition and $80\times$ shallower derivative. Sylow torsion does not shift the crash — it **smears** the transition.

**(4) Sylow torsion defects.** The mod-4 partition used to define odd/even blocks is **not** a multiplicative character on $(\mathbb{Z}/m\mathbb{Z})^*$ when $4 \nmid m$ (which holds for all primorials). This non-multiplicativity means $q$-torsion defects — asymmetries in how residues with $q$-adic multiplicative order distribute across the mod-4 blocks — are not washed out by primorial expansion, contrary to a CRT-based algebraic prediction. At $m = 2310$, the 5-adic defect $\Delta v_5 = 4$ (from $11 \equiv 1 \pmod{5}$) is confirmed by exact computation. At $m = 510510$, residual defects $\Delta v_3 = \Delta v_5 = 16$ persist despite the terminal prime 17 injecting no torsion, falsifying the "Markovian washing" hypothesis (Section 8).

**(5) Selective screening of prime-periodic perturbations.** Treating $H(\alpha)$ as a tight-binding Hamiltonian and measuring the time-averaged Loschmidt echo $\langle L \rangle_\infty$ under perturbation $H_\varepsilon = H(\alpha) + \varepsilon V$, we find that the flat band is fragile against random Hermitian noise ($35.6\%$ survival at $\varepsilon = 0.03$) but selectively robust against prime-periodic diagonal potentials $V_p = \text{diag}(\cos(2\pi r/p))$. At $m = 2310$ and $\alpha = \alpha_c$, in-factorization primes ($p \mid m$) retain $95.0\%$ of the echo at $\varepsilon = 1.0$, versus $70.0\%$ for missing primes ($p \nmid m$) and $1.4\%$ for random noise. The primes $p = 2$ and $p = 3$ are perfectly screened ($100\%$ survival at all $\varepsilon$) because the coprime sieve eliminates those Fourier components from the Hilbert space. The IN/OUT screening ratio amplifies from $1.15$ at $\alpha = 0$ to $1.36$ at $\alpha_c$, establishing that the critical shear enhances arithmetic discrimination (Section 10.1).

**(6) CRT block architecture.** The cosine perturbation operators $V_p = \text{diag}(\cos(2\pi r/p))$ satisfy $V_2 = -I$ and $V_3 = -\tfrac{1}{2}I$ on the coprime residues — they act as exact scalars, a consequence of $\cos(2\pi r/p)$ being constant on coprimes only for $p = 2, 3$ (the *$p = 3$ boundary*). The remaining in-factorization operators $V_5, V_7, V_{11}$ have non-degenerate spectra whose joint eigenspaces define a CRT decomposition of the 480-dimensional Hilbert space into $30$ blocks of dimension $16$. Within each block, $V_5$ spread $= 3.2 \times 10^{-13}$ (numerical zero, eigensolver tolerance) — the protection is exact by construction, not approximate. The missing-prime operators $V_{13}$ (noise $= 0.008$), $V_{17}$ ($0.014$), $V_{19}$ ($0.007$) act as non-trivial gates within each 16D block, with maximum Rabi transfer $0.41$. Because the CRT blocking operators $V_p = \text{diag}(\cos(2\pi r/p))$ satisfy $\cos(2\pi(m-r)/p) = \cos(2\pi r/p)$ for $p \mid m$, they commute exactly with $P_\sigma$, so each 16D block bisects into two 8D halves under $\sigma$-parity. Within the active $\sigma$-even register, the CRT sublevels are 8-dimensional, giving an **8-level qudit** per qubit state. Block leakage under $H(\alpha_c)$ evolution is $7.4\%$ at $t = 200$ ($92.6\%$ contained). The CRT blocks are not the computational basis — they are the Zeeman fine structure within each qubit level (Section 11).

**(7) Exact $\sigma$-parity conservation and the arithmetic qubit.** The palindromic involution $P_\sigma: r \mapsto m - r$ commutes exactly with $H(\alpha)$ for all $\alpha$: $[H(\alpha), P_\sigma] = 0$ (algebraic identity, not numerical approximation; error $= 0.00 \times 10^{0}$ at all tested moduli and $\alpha$-values). Consequently, every eigenstate has exact $\sigma$-parity $\pm 1$, and the Hilbert space splits into two impenetrable $\varphi(m)/2$-dimensional sectors ($\mathcal{H}_+$ and $\mathcal{H}_-$). Within each sector, the Hamiltonian acts as a two-level system coupling $\tau$-even and $\tau$-odd subspaces through the active transport commutator. The phase transition at $\alpha_c$ is a Rabi crossing — the resonant balance point where the $\tau$-off-diagonal coupling strength (cross-irrep Frobenius share $= 50.7\%$ of the Hamiltonian at $\alpha_c$) matches the $\tau$-diagonal detuning, producing maximum mixing between the $\tau$-even (additive standing wave) and $\tau$-odd (multiplicative current) characters. The commutator-to-$D_{\text{sym}}$ ratio in $\tau$-off-diagonal coupling is exactly $4.00$. This establishes a quantum number hierarchy: $\sigma$-parity (exact, conserved) $\supset$ $\tau$-parity (dynamical Rabi variable) $\perp$ CRT block (spectral fine structure, transverse to $\tau$). The system is an arithmetic qubit with four layers of error protection: algebraic $\sigma$-conservation, number-theoretic coprime sieve, spectral flat-band averaging, and dimensional collective encoding (Section 11).

As a consistency check, the computation independently recovers the $1/\sqrt{2}$ commutator ratio:

$$\frac{\|[D_{\text{sym}},\, P_\tau]\|_F}{\|D_{\text{sym}}\|_F} \longrightarrow \frac{1}{\sqrt{2}}$$

which is a derived consequence of two results from the companion paper [1]: the irreconcilability coupling constant $\|[D,\tau]\|_F / \lambda_P \to \sqrt{2/3}$ and the Frobenius-to-Perron ratio $\|D\|_F / \lambda_P \to 2/\sqrt{3}$. The identity $1/\sqrt{2} = \sqrt{2/3} \,/\, (2/\sqrt{3})$ holds exactly.

The computation also falsifies three predictions (Section 7): (i) no spectral phase transition occurs at $\alpha = \sqrt{2/3}$; (ii) chirality balance ($|n_+ - n_-| = 0$) does not control localization — the balanced modulus $m = 30$ is the *most* dispersive, not the most localized; and (iii) the conjectured exact critical value $\alpha_c = \sqrt{3/2} = 1/\sqrt{2/3}$ is ruled out by micro-sweep data showing the crash at $+1.16\%$ above $\sqrt{3/2}$, converging instead toward $\sqrt{135/88}$.

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Construction of the Active Transport Operator](#2-construction-of-the-active-transport-operator)
3. [The $1/\sqrt{2}$ Commutator Ratio](#3-the-1sqrt2-commutator-ratio)
4. [Flat-Band Condensation](#4-flat-band-condensation)
5. [Perron Gap Universality](#5-perron-gap-universality)
6. [The Gap Position Phase Transition](#6-the-gap-position-phase-transition)
7. [Falsified Predictions](#7-falsified-predictions)
8. [Sylow Torsion Defects and the Mod-4 Partition](#8-sylow-torsion-defects-and-the-mod-4-partition)
9. [Discussion](#9-discussion)
10. [Perturbation Response and the CRT Block Architecture](#10-perturbation-response-and-the-crt-block-architecture)
11. [Exact σ-Parity Conservation and the Arithmetic Qubit](#11-exact-σ-parity-conservation-and-the-arithmetic-qubit)
12. [Open Problems](#12-open-problems)
13. [Appendix A: Computational Details](#13-appendix-a-computational-details)
14. [Script Inventory](#14-script-inventory)
15. [Summary](#15-summary)
16. [References](#16-references)

---

## 1. Introduction

### 1.1 Background and Motivation

The palindromic distance matrix $D_{\text{sym}}$ on coprime residues modulo a primorial $m = \prod_{i=1}^k p_i$ [4] is a real symmetric non-negative matrix of dimension $\varphi(m) \times \varphi(m)$, with entries $D_{ij} = \min(|r_i - r_j|, m - |r_i - r_j|)$ for coprime residues $r_i$. Its spectral properties, proved in [1] and [2], include: (i) a polynomial hierarchy of determinant ratios terminating at $m = 30030$ [2]; (ii) a Perron eigenvalue growing as $\tfrac{1}{2}\prod p(p-1)$ with base constant $C = 1/2$ [1], consistent with the Perron–Frobenius theorem for non-negative matrices [5]; (iii) a boundary saturation phenomenon at $-1/\sqrt{2}$ [1]; and (iv) the irreconcilability coupling constant $\|[D,\tau]\|_F / \lambda_P \to \sqrt{2/3}$ [1].

The underlying mechanism is the *Permanent Irreconcilability Theorem* [1]: the additive palindromic involution $\sigma(r) = m - r$ and the multiplicative inverse $\tau(r) = r^{-1} \bmod m$ can never agree on any coprime residue when $3 \mid m$, because $r^2 \equiv -1 \pmod{m}$ has no solution. The commutator $[D_{\text{sym}}, P_\tau]$ measures this irreducible tension.

The *Double Helix Parity Law* [1] further reveals that the odd prime factors of $m$ organize into two chiral strands indexed by their quadratic residue class modulo 4: primes $p \equiv 1 \pmod{4}$ (Strand $+$) and primes $p \equiv 3 \pmod{4}$ (Strand $-$). In the FFT power spectrum of the eigenvalue ratio, same-parity ancestral shadows dominate opposite-parity shadows by factors of $2.3\times$ to $17.5\times$.

### 1.2 The Chirality Balance Hypothesis

These structural results motivate a natural hypothesis:

> *When the chirality imbalance $I = |n_+ - n_-|$ vanishes, the multiplicative inverse operator should perfectly annihilate additive spatial turbulence, producing a topological phase transition at shear amplitude $\alpha = \sqrt{2/3}$.*

Specifically, the hypothesis predicts:

1. A *spectral quench* — block-diagonalization or flat-band collapse — at $\alpha = \sqrt{2/3}$ for chirally balanced moduli ($I = 0$).
2. The balanced modulus $m = 30$ ($I = 0$) should exhibit *stronger* localization than unbalanced moduli ($m = 210, 2310, 30030$, all with $I \geq 1$).
3. The chirality table:

| $m$ | Active odd primes | $p \equiv 1 \pmod{4}$ | $p \equiv 3 \pmod{4}$ | $I$ |
|---|---|---|---|---|
| 30 | $\{3, 5\}$ | $\{5\}$ | $\{3\}$ | 0 |
| 210 | $\{3, 5, 7\}$ | $\{5\}$ | $\{3, 7\}$ | 1 |
| 2310 | $\{3, 5, 7, 11\}$ | $\{5\}$ | $\{3, 7, 11\}$ | 2 |
| 30030 | $\{3, 5, 7, 11, 13\}$ | $\{5, 13\}$ | $\{3, 7, 11\}$ | 1 |

The hypothesis is computationally falsifiable: construct $H(\alpha)$, sweep $\alpha$ from 0 to 1.5, and search for a spectral reorganization at $\alpha = \sqrt{2/3} \approx 0.8165$.

### 1.3 What the Computation Found

The test was executed for all four primorials, requiring full eigendecomposition of Hermitian matrices up to $5760 \times 5760$ at 62 shear values each. The results:

1. **No phase transition at $\alpha = \sqrt{2/3}$.** All four spectral indicators (block separation, gap ratio, eigenvalue variance, gap position) pass smoothly through $\alpha = 0.8165$ without any sharp feature.

2. **Chirality balance does not control localization.** The balanced modulus $m = 30$ ($I = 0$) is the *most dispersive* (highest eigenvalue variance, lowest IPR ratio), not the most localized. Localization scales with $\varphi(m)$, not chirality.

3. **The actual transition lives at $\alpha_c \approx \sqrt{135/88} \approx 1.2386$.** The gap position undergoes a sharp crash from 1.0 to 0.0 at $\alpha_c$, where the Perron gap is overtaken by an internal bulk gap. A micro-sweep at $m = 2310$ resolves the crash at $\alpha_c = 1.2389219$ with width $\sim 4 \times 10^{-6}$. The $m = 30030$ micro-sweep (503 phase-B points) sharpens this to $\alpha_c = 1.2387666$ with width $1.0 \times 10^{-7}$ and derivative $-4{,}989{,}133$. The candidate closed form $\sqrt{135/88}$ matches to $+0.015\%$ error, with convergence decelerating sharply (local exponent drops from $\sim 2.0$ to $\sim 0.24$ across successive primorial pairs). A non-primorial control ($m = 385$) confirms the transition is universal (Section 6).

4. **Three new spectral universalities emerge:** flat-band condensation ($\text{Var} \propto 1/\varphi$), Perron gap convergence ($\Delta_P \to 0.8146$), and the $1/\sqrt{2}$ commutator ratio — the last being a derived verification of the companion paper's results.

### 1.4 Overview

Section 2 constructs the active transport operator $H(\alpha)$ and establishes its Hermiticity, including the derivation of the chirality identity chain from Spectral Isotropy (Section 2.5). Section 3 derives the $1/\sqrt{2}$ commutator ratio from known results and verifies it at four primorial levels. Section 4 proves flat-band condensation and presents the torsion crossover mechanism via isospectral non-conjugacy (Section 4.5). Section 5 establishes Perron gap universality. Section 6 identifies the gap position phase transition. Section 7 documents the falsified predictions. Section 8 discusses Sylow torsion defects and the mod-4 partition. Section 9 collects spectral architecture at $\alpha_c$, the Landau-Zener connection, the nozzle test, and plasma incompatibility. Section 10 presents the perturbation response experiments and the CRT block architecture. Section 11 presents the exact $\sigma$-parity conservation theorem that reveals the true computational basis: an arithmetic qubit whose states are $\tau$-parity within a $\sigma$-sector, protected by four layers of number-theoretic error resistance. Section 12 collects open problems.

**Guide for different readers.** Readers interested in the algebraic structure should focus on Sections 2–3, 8, and 11. Those interested in spectral phenomenology should prioritize Sections 4–6 and 9. The perturbation response and CRT architecture (Section 10.1–10.3) are largely self-contained. The central result — the arithmetic qubit and its four-layer protection — is in Section 11. Open problems (Section 12) are grouped by theme: spectral constants, torsion, perturbation response, chirality, and the $\sigma$-$\tau$ architecture. Reproducibility details are in Sections 13–14.

---

## 2. Construction of the Active Transport Operator

### 2.1 The Palindromic Distance Matrix

For a primorial $m = \prod_{i=1}^k p_i$ with coprime residues $\mathcal{R} = \{r : 1 \leq r < m, \gcd(r,m) = 1\}$ of cardinality $n = \varphi(m)$, the palindromic distance matrix $D_{\text{sym}} \in \mathbb{R}^{n \times n}$ has entries

$$D_{ij} = \min(|r_i - r_j|, \, m - |r_i - r_j|)$$

for residues $r_i, r_j \in \mathcal{R}$, ordered by magnitude. The matrix is real, symmetric, non-negative, with zeros on the diagonal.

### 2.2 The Multiplicative Inversion Permutation

The multiplicative inverse map $\tau : r \mapsto r^{-1} \bmod m$ is a bijection on $\mathcal{R}$ (since $\gcd(r,m) = 1$ implies $r^{-1}$ exists and $\gcd(r^{-1},m) = 1$). Its matrix representation $P_\tau \in \{0,1\}^{n \times n}$ is a permutation matrix: $(P_\tau)_{ij} = \mathbf{1}[r_j \equiv r_i^{-1} \pmod{m}]$.

Key properties:
- $P_\tau^2 = I$ (since $\tau(\tau(r)) = r$), so $P_\tau$ is an involution.
- $P_\tau = P_\tau^\top$ (the permutation matrix of an involution is symmetric).
- $P_\tau$ is orthogonal: $P_\tau P_\tau^\top = I$.

### 2.3 The Commutator and Its Skew-Symmetry

The commutator $[D_{\text{sym}}, P_\tau] = D_{\text{sym}} P_\tau - P_\tau D_{\text{sym}}$ has a direct entry-wise formula: since $(DP_\tau)_{ij} = D_{i,\tau(j)}$ (permuting columns) and $(P_\tau D)_{ij} = D_{\tau(i),j}$ (permuting rows),

$$[D, P_\tau]_{ij} = d(r_i, r_{\tau(j)}) - d(r_{\tau(i)}, r_j).$$

**Proposition 2.1** (Skew-symmetry). $[D_{\text{sym}}, P_\tau]^\top = -[D_{\text{sym}}, P_\tau]$.

*Proof.* $[D, P_\tau]_{ji} = d(r_j, r_{\tau(i)}) - d(r_{\tau(j)}, r_i) = -(d(r_{\tau(j)}, r_i) - d(r_j, r_{\tau(i)}))$. Since palindromic distance is symmetric in both arguments ($d(a,b) = d(b,a)$) and relabeling gives $d(r_{\tau(j)}, r_i) = d(r_i, r_{\tau(j)})$ and $d(r_j, r_{\tau(i)}) = d(r_{\tau(i)}, r_j)$, we have $[D, P_\tau]_{ji} = -(d(r_i, r_{\tau(j)}) - d(r_{\tau(i)}, r_j)) = -[D, P_\tau]_{ij}$. $\square$

**Verification.** The skew-symmetry residual $\|[D,P_\tau] + [D,P_\tau]^\top\|_F / \|[D,P_\tau]\|_F$ is exactly $0.00 \times 10^{0}$ at all four primorial levels.

### 2.4 The Hermitian Operator $H(\alpha)$

**Definition 2.2.** The *active transport operator* at shear amplitude $\alpha \geq 0$ is

$$H(\alpha) = \frac{D_{\text{sym}}}{\lambda_P} + i\alpha \, \frac{[D_{\text{sym}},\, P_\tau]}{\lambda_P}$$

where $\lambda_P = \lambda_{\max}(D_{\text{sym}})$ is the Perron eigenvalue.

**Proposition 2.3** (Hermiticity). $H(\alpha)$ is Hermitian for all $\alpha \in \mathbb{R}$, and therefore has purely real eigenvalues.

*Proof.* Write $H = A + i\alpha B$ where $A = D_{\text{sym}}/\lambda_P$ (real symmetric) and $B = [D_{\text{sym}}, P_\tau]/\lambda_P$ (real skew-symmetric by Proposition 2.1). Then $H^\dagger = A^\top - i\alpha B^\top = A + i\alpha B = H$. $\square$

**Remark 2.4** (Physical interpretation). The operator $H(\alpha)$ models the interaction between additive spatial diffusion ($D_{\text{sym}}$, the passive cage geometry) and multiplicative helical resonances ($[D, P_\tau]$, the active shear). The parameter $\alpha$ sweeps from pure diffusion ($\alpha = 0$) to shear-dominated dynamics ($\alpha \gg 1$). The normalization by $\lambda_P$ ensures $H(0)$ has Perron eigenvalue exactly 1.

### 2.5 Derivation from Spectral Isotropy

The *Spectral Isotropy* identity [1, Section 9.7] $D + D^\top = m(J + I)$ for the directed distance matrix $D$ allows $H(\alpha)$ to be rewritten entirely in terms of the antisymmetric part of $D$. Define $A = (D - D^\top)/2$ (integer-valued for even $m$, half-integer-valued for odd $m$).

**Identity 2.5** (Symmetric folding). 
$$D_{\text{sym}} = \tfrac{m}{2}(J - I) - |A|$$
where $|A|$ denotes the Hadamard (element-wise) absolute value, $|A|_{ij} = |A_{ij}|$ (written $|A|_{\text{ew}}$ when emphasis on the element-wise nature is needed), and $J$ is the all-ones matrix. (This is *not* the matrix absolute value $(A^\top A)^{1/2}$; the element-wise folding is essential — it preserves additive geometry while scattering the cotangent eigenphases into the flat band.)

*Derivation.* From $D + D^\top = m(J + I)$ and $D_{\text{sym}} = (D + D^\top)/2 - \text{diag}$, one obtains $D_{\text{sym}} = m/2 \cdot (J + I) - D = m/2 \cdot (J - I) - (D - m/2 \cdot I)$. The palindromic symmetry $d(r, m-r) = m/2$ forces $D_{\text{sym}}$ to carry exactly the spectral content of $|A|$, subtracted from the constant background $m/2 \cdot (J - I)$.

| $m$ | $\varphi(m)$ | $\|D_{\text{sym}} - [m/2 \cdot (J-I) - |A|]\|_F$ | $A$ entry type |
|---|---|---|---|
| 30 | 8 | $0.00 \times 10^{0}$ | integer |
| 105 | 48 | $0.00 \times 10^{0}$ | half-integer |
| 210 | 48 | $0.00 \times 10^{0}$ | integer |
| 1155 | 480 | $0.00 \times 10^{0}$ | half-integer |
| 2310 | 480 | $0.00 \times 10^{0}$ | integer |

**Identity 2.6** (Commutator reduction). Since the constant background $m/2 \cdot (J - I)$ commutes with every matrix (its column sums are constant),
$$[D_{\text{sym}}, P_\tau] = -[|A|, P_\tau].$$

**Verification.** $\|[D_{\text{sym}}, P_\tau] + [|A|, P_\tau]\|_F = 0.00 \times 10^{0}$ at all five moduli.

**Corollary 2.7** (Canonical decomposition of $H(\alpha)$). Substituting Identities 2.5 and 2.6 into Definition 2.2:
$$H(\alpha) = \frac{m/2 \cdot (J - I) - |A|}{\lambda_P} - i\alpha\,\frac{[|A|, P_\tau]}{\lambda_P}.$$

Both the real and imaginary parts of $H(\alpha)$ are determined by a single object: the element-wise absolute value $|A|$ of the antisymmetric part of the directed distance matrix. The constant background $m/2 \cdot (J - I)$ contributes only a uniform shift to the real part.

**Remark 2.8** (Connection to cotangent eigenphase formula). The Spectral Isotropy paper [3] establishes that $\operatorname{Im}(\lambda_k) = -R \cdot \cot(k\pi/\varphi(m))$ for the eigenvalues of $D$. Identity 2.5 connects this to $H(\alpha)$: the cotangent eigenphases of $D$ become, through the folding $A \to |A|$, the spectral content that drives localization in the active transport operator.

### 2.6 The Unexplained Isospectrality of $|A|/m$

The element-wise folding $A \to |A|$ (where $|A|_{ij} = |A_{ij}|$) is a non-linear operation on matrix entries. Standard matrix analysis provides no general mechanism for such operations to preserve spectra across matrices with different spatial wirings. Yet computation shows that $|A|/m$ has *identical* eigenvalues for matched pairs (e.g., $m = 1155$ and $m = 2310$, both $\varphi = 480$) to machine precision ($< 10^{-13}$), despite fundamentally different coprime gap structures (odd gaps allowed at odd $m$, only even gaps at even $m$). The eigenvectors of $|A|/m$ are NOT related by conjugation — the Frobenius non-conjugacy error grows with $\varphi$ — so this is not unitary equivalence. The isospectrality may be a consequence of a hidden symmetry in the singular value structure of $A$, but no algebraic proof exists. This is listed as Problem 15.

---

## 3. The $1/\sqrt{2}$ Commutator Ratio

### 3.1 Computational Evidence

The Frobenius norm ratio $\|[D,P_\tau]\|_F / \|D\|_F$ was computed at four primorial levels:

| $m$ | $\varphi(m)$ | $\|[D,P_\tau]\|_F / \|D\|_F$ | $(\text{ratio})^2$ | Error from $1/2$ |
|---|---|---|---|---|
| 30 | 8 | 0.528498 | 0.279310 | $-0.2207$ |
| 210 | 48 | 0.697550 | 0.486576 | $-0.0134$ |
| 2310 | 480 | 0.706280 | 0.498831 | $-0.0012$ |
| 30030 | 5760 | 0.706982 | 0.499823 | $-0.00018$ |

The ratio converges to $1/\sqrt{2} \approx 0.707107$, equivalently the squared ratio converges to $1/2$, with the error dropping by approximately an order of magnitude at each step.

### 3.2 Derivation from Known Results

The $1/\sqrt{2}$ limit is not a new constant but a derived consequence of two results from the companion paper [1]:

**Result A** (Irreconcilability coupling constant, [1] Theorem 9.7):

$$\frac{\|[D, P_\tau]\|_F}{\lambda_P} \longrightarrow \sqrt{\frac{2}{3}}$$

**Result B** (Frobenius-to-Perron ratio, [1] Section 9.7 Step 4):

$$\frac{\|D\|_F}{\lambda_P} \longrightarrow \frac{2}{\sqrt{3}}$$

Dividing Result A by Result B:

$$\frac{\|[D, P_\tau]\|_F}{\|D\|_F} = \frac{\|[D, P_\tau]\|_F / \lambda_P}{\|D\|_F / \lambda_P} \longrightarrow \frac{\sqrt{2/3}}{2/\sqrt{3}} = \frac{\sqrt{2}}{\sqrt{3}} \cdot \frac{\sqrt{3}}{2} = \frac{\sqrt{2}}{2} = \frac{1}{\sqrt{2}}$$

Squaring:

$$\frac{\|[D, P_\tau]\|_F^2}{\|D\|_F^2} \longrightarrow \frac{1}{2}$$

### 3.3 The $3/4$ Correlation

Equivalently, the $1/\sqrt{2}$ result restates the alignment theorem from [1]. Using the identity $\|[D,P_\tau]\|_F^2 = 2\|D\|_F^2 - 2\,\text{tr}(P_\tau^\top D P_\tau D)$ (see [1] Section 9.7 Step 1), the squared ratio becomes:

$$\frac{\|[D,P_\tau]\|_F^2}{\|D\|_F^2} = 2 - 2\,\frac{\text{tr}(P_\tau^\top D P_\tau D)}{\|D\|_F^2} = 2(1 - A)$$

where $A = \text{tr}(P_\tau^\top D P_\tau D) / \|D\|_F^2$ is the *alignment* of $D$ under conjugation by $P_\tau$.

The squared ratio converging to $1/2$ is equivalent to $A \to 3/4$: the palindromic distance matrix "remembers" exactly three-quarters of its Hilbert–Schmidt structure when conjugated by the multiplicative inversion. The remaining one-quarter — the *irreducible non-commutativity* — is what the commutator captures.

The present computation provides an independent verification of $A \to 3/4$ through a completely different code path (constructing $H(\alpha)$ and measuring its spectral properties) than the original proof in [1] (which computed the alignment directly via entry-wise correlation analysis).

---

## 4. Flat-Band Condensation

### 4.1 Eigenvalue Variance Scaling

At any fixed $\alpha > 0$, the variance of the eigenvalues of $H(\alpha)$ decays as $\mathcal{O}(1/\varphi(m))$. We demonstrate this at $\alpha = \sqrt{2/3}$:

| $m$ | $\varphi(m)$ | $\text{Var}(\lambda)$ | $\text{Var} \cdot \varphi$ | Relative change |
|---|---|---|---|---|
| 30 | 8 | 0.203022 | 1.624 | — |
| 210 | 48 | 0.036839 | 1.768 | $+8.9\%$ |
| 2310 | 480 | 0.003702 | 1.777 | $+0.5\%$ |
| 30030 | 5760 | 0.000309 | 1.780 | $+0.2\%$ |

The product $\text{Var}(\lambda) \cdot \varphi(m)$ converges to a constant $V_\infty \approx 1.78$. The scaling is:

$$\text{Var}(\lambda) \sim \frac{V_\infty}{\varphi(m)}, \qquad V_\infty \approx 1.780 \pm 0.002.$$

**Remark 4.1.** The flat-band condensation is structurally guaranteed by the normalization $H(0)$'s Perron eigenvalue being 1. As $\varphi(m) \to \infty$, the matrix $D/\lambda_P$ concentrates its spectral mass: the Perron eigenvector (which approaches uniform by equidistribution of coprime residues [1, Section 9.7]) captures the leading eigenvalue, and all remaining $\varphi(m) - 1$ eigenvalues are $o(1)$. The additive shear term $i\alpha [D, P_\tau]/\lambda_P$ perturbs these into a band of width $\mathcal{O}(1/\sqrt{\varphi})$, producing $\text{Var} = \mathcal{O}(1/\varphi)$. The condensation *per se* is therefore a consequence of Perron–Frobenius concentration for dense non-negative matrices [5], not a topological phenomenon. What is arithmetically specific is (a) the convergence to a definite constant $V_\infty \approx 1.780$, which a random dense matrix would not produce; (b) the *structure* that survives in the flat band — the CRT block decomposition (Section 10.3) with exact scalar operators, which is absent from generic concentrated spectra; and (c) the sharp phase transition at $\alpha_c$, which requires the interplay of additive palindromic symmetry and multiplicative inversion.

### 4.2 The Three-Point Spectral Measure

Examination of the sorted eigenvalue spectra at $\alpha = \sqrt{2/3}$ reveals the limiting spectral structure:

| $m$ | $\lambda_{\max}$ (Perron) | $\lambda_{\min}$ (anti-Perron) | Fraction of $\lambda_i$ in $[-0.05, 0.05]$ |
|---|---|---|---|
| 30 | $+1.0000$ | $-0.6341$ | 37.5% (3 of 8) |
| 210 | $+1.0000$ | $-0.6032$ | 87.5% (42 of 48) |
| 2310 | $+1.0000$ | $-0.5917$ | 98.3% (472 of 480) |
| 30030 | $+1.0000$ | $-0.5908$ | 99.9% (5752 of 5760) |

The spectrum progressively collapses onto three isolated points: $\{+1, \approx -0.59, 0\}$. The zero mode is massively degenerate from $\varphi(m) = 480$ onward. The anti-Perron eigenvalue converges to $\lambda_- \approx -0.591 \pm 0.001$.

**Remark 4.2.** The convergence of $\lambda_-$ toward $\approx -0.591$ raises the question of whether this is a clean algebraic number. The value $-1 + \sqrt{2/3} \approx -0.1835$ does not match; $-\sqrt{2/3} \approx -0.8165$ does not match; $1/\sqrt{2} - \sqrt{2} \approx -0.707$ does not match. We leave the identification of $\lambda_-$ as an open problem.

### 4.3 Inverse Participation Ratio

The inverse participation ratio (IPR) [6] measures eigenvector localization. For an eigenvector $\mathbf{v}$ with $\|\mathbf{v}\| = 1$, the IPR is $\sum_i |v_i|^4$. A perfectly delocalized state (all components equal) has $\text{IPR} = 1/n$; a fully localized state has $\text{IPR} = 1$.

| $m$ | Type | $\varphi(m)$ | Mean IPR | Delocalized ref ($1/\varphi$) | IPR ratio |
|---|---|---|---|---|---|
| 30 | primorial | 8 | 0.178971 | 0.125000 | $1.43\times$ |
| 210 | primorial | 48 | 0.039937 | 0.020833 | $1.92\times$ |
| **385** | **non-primorial** | **240** | **0.008765** | **0.004167** | $\mathbf{2.10\times}$ |
| 2310 | primorial | 480 | 0.005115 | 0.002083 | $2.46\times$ |
| **1155** | **non-primorial** | **480** | **0.005139** | **0.002083** | $\mathbf{2.47\times}$ |
| **15015** | **non-primorial** | **5760** | **0.000464** | **0.000174** | $\mathbf{2.67\times}$ |
| 30030 | primorial | 5760 | 0.000501 | 0.000174 | $2.88\times$ |

The IPR ratio increases monotonically with system size across both primorial and non-primorial moduli. Two matched-dimension pairs isolate the torsion effect:

| $\varphi$ | Primorial | IPR ratio | Non-primorial | IPR ratio | Torsion effect |
|---|---|---|---|---|---|
| 480 | $m = 2310$ | $2.46\times$ | $m = 1155$ | $2.47\times$ | $-0.4\%$ (negligible) |
| 5760 | $m = 30030$ | $2.88\times$ | $m = 15015$ | $2.67\times$ | $+7.9\%$ |

At $\varphi = 480$, torsion has no measurable effect on eigenvector localization. At $\varphi = 5760$, the primorial exceeds the non-primorial by $7.9\%$. The torsion–localization coupling is **dimension-dependent** and strengthens with $\varphi$.

**Observation 4.3.** The IPR ratio sequence is consistent with $\text{IPR ratio} \approx 0.38 \ln(\varphi)$, though seven data points cannot distinguish logarithmic growth from alternative scaling forms (e.g., a small power law $\varphi^\beta$ with $\beta \ll 1$). The non-primorial controls ($m = 385, 1155, 15015$) fall on or slightly below the primorial trend line; the torsion–localization gap widens from $< 1\%$ at $\varphi = 480$ to $\sim 8\%$ at $\varphi = 5760$.

**Observation 4.4** (Perron eigenvector uniformity). Direct extraction of the Perron eigenvector at $\alpha_c(30030)$ confirms that it is perfectly flat: Fourier analysis shows $100\%$ DC power and $0.0\%$ spectral power at every tested period ($3, 5, 7, 11, 13, 17, 24, 48, 480, 5760$). The mod-4 chiral asymmetry is $\mathcal{A} = -2.5 \times 10^{-6}$ with a $50.000\% / 50.000\%$ mass split, constant to $\sigma = 0$ across the crash window. This uniformity is a theorem-level consequence of flat-band condensation: as $\varphi \to \infty$, the normalized matrix $D_{\text{sym}}/\lambda_P$ is a rank-1 perturbation of the flat band, forcing the Perron eigenvector toward $\mathbf{1}/\sqrt{\varphi}$. Any hypothesis requiring spatial structure in the dominant eigenvector (periodic modulation, chiral asymmetry, coil winding patterns) is ruled out by the flat-band mechanism itself.

### 4.5 Torsion Crossover: Isospectral Non-Conjugacy

The torsion–localization gap between matched pairs (same $\varphi$, different 2-torsion) reverses sign as $\varphi$ increases. Three experiments identify the mechanism.

*Experiment 1: Counterfactual $P_\tau$ swap.* For the matched pair $m = 1155$ (odd) and $m = 2310$ (even), both with $\varphi = 480$, we construct four hybrid Hamiltonians by pairing each $D_{\text{sym}}$ with each $P_\tau$:

| Hybrid | $D_{\text{sym}}$ from | $P_\tau$ from | Mean IPR ratio | Max eigenvalue |
|---|---|---|---|---|
| Odd $D$ + Odd $P$ | $m = 1155$ | $m = 1155$ | $2.4669\times$ | natural |
| Odd $D$ + Even $P$ | $m = 1155$ | $m = 2310$ | $2.4669\times$ | identical |
| Even $D$ + Odd $P$ | $m = 2310$ | $m = 1155$ | $2.2802\times$ | natural |
| Even $D$ + Even $P$ | $m = 2310$ | $m = 2310$ | $2.2802\times$ | identical |

Swapping $P_\tau$ changes neither the eigenvalues (to machine precision) nor the IPR. $D_{\text{sym}}$ alone determines the localization.

*Experiment 2: Isospectral non-conjugacy of $|A|/m$.* The matrix $|A|/m$ has **identical eigenvalues** for matched pairs (max difference $< 10^{-13}$), but it is **not conjugate** under the coprime bijection $\sigma: \mathcal{R}_{\text{odd}} \to \mathcal{R}_{\text{even}}$: the Frobenius error $\||A|_{\text{odd}}/m - \sigma |A|_{\text{even}}/m\, \sigma^\top\|_F$ is $1.5$ at $\varphi = 8$ and $9.8$ at $\varphi = 48$ (the two smallest matched-pair sizes tested; the $\varphi = 480$ value remains to be computed). In contrast, $P_\tau$ is **exactly conjugate** under $\sigma$ (Frobenius error $= 0.0$ at all sizes).

Same frequencies, different spatial wiring — the eigenvectors of $|A|/m$ (and therefore $D_{\text{sym}}$, by Identity 2.5) are shaped differently at odd and even $m$.

*Experiment 3: Kurtosis explosion.* The eigenvector weight distributions diverge dramatically:

| $m$ | $\varphi$ | Max IPR / $n$ | Kurtosis | Skewness |
|---|---|---|---|---|
| 1155 (odd) | 480 | $21.85\times$ | 117.0 | 8.25 |
| 2310 (even) | 480 | $6.26\times$ | 6.9 | 1.34 |

At odd $m$, the coprime gaps include odd values $\{1, 2, 3, 4, 5, \ldots\}$, allowing isolated coprime clusters that host topological hotspots ($21.85\times$ delocalized reference). At even $m$, coprime gaps are restricted to even values $\{2, 4, 6, 8, \ldots\}$, spacing the residues more uniformly and capping the hotspot intensity ($6.26\times$).

The mean IPR (averaged over all eigenvectors) nonetheless converges because flat-band condensation — which carries $> 98\%$ of eigenvectors to near-delocalized — overwhelms the tail. The **sign** of the torsion gap (odd more localized vs. even more localized) depends on which effect dominates the mean: the coprime-gap hotspots in $D_{\text{sym}}$ or the flat-band dilution. At $\varphi \lesssim 480$, the hotspot contribution is small and the sign fluctuates with factorization; at $\varphi \gtrsim 5000$, it stabilizes at $+7$–$8\%$ in favor of the even (primorial) modulus.

**Required condition for crossover existence.** The flat-band eigenvalue variance decays as $\text{Var}(\lambda) \sim V_\infty / \varphi$ (Section 4.1), which drives the mean IPR toward the delocalized limit as the flat-band subspace expands. Simultaneously, the coprime-gap kurtosis $\kappa(\varphi)$ grows at odd $m$ as topological hotspots (Section 4.5, Experiment 3) concentrate eigenvector weight. The crossover demands that the dilution rate of the mean IPR — driven by the expanding flat-band dimensionality — strictly outpaces the growth rate of the spatial kurtosis driven by odd-gap hotspots. If $\kappa$ grew faster than the flat-band dilution, the odd-modulus kurtosis advantage would dominate at all scales, and no crossover would occur. The observed crossover between $\varphi = 480$ and $\varphi = 1440$ implies that the flat-band condensation rate is the dominant asymptotic effect. Proving this — or determining the exact scaling of $\kappa(\varphi)$ — is Problem 16.

**Observation 4.6 (Klein Four-Group Orbit Structure and the $O(2)$ Leaky Trap).** The coprime residues $r \in (\mathbb{Z}/m\mathbb{Z})^\times$ partition into orbits under the Klein four-group $\langle \sigma, \tau \rangle \cong \mathbb{Z}_2 \times \mathbb{Z}_2$, where $\sigma(r) = m - r$ (additive palindrome involution) and $\tau(r) = r^{-1} \bmod m$ (multiplicative inverse involution). Each orbit has the form $O(r) = \{r,\; m{-}r,\; r^{-1},\; (m{-}r)^{-1}\}$ (size 4 generically, size 2 when $r = r^{-1}$ or $r = m{-}r$).

At $m = 1155$ (odd primorial), $\varphi = 480$ residues decompose into 124 Klein-4 orbits (116 of size 4, 8 of size 2). The orbit $O(2) = \{2, 1153, 578, 577\}$ exists and acts as a *leaky trap* for flat-band eigenvectors: the most localized mode ($\text{IPR} \times \varphi = 21.8$) concentrates $41.3\%$ of its mass inside $O(2)$, spanning effectively $\sim 5.5$ orbits total (Orbit Participation Ratio, OPR $= 5.53$ out of 124; analogous to IPR but in orbit space). Four of the top-10 localized modes have $O(2)$ as their dominant orbit, with $O(2)$'s mass fraction $8\times$ the background mean ($5.26\%$).

The trapping mechanism is a boundary effect: residue $r = 2$ is the smallest non-trivial coprime, sitting at the edge of the lattice where the coprime gap to $r = 1$ is exactly 1 (an odd gap). The additive component $D_{\text{sym}}$ creates a steep potential well at this boundary, while the multiplicative commutator $-i\alpha[|A|, P_\tau]$ pumps amplitude through the orbit network — hence the leakage into $\sim 16$ additional orbits.

At the even primorial $m = 2310$, residue 2 is deleted by the coprime sieve ($\gcd(2, 2310) = 2 \neq 1$), erasing $O(2)$ entirely. Without this boundary trap, orbit-space localization collapses: the most localized mode has OPR $= 26.1$ (vs $5.5$), peak orbit mass drops to $10.9\%$ (vs $41.3\%$), and eigenvector kurtosis falls from 114 to 3.9.

The $\tau$-parity expectation values $\langle v | P_\tau | v \rangle$ for the top-10 flat-band modes at odd $m$ range from $0.005$ to $0.49$ — confirming that eigenvectors of $H(\alpha)$ do not diagonalize $P_\tau$. The commutator $[|A|, P_\tau]$ generates rotations between $\tau$-eigenspaces rather than preserving them. At even $m$, without the $O(2)$ funnel disrupting $\tau$-symmetry, modes can settle into near-$\tau$-eigenstates (values up to $0.82$).

The flat band is therefore not a collection of decoupled Klein-4 rings, but a single interconnected orbit network where additive boundary traps create semi-localized attractor states. The kurtosis explosion (Section 4, Experiment 3) is the $O(2)$ funnel effect: the deepest additive trap coincides with the tightest Klein-4 orbit at the lattice boundary, and its deletion by the coprime sieve at even $m$ removes the dominant localization mechanism.

---

## 5. Perron Gap Universality

The *Perron gap* $\Delta_P = \lambda_1 - \lambda_2$ is the difference between the largest eigenvalue ($+1$ by normalization) and the second-largest eigenvalue of $H(\alpha)$.

At $\alpha = \sqrt{2/3}$:

| $m$ | $\varphi(m)$ | Perron gap $\Delta_P$ | Spectral width $W$ | $\Delta_P / W$ |
|---|---|---|---|---|
| 30 | 8 | 0.856826 | 1.634124 | 0.52433 |
| 210 | 48 | 0.813760 | 1.603232 | 0.50757 |
| 2310 | 480 | 0.814621 | 1.591708 | 0.51179 |
| 30030 | 5760 | 0.814642 | 1.590752 | 0.51211 |

Three features emerge:

1. **Perron gap convergence.** $\Delta_P \to 0.8146$ with oscillatory approach: the sequence overshoots at $m = 30$, undershoots at $m = 210$, then converges from below. The gap between $m = 2310$ and $m = 30030$ is only $0.000021$, indicating rapid convergence.

2. **Spectral width convergence.** $W \to 1.591$ with error $< 0.001$ by $m = 2310$.

3. **Gap ratio convergence.** $\Delta_P / W \to 0.512$, tantalizingly close to $1/2$ but not exactly $1/2$ (the error at $m = 30030$ is $+0.024$).

**Remark 5.1.** The peak gap ratio (supremum over $\alpha$) occurs at $\alpha = 0$ for all four moduli, where it takes the value $\approx 0.711$, strikingly close to $1/\sqrt{2} \approx 0.7071$. The sequence $0.688 \to 0.707 \to 0.711 \to 0.712$ suggests convergence to a value *slightly above* $1/\sqrt{2}$. This may be a finite-size artifact.

---

## 6. The Gap Position Phase Transition

### 6.1 The Transition

The *gap position* is defined as $g = j_{\max} / (n - 2)$ where $j_{\max}$ is the index of the largest consecutive gap $\lambda_{j+1} - \lambda_j$ in the sorted eigenvalue sequence. The normalization maps $g \in [0, 1]$ with $g = 0$ meaning the gap is between the two smallest eigenvalues (bulk/bottom), $g = 0.5$ meaning the gap is at the spectral center, and $g = 1$ meaning the gap is between the two largest (Perron gap).

For $\alpha < 1.2$, all four moduli have $g = 1.0$ — the Perron gap is always the largest gap in the spectrum. This means the additive diffusion ($D_{\text{sym}}$) controls the spectral structure at moderate shear.

At a critical shear amplitude $\alpha_c$, the gap position crashes catastrophically from $1.0$ to $\approx 0.0$:

| $m$ | $I$ | Approx. $\alpha_c$ (onset of crash) |
|---|---|---|
| 30 | 0 | $\sim 1.3$ |
| 210 | 1 | $\sim 1.2$ |
| 2310 | 2 | $\sim 1.2$ |
| 30030 | 1 | $\sim 1.2$ |

The transition sharpens with increasing $\varphi(m)$. At $m = 30$ (only 8 residues), the crash is gradual. At $m = 30030$ (5760 residues), it is essentially a step function.

### 6.2 Fine-Resolution Sweep

To pinpoint $\alpha_c$, a high-resolution sweep was conducted with 222 $\alpha$-values: 200 fine points in $[1.0, 1.45]$ (spacing $\Delta\alpha \approx 0.00226$), 20 context points in $[0.5, 1.0]$, and mandatory probes at $\sqrt{2/3}$ and $\sqrt{3/2}$. The crash point was defined as the interpolated $\alpha$ where the gap position crosses $0.5$.

| $m$ | $I$ | $\varphi$ | Crash $\alpha$ (coarse) | Error from $\sqrt{3/2}$ | Steepest derivative |
|---|---|---|---|---|---|
| 30 | 0 | 8 | N/A (no crash) | — | 0.00 |
| 210 | 1 | 48 | 1.27025 | $+3.72\%$ | $-433$ |
| 2310 | 2 | 480 | 1.23857 | $+1.13\%$ | $-441$ |
| 30030 | 1 | 5760 | 1.23857 | $+1.13\%$ | $-442$ |

The transition is a first-order spectral reorganization: the derivative $d(\text{gap\_pos})/d\alpha$ peaks at $\sim -442$, essentially a delta function.

The apparent plateau between $m = 2310$ and $m = 30030$ (identical crash values to five decimal places) initially appeared to falsify $\sqrt{3/2}$. This was a resolution artifact, as explained in Section 6.3.

### 6.3 Resolution Artifact and the Micro-Sweep

The coarse grid (spacing $\Delta\alpha \approx 0.00226$) created a critical methodological artifact: the **transition width** at $m = 2310$ was found to be approximately $3.8 \times 10^{-6}$ — over 500 times narrower than the grid spacing. Both $m = 2310$ and $m = 30030$ crashed between the *same two adjacent grid points* ($\alpha \approx 1.23744$ and $\alpha \approx 1.23970$), producing identical interpolated crash values despite their true crash points differing.

A micro-sweep at $m = 2310$ with 2003 $\alpha$-values in a window of width $0.015$ (spacing $\Delta\alpha \approx 7.5 \times 10^{-6}$) resolved the true transition:

**Table 6.2.** Micro-sweep results for $m = 2310$ ($\varphi = 480$).

| Quantity | Coarse grid | Micro-sweep |
|---|---|---|
| Crash $\alpha$ | 1.23857 | 1.2389219 |
| Transition width (75% → 25%) | unresolved | $3.76 \times 10^{-6}$ |
| Steepest $d/d\alpha$ | $-442$ | $-132{,}988$ |
| Error from $\sqrt{3/2}$ | $+1.13\%$ | $+1.16\%$ |
| Error from $\sqrt{135/88}$ | — | $+0.027\%$ |

The true transition derivative is 300 times steeper than reported by the coarse grid. The transition is extraordinarily sharp — the gap position drops from 0.75 to 0.25 in an interval of less than $4 \times 10^{-6}$.

### 6.4 The $\sqrt{3/2}$ Conjecture: Falsified

The coarse-grid "plateau" was an artifact, but the $\sqrt{3/2}$ conjecture is nonetheless falsified. The micro-sweep gives a precise crash point of $\alpha_c(2310) = 1.2389219$, which is $+1.16\%$ above $\sqrt{3/2} \approx 1.22474$. The product $\kappa \cdot \alpha_c \approx \sqrt{2/3} \times 1.2389 \approx 1.0116$, not $1.0$.

### 6.5 The $\sqrt{135/88}$ Conjecture

The micro-sweep crash value $1.2389219$ is remarkably close to $\sqrt{135/88} \approx 1.2385842$:

$$\sqrt{\frac{135}{88}} = \sqrt{\frac{3^3 \cdot 5}{2^3 \cdot 11}} \approx 1.238584235767$$

Error from $\sqrt{135/88}$: $+0.027\%$ at $m = 2310$ ($\varphi = 480$).

The factorization is striking: the numerator $135 = 3^3 \cdot 5$ and denominator $88 = 2^3 \cdot 11$ are built entirely from the first five primes $\{2, 3, 5, 11\}$, with 7 absent — the same primes that divide the primorial $m = 2310 = 2 \cdot 3 \cdot 5 \cdot 7 \cdot 11$.

**Error scaling across moduli.** A non-primorial control ($m = 385 = 5 \cdot 7 \cdot 11$, $\varphi = 240$), which also exhibits the phase transition, permits a convergence table:

| $m$ | Type | $\varphi$ | Crash $\alpha$ | Error from $\sqrt{135/88}$ | Width | Steepest $dg/d\alpha$ |
|---|---|---|---|---|---|---|
| 210 | primorial | 48 | 1.27025 | $+2.56\%$ | — | — |
| 385 | non-primorial | 240 | 1.24227 | $+0.30\%$ | $6.0 \times 10^{-6}$ | $-82{,}817$ |
| 2310 | primorial | 480 | 1.23892 | $+0.027\%$ | $3.8 \times 10^{-6}$ | $-133{,}000$ |
| 30030 | primorial | 5760 | 1.23877 | $+0.015\%$ | $1.0 \times 10^{-7}$ | $-4{,}989{,}133$ |

The error decreases monotonically with $\varphi$, but the convergence is sharply decelerating — inconsistent with any stable power law. The local scaling exponent $\beta$ (defined by $\text{error} \propto 1/\varphi^\beta$ between successive primorial pairs) drops dramatically: $\beta \approx 2.0$ from $\varphi = 48 \to 480$ (error $2.56\% \to 0.027\%$, a $95\times$ reduction for $10\times$ growth in $\varphi$), but only $\beta \approx 0.24$ from $\varphi = 480 \to 5760$ (error $0.027\% \to 0.015\%$, a $1.8\times$ reduction for $12\times$ growth). This deceleration is consistent with logarithmic convergence or a finite-size crossover, not an algebraic scaling law. The transition width simultaneously collapsed $38\times$ (from $3.8 \times 10^{-6}$ to $1.0 \times 10^{-7}$) and the steepest derivative grew $37\times$ (from $-133{,}000$ to $-4{,}989{,}133$), indicating a sharpening first-order transition.

**Conjecture 6.1** (Critical shear). The spectral phase transition occurs at

$$\alpha_c = \sqrt{\frac{135}{88}} = \sqrt{\frac{3^3 \cdot 5}{2^3 \cdot 11}}$$

in the limit $\varphi(m) \to \infty$, with finite-size corrections scaling sub-quadratically.

**Caveat.** This is a numerical conjecture, not a proven identity: $\sqrt{135/88}$ is the best-fit algebraic approximant to the observed crash values at $m \leq 30030$. Four data points cannot distinguish an algebraic limit from a nearby transcendental with logarithmic corrections. The convergence is sharply decelerating: the local exponent drops from $\beta \approx 2.0$ ($\varphi = 48 \to 480$) to $\beta \approx 0.24$ ($\varphi = 480 \to 5760$), which is inconsistent with any stable power law and may indicate logarithmic convergence. The conjecture rests on the suggestive factorization $3^3 \cdot 5 / (2^3 \cdot 11)$ and monotone convergence from above, not on a derivation from the group structure. **Falsification criterion:** The $\sqrt{3/2}$ candidate was ruled out because the observed crash values converge *away* from it (error $+1.16\%$ at $m = 2310$, stable under refinement). In contrast, $\sqrt{135/88}$ shows monotone convergence *toward* itself ($+2.56\% \to +0.30\% \to +0.027\% \to +0.015\%$). A candidate is falsified when the error stabilizes or grows at higher $\varphi$; $\sqrt{135/88}$ is retained because the error decreases at every step, though the deceleration leaves open the possibility of a nearby transcendental limit. We cannot currently explain why these specific prime powers appear, nor why 7 — a factor of the primorial — is absent from the ratio. Resolution requires either an analytic derivation or computation at $m = 510510$ ($\varphi = 92160$), which would test to $\sim 0.005\%$ precision.

The $m = 30030$ micro-sweep (503 points, spacing $\sim 1.4 \times 10^{-8}$) resolves the crash at $\alpha_c(30030) = 1.2387666$, with error $+0.015\%$ from $\sqrt{135/88}$, transition width $1.0 \times 10^{-7}$, and steepest derivative $-4{,}989{,}133$.

Convergence from $m = 2310$ to $m = 30030$ is **confirmed but sharply decelerating**: error dropped from $+0.027\%$ to $+0.015\%$ (a $1.8\times$ reduction) while $\varphi$ grew from 480 to 5760 (a $12\times$ increase), giving a local exponent $\beta \approx 0.24$ — far below the $\beta \approx 2.0$ measured at the previous step. The deceleration rules out a stable algebraic scaling law. Correspondence to $\sqrt{135/88}$ now stands at five significant figures ($1.23877$ vs $1.23858$), but whether the limit is exactly algebraic or involves logarithmic corrections remains open.

The transition dynamics, however, are unambiguously first-order: the width collapsed $38\times$ and the derivative grew $37\times$ between $m = 2310$ and $m = 30030$. The related products:

$$\kappa \cdot \alpha_c = \sqrt{\frac{2}{3}} \cdot \sqrt{\frac{135}{88}} = \sqrt{\frac{270}{264}} = \sqrt{\frac{45}{44}} \approx 1.01130$$

This is close to but not equal to 1.

### 6.6 Non-Primorial Control: $m = 385$

The modulus $m = 385 = 5 \cdot 7 \cdot 11$ (non-primorial, $\varphi = 240$) provides a critical control test: *does the phase transition require primorial structure, or only a sufficient number of coprime residues?*

A micro-sweep at $m = 385$ reveals a clear phase transition at $\alpha_c(385) = 1.24227$, with transition width $6.0 \times 10^{-6}$ and steepest derivative $-82{,}817$. The error from $\sqrt{135/88}$ is $+0.30\%$.

**Conclusion.** The spectral phase transition is not specific to primorials. It occurs for any squarefree modulus with sufficient Euler totient, converging toward the same critical shear. This rules out explanations that depend on primorial-specific properties (such as the consecutive-prime factorization) and suggests the transition is a universal feature of the active transport operator on coprime residue lattices.

### 6.7 Post-Transition Structure

Beyond the crash, the gap position remains at $g \approx 0$ out to $\alpha = 2.5$ with no second transition. The Perron gap decays monotonically from $\Delta_P \approx 0.66$ at $\alpha = 1.24$ to $\Delta_P \approx 0.17$ at $\alpha = 2.5$, but the bulk gap always dominates. The transition at $\alpha_c$ is the unique spectral reorganization point.

### 6.8 IPR Across the Transition

The inverse participation ratio was computed at $\alpha = \sqrt{3/2}$ (i.e., near but not exactly at the true critical value) for all four primorials:

| $m$ | $\varphi(m)$ | IPR ratio at $\alpha = \sqrt{3/2}$ |
|---|---|---|
| 30 | 8 | $1.45\times$ |
| 210 | 48 | $1.88\times$ |
| 2310 | 480 | $2.28\times$ |
| 30030 | 5760 | $2.63\times$ |

Localization grows monotonically with system size across the transition, with no spike at $\alpha_c$. The transition rearranges spectral *gaps*, not eigenvector *structure*.

### 6.9 Physical Interpretation

Below $\alpha_c$: the additive distance structure (the "cage") controls the spectral gap. Perturbations from the multiplicative inversion permutation are subdominant. The Perron mode is isolated.

Above $\alpha_c$: the multiplicative commutator term overwhelms the additive structure. Internal bulk gaps, generated by the complex-valued off-diagonal entries of $i\alpha[D,P_\tau]/\lambda_P$, grow larger than the Perron gap. The spectral structure "switches authority" from the additive to the multiplicative sector.

---

## 7. Falsified Predictions

### 7.1 The $\sqrt{2/3}$ Quench Hypothesis

**Prediction:** A spectral phase transition (block-diagonalization, flat-band collapse, or sharp gap feature) occurs at $\alpha = \sqrt{2/3} \approx 0.8165$.

**Observation:** All four spectral indicators (block separation, gap ratio, eigenvalue variance relative to $\alpha = 0$, gap position) pass smoothly through $\alpha = 0.8165$ without discontinuity or inflection. The value $\alpha = \sqrt{2/3}$ is unremarkable in the spectral sweep data.

**Why it fails:** The irreconcilability coupling constant $\sqrt{2/3}$ describes the *asymptotic norm* of the commutator relative to $\lambda_P$ — a static property of the commutator's Frobenius norm, not a dynamical threshold. There is no reason for the eigenvalue structure of $H(\alpha)$ to exhibit a singularity at the shear value where the commutator norm equals one particular normalization of $D$. The relevant dynamical threshold is instead $\alpha_c \approx \sqrt{135/88}$ (Section 6.5), which marks the point where the commutator's contribution to $H(\alpha)$ becomes spectrally dominant over the base term.

### 7.2 The Chirality Balance Hypothesis

**Prediction:** Chirally balanced moduli ($I = 0$) exhibit stronger spectral localization (tighter block structure, larger gaps) than unbalanced moduli.

**Observation:** At $\alpha = \sqrt{2/3}$:

| $m$ | $I$ | Eigenvalue variance | IPR ratio |
|---|---|---|---|
| 30 | 0 (balanced) | 0.203 (highest) | 1.43 (lowest) |
| 210 | 1 | 0.037 | 1.92 |
| 2310 | 2 | 0.004 | 2.46 |
| 30030 | 1 | 0.0003 (lowest) | 2.88 (highest) |

The balanced modulus $m = 30$ has the *highest* eigenvalue variance and the *lowest* IPR ratio — it is the most dispersive, not the most localized. Localization increases monotonically with $\varphi(m)$ regardless of chirality. The correlation between $I$ and any spectral indicator is zero after controlling for $\varphi$.

**Why it fails:** The flat-band condensation (Section 4) is a dimensionality effect: as $\varphi(m) \to \infty$, the normalized operator $H(\alpha)$ has $\mathcal{O}(\varphi)$ eigenvalues that must fit between the Perron eigenvalue and the anti-Perron eigenvalue, forcing the variance to scale as $1/\varphi$. Chirality balance affects the *convergence rate* of quantities like the alignment $A \to 3/4$ (as shown by the Chebyshev Cooling Law in [1]), but at fixed $\varphi$, the dominant effect is the $1/\varphi$ scaling, which dwarfs any chirality correction.

### 7.3 The $\sqrt{3/2}$ Conjecture

**Prediction:** The critical shear amplitude is exactly $\alpha_c = \sqrt{3/2} = 1/\kappa$, where $\kappa = \sqrt{2/3}$ is the irreconcilability coupling constant, so that $\kappa \cdot \alpha_c = 1$ exactly.

**Observation:** A fine-resolution sweep (222 $\alpha$-values) showed the crash point appearing to plateau between $m = 2310$ and $m = 30030$ at $+1.13\%$ above $\sqrt{3/2}$. However, a micro-sweep (Section 6.3) revealed this "plateau" was a grid resolution artifact — the transition width ($\sim 4 \times 10^{-6}$) was 500 times narrower than the grid spacing, causing both moduli to crash between the same two grid points.

The true crash at $m = 2310$ (micro-sweep) is $\alpha_c = 1.2389219$, which is $+1.16\%$ above $\sqrt{3/2}$ — still decisively falsified. The candidate closed form is now $\sqrt{135/88}$ (Conjecture 6.1), with error $+0.027\%$ at $m = 2310$ and $+0.015\%$ at $m = 30030$.

**Why it fails:** The product $\kappa \cdot \sqrt{135/88} = \sqrt{45/44} \approx 1.0114$, not $1.0$. The coupling constant $\kappa = \sqrt{2/3}$ is not the simple reciprocal of the critical shear. The critical value encodes additional arithmetic structure — the factorization $135/88 = 3^3 \cdot 5 / (2^3 \cdot 11)$ involves prime powers beyond what appears in $\kappa = \sqrt{2/3}$ alone.

### 7.4 What Worked

Despite the failures, the chirality balance hypothesis got three structural elements right:

1. **The operator $H(\alpha)$ is well-defined and Hermitian.** The construction is mathematically sound: the skew-symmetry of $[D, P_\tau]$ guarantees purely real eigenvalues.
2. **The flat-band condensation is real.** The prediction that eigenvalues would collapse onto a small number of modes was correct — it's a genuine spectral feature of the primorial limit.
3. **The motivating question was productive.** Asking "where is the spectral phase transition?" led to the discovery of the $\alpha_c \approx \sqrt{135/88}$ transition, even though it was at the wrong $\alpha$.

---

## 8. Sylow Torsion Defects and the Mod-4 Partition

The falsified predictions of Section 7 leave a residual question: the Double Helix Parity Law [1] partitions primes by their class modulo 4, and the chirality balance hypothesis built on that partition. The hypothesis failed — but does the mod-4 partition itself carry any genuine arithmetic signature in the spectrum of $H(\alpha)$?

### 8.1 The Odd/Even Block Partition

Throughout this paper, the coprime residues $\mathcal{R}$ are implicitly partitioned by the map $r \mapsto r \bmod 4$: the "odd block" consists of residues $r \equiv 1 \pmod{4}$ and the "even block" consists of residues $r \equiv 3 \pmod{4}$. (Since $2 \mid m$ for all primorials $m \geq 6$, every coprime residue is odd, and the mod-4 partition is exhaustive.)

This partition plays a central role in the Double Helix Parity Law [1] and in the block structure of $D_{\text{sym}}$.

### 8.2 Non-Multiplicativity of the Mod-4 Partition

A critical subtlety: the map $\chi_{-1} : r \mapsto (-1)^{(r-1)/2}$ (equivalently, $+1$ if $r \equiv 1 \pmod{4}$, $-1$ if $r \equiv 3 \pmod{4}$) is **not** a multiplicative Dirichlet character [7] on $(\mathbb{Z}/m\mathbb{Z})^*$ when $4 \nmid m$.

For all primorials $m = 2 \cdot 3 \cdot 5 \cdots p_k$, we have $m \equiv 2 \pmod{4}$, so $4 \nmid m$. Direct computation confirms non-multiplicativity:

| $m$ | Counterexample |
|---|---|
| 30 | $7 \times 7 \equiv 19 \pmod{30}$, but $\chi_{-1}(7)^2 = 1 \neq \chi_{-1}(19) = -1$ |
| 210 | $11 \times 23 \equiv 43 \pmod{210}$, but $\chi_{-1}(11)\chi_{-1}(23) = 1 \neq \chi_{-1}(43) = -1$ |
| 2310 | $29 \times 83 \equiv 97 \pmod{2310}$, but $\chi_{-1}(29)\chi_{-1}(83) = -1 \neq \chi_{-1}(97) = +1$ |

This non-multiplicativity has a consequence: arguments that rely on the CRT factorization $\chi_{-1} = \chi_{-1}^{(m_{\text{old}})} \otimes \chi_{-1}^{(p)}$ to study block structure when passing from $m$ to $m \cdot p$ are invalid. The mod-4 partition is not a character of the group — it is a "phantom character" that depends on the integer representative, not just the residue class.

### 8.3 Sylow Torsion and Bragg Scattering

The multiplicative orders $\text{ord}(r, m)$ of coprime residues carry $q$-adic structure. For a prime $q$, define $v_q(\text{ord}(r,m))$ as the $q$-adic valuation of the order. The $q$-torsion defect is

$$\Delta v_q(m) = \sum_{\substack{r \in \mathcal{R} \\ v_q(\text{ord}(r)) > 0 \\ r \equiv 1 \bmod 4}} 1 \;-\; \sum_{\substack{r \in \mathcal{R} \\ v_q(\text{ord}(r)) > 0 \\ r \equiv 3 \bmod 4}} 1$$

i.e., the difference in the count of $q$-torsion-carrying residues between the odd and even blocks.

**Computed defects:**

| $m$ | Terminal prime | $\Delta v_3$ | $\Delta v_5$ | $\Delta v_7$ |
|---|---|---|---|---|
| 30 | 5 | 0 | — | — |
| 210 | 7 | $-4$ | 0 | — |
| 2310 | 11 | 0 | $+4$ | 0 |
| 30030 | 13 | 0 | 0 | 0 |

The defect injection obeys a **Bragg scattering condition**: a new prime $p$ entering the primorial injects $q$-torsion asymmetrically into the mod-4 blocks if and only if $p \equiv 1 \pmod{q}$.

At $m = 2310$ (terminal prime 11), the condition $11 \equiv 1 \pmod{5}$ triggers a 5-adic defect $\Delta v_5 = +4$. The full order-statistic histogram:

| $v_5(\text{ord})$ | Odd block | Even block | Difference |
|---|---|---|---|
| 0 | 46 | 50 | $-4$ |
| 1 | 194 | 190 | $+4$ |

This is exact arithmetic on the group structure — no statistical sampling.

### 8.4 The Washing Conjecture: Partially Falsified

An algebraic argument (based on the parity cross-product of CRT components) predicted that expanding a primorial from $m$ to $m \cdot p$ should perfectly wash out all historical $q$-torsion defects, leaving only freshly injected defects from the terminal prime.

**Prediction for $m = 510510$ ($= 2 \cdot 3 \cdot 5 \cdot 7 \cdot 11 \cdot 13 \cdot 17$, terminal prime 17):** Since $17 \not\equiv 1 \pmod{q}$ for any $q \in \{3, 5, 7, 11, 13\}$, the prediction was $\Delta v_q = 0$ for all odd primes $q$ — a "perfectly balanced" primorial.

**Computed result ($\varphi = 92{,}160$, exact arithmetic):**

| $q$ | $\Delta v_q$ predicted | $\Delta v_q$ actual |
|---|---|---|
| 3 | 0 | $+16$ |
| 5 | 0 | $+16$ |
| 7 | 0 | 0 |
| 11 | 0 | 0 |
| 13 | 0 | 0 |

The 3-adic and 5-adic defects persist with $\Delta v_q = 16$ despite 17 injecting no new torsion. The washing conjecture is **falsified**.

**Root cause.** The washing argument assumes $\chi_{-1}$ is a multiplicative character that factors via CRT. Since $\chi_{-1}$ is the mod-4 partition and $4 \nmid m$, it does NOT factor (Section 8.2). When $m \to m \cdot p$, each parent residue $a$ (mod $m$) spawns $p-1$ surviving offspring $a + km$ for $k \in \{0, \ldots, p-1\} \setminus \{k^*\}$, where $k^*(a) \equiv -m^{-1} a \pmod{p}$ is the killed offspring. Since $m \equiv 2 \pmod{4}$, even-$k$ offspring preserve the parent's mod-4 class while odd-$k$ offspring flip it. The parity of $k^*(a)$ depends on $a \bmod p$, creating a systematic tilt: parents in one half-plane of $(\mathbb{Z}/p\mathbb{Z})^*$ preferentially lose same-parity offspring, while parents in the other half lose opposite-parity offspring.

The defect magnitude $16 = p - 1 = 17 - 1$ is the size of the CRT fiber.

### 8.5 Connection to the Phase Transition

The overnight batch (Section 13) resolves the question of whether Sylow torsion defects affect the critical shear $\alpha_c$. Comparing $m = 15015 = 3 \cdot 5 \cdot 7 \cdot 11 \cdot 13$ with $m = 30030 = 2 \cdot 3 \cdot 5 \cdot 7 \cdot 11 \cdot 13$ (both $\varphi = 5760$, different Sylow structure):

- **The crash positions are effectively identical**: $\alpha_c(15015) = 1.2387689$ vs $\alpha_c(30030) = 1.2387666$ ($\delta = +0.00018\%$).
- **The transition widths differ by $80\times$**: $8.0 \times 10^{-6}$ (non-primorial) vs $1.0 \times 10^{-7}$ (primorial).

Sylow torsion does not shift the transition — it **smears** it. The critical shear $\alpha_c$ is controlled by the matrix dimension $\varphi$, while the transition sharpness encodes arithmetic quality. Primorials, which include every small prime and have the most uniform coprime residue distribution, produce the narrowest transitions.

**Remark 8.1** (Torsion defect coincidence: falsified). At $m = 30030$, the 3-adic torsion defect $\Delta v_3 = 1$ produces a $2$:$1$ partition of residues between mod-4 classes. This coincidentally matches the $T_i/T_e = 2$ ion-to-electron pedestal temperature ratio observed at the JET L-H transition (Andrew et al., 2003 [9]). However, the coincidence is falsified: Ryter et al. (2014) [10] showed that in pure ECRH plasmas on ASDEX Upgrade (electron heating only, no NBI), $T_e > T_i$ at the L-H transition. The $T_i/T_e = 2$ ratio is an NBI heating artifact, not a fundamental property of the confinement transition. The numerical match between the torsion defect and the temperature ratio is therefore accidental.

---

## 9. Discussion

### 9.1 The Active Transport Operator as a New Object

To our knowledge, $H(\alpha) = D/\lambda_P + i\alpha [D, P_\tau]/\lambda_P$ has not appeared previously in the number-theoretic or spectral theory literature. The operator combines the additive geometry of the coprime residues (via $D_{\text{sym}}$) with their multiplicative structure (via the inversion permutation $P_\tau$) into a single Hermitian matrix whose spectrum encodes the interaction between these two fundamental arithmetic operations.

The flat-band condensation and the gap position phase transition suggest that $H(\alpha)$ has a well-defined thermodynamic limit as $m \to \infty$ through primorials. The spectral measure should converge weakly to a distribution with atoms at $\{+1, \lambda_-, 0\}$ with weights $\{1/\varphi, 1/\varphi, 1 - 2/\varphi\}$.

**Guide to this section.** Section 9.2 documents the full spectral architecture at $\alpha_c$, including the stable Perron–anti-Perron gap $\Delta\lambda = 1.7441 \pm 3.7 \times 10^{-4}$ (203 $\alpha$-values). Section 9.3 presents the Landau-Zener connection — a dynamical interpretation of the narrow transition width. Section 9.4 tests whether the anti-Perron chiral asymmetry is genuine via the even/odd modulus nozzle test. Section 9.5 discusses the relation to plasma physics, distinguishing falsified predictions from open questions. The perturbation response experiments and the CRT block architecture are presented separately in Section 10. Open problems are collected in Section 12.

### 9.2 Spectral Architecture at the Crash Point

Full eigendecomposition at $\alpha_c(30030) = 1.2387666$ reveals a three-component spectral architecture:

| Component | Eigenvalue | Count | Spectral weight | Chiral asymmetry $\mathcal{A}$ |
|---|---|---|---|---|
| Perron mode | $\lambda_P = +1.0000$ | 1 | $1/\varphi$ | $-2.5 \times 10^{-6}$ (balanced) |
| Anti-Perron mode | $\lambda_- = -0.7441$ | 1 | $1/\varphi$ | $-0.00296$ (even-dominant) |
| Dispersive bulk | $0.01 < \lvert\lambda\rvert < 0.34$ | 13 | $0.2\%$ | — |
| Flat band | $\lvert\lambda\rvert < 0.01$ | 5745 | $99.7\%$ | — |

The four components sum to $1 + 1 + 13 + 5745 = 5760 = \varphi(30030)$, exhausting the full Hilbert space. The 13 dispersive modes include the $\sigma$-odd sector boundary states (max eigenvalue $0.337$, min $-0.742$; Section 11.2) that sit between the flat-band edge and the anti-Perron mode.

**Remark 9.1** (Anti-Perron eigenvalue is $\alpha$-dependent). The value $\lambda_- = -0.7441$ at $\alpha_c$ differs from the $\lambda_- \approx -0.591$ reported at $\alpha = \sqrt{2/3}$ in Section 4.2. The anti-Perron eigenvalue shifts continuously as the shear parameter increases; the gap $\Delta\lambda = \lambda_P - \lambda_-$ widens from $1.591$ at $\alpha = \sqrt{2/3}$ to $1.7441$ at $\alpha_c$.

The Perron gap (distance from $\lambda_P$ to the next eigenvalue) is $\Delta_P = 0.661$, while the anti-Perron gap (distance from $\lambda_-$ to the next eigenvalue above it) is only $0.0003$ — the anti-Perron mode is nearly degenerate with the flat-band edge.

**Observation 9.2** (Perron–anti-Perron gap stability). The eigenvalue gap between the Perron and anti-Perron modes is a **structural invariant of $\varphi$, independent of primorial structure**:

| $m$ | Type | $\varphi$ | $\Delta\lambda = \lambda_P - \lambda_-$ | Stability ($\sigma$ across crash window) |
|---|---|---|---|---|
| 30030 | primorial | 5760 | $1.7441 \pm 3.7 \times 10^{-4}$ | 203 $\alpha$ values across $[\sqrt{3/2},\, \alpha_c]$; $0.021\%$ relative variation |
| 15015 | non-primorial | 5760 | $1.7441569 \pm 1.5 \times 10^{-5}$ | 5 points, slight drift |

The gap differs by $0.0012\%$ between moduli — both within rounding of a common value $\approx 1.7441$. The gap is dimension-determined ($\varphi = 5760$), not torsion-sensitive.

**Observation 9.3** (Active parity compensation in the Perron state). The Perron eigenvector is nearly uniform — a theorem-level consequence of flat-band condensation, since $D_{\text{sym}}/\lambda_P$ is a rank-1 perturbation of the flat band. However, it is not *perfectly* flat. A perfectly uniform vector $\mathbf{1}/\sqrt{\varphi}$ on the $m = 30030$ lattice ($|\mathcal{R}_1| = 2879$, $|\mathcal{R}_3| = 2881$) would produce chiral asymmetry $\mathcal{A}_{\text{flat}} = -2/5760 = -3.47 \times 10^{-4}$. The measured value is $\mathcal{A} = -2.5 \times 10^{-6}$ — a factor of $\sim 140$ smaller. Since float64 eigensolvers carry precision to $\sim 10^{-15}$, this discrepancy is not numerical noise; it is a signal that the Perron state **actively modulates its amplitudes** to compensate for the coprime sieve's chiral imbalance. The eigenvector slightly increases its weight on the minority $1$-block residues, driving $\mathcal{A}$ two orders of magnitude closer to zero than geometric uniformity would predict. Flat-band condensation forces *near*-uniformity, but within that near-uniform envelope, the Perron state actively resists the lattice's geometric defect to preserve approximate parity balance. The mechanism is Perron--Frobenius: the leading eigenvector of a positive matrix maximizes $\langle \mathbf{1} | v \rangle$, which penalizes deviations from uniformity. The residual asymmetry $\mathcal{A} \sim 10^{-6}$ reflects the finite gap between the Perron eigenvalue and the flat band, which controls how tightly the eigenvector is pinned to uniformity.

### 9.3 The Landau-Zener Connection

The narrow transition width ($\delta \sim 10^{-7}$) and the stable gap ($\Delta\lambda = 1.7441$) invite a dynamical interpretation via the Landau-Zener formula [12]. If $\alpha$ is promoted from a static parameter to a time-dependent sweep $\alpha(t) = \alpha_0 + vt$, and if a physical system exists whose Hamiltonian is $H(\alpha(t))$, then the non-adiabatic transition probability at the avoided crossing is

$$P_{LZ} = \exp\left(-\frac{2\pi \delta^2}{\hbar \, |v|}\right)$$

where $\delta$ is the minimum gap at the crossing and $v = d\alpha/dt$ is the sweep rate. With $\delta \sim 10^{-7}$, even extremely slow sweeps ($|v| \gg \delta^2$) yield $P_{LZ} \to 1$ — near-complete non-adiabatic transfer from the symmetric Perron state to the bulk/anti-Perron manifold.

The anti-Perron asymmetry $\mathcal{A}_{\text{anti}} = -0.003$ would then determine the spatial bias of the transferred population. This is the honest version of the "topological exhaust" mechanism: not a resonant pump, but a parametric sweep through a sharp avoided crossing. The sweep-rate insensitivity (any $v$ works, because $\delta$ is so small) would make the mechanism robust to engineering tolerances.

**Critical caveat.** This interpretation requires:
1. Identifying $\alpha$ with a time-dependent physical parameter.
2. Identifying the physical system whose Hamiltonian *is* $H(\alpha)$, not merely *resembles* it.
3. Neither has been established. The Landau-Zener connection is stated as a mathematical observation about the parametric structure of $H(\alpha)$, not as a physical claim.

### 9.4 The Nozzle Test: Anti-Perron Asymmetry at Even vs. Odd Modulus

The anti-Perron chiral asymmetry $\mathcal{A}_{\text{anti}} = -0.00296$ at $m = 30030$ (Section 9.2) measures mass imbalance between residues $\equiv 1 \pmod{4}$ and residues $\equiv 3 \pmod{4}$. At $m = 30030$ (even), every coprime residue is odd, so the mod-4 partition $\{r \equiv 1\} \cup \{r \equiv 3\}$ covers **all** $\varphi = 5760$ residues.

At $m = 15015$ (odd), coprime residues can be even. The partition structure differs fundamentally:

| Class | $m = 30030$ (even) | $m = 15015$ (odd) |
|---|---|---|
| $r \equiv 0 \pmod{4}$ | 0 (impossible) | 1440 (25.0%) |
| $r \equiv 1 \pmod{4}$ | 2879 (50.0%) | 1440 (25.0%) |
| $r \equiv 2 \pmod{4}$ | 0 (impossible) | 1440 (25.0%) |
| $r \equiv 3 \pmod{4}$ | 2881 (50.0%) | 1440 (25.0%) |

The initial nozzle test (`nozzle_test.py`) showed a sign flip: $\mathcal{A}_{\text{anti}} = -0.00296$ at $m = 30030$ vs $+0.00227$ at $m = 15015$. The critical concern was whether this reflected a genuine change or a partition artifact, since the mod-4 partition covers 100% of residues at even $m$ but only 50% at odd $m$.

#### 9.4.1 4-Class Resolution

The full 4-class partition analysis (`partition_4class.py`) resolves this definitively. Anti-Perron eigenvector mass across all four mod-4 classes:

| Class | $m = 30030$ count | $m = 30030$ mass | $m = 15015$ count | $m = 15015$ mass | $m = 15015$ enrichment |
|---|---|---|---|---|---|
| $r \equiv 0 \pmod{4}$ | 0 | — | 1440 | 24.886% | 0.9955 |
| $r \equiv 1 \pmod{4}$ | 2880 | 49.852% | 1440 | 25.114% | 1.0045 |
| $r \equiv 2 \pmod{4}$ | 0 | — | 1440 | 25.114% | 1.0045 |
| $r \equiv 3 \pmod{4}$ | 2880 | 50.148% | 1440 | 24.886% | 0.9955 |

**The sign flip is genuine.** Within the $\{1, 3\}$ sector alone — the only comparison available at both moduli — the sector-normalized asymmetry flips:

| | $m = 30030$ | $m = 15015$ |
|---|---|---|
| $\mathcal{A}$ (total-normalized) | $-0.00296$ | $+0.00227$ |
| $\mathcal{A}$ (sector-normalized, within $\{1,3\}$) | $-0.00296$ | $+0.00454$ |
| $\{1,3\}$-sector coverage | 100% | 50% |

The sector normalization actually **amplifies** the sign flip (from $+0.002$ to $+0.005$), ruling out the dilution artifact hypothesis.

#### 9.4.2 Klein Four-Group Mirror Symmetry

The 4-class data at $m = 15015$ exhibits an exact pairing:

$$\text{enrichment}(c) = \text{enrichment}(c + 2 \bmod 4)$$

i.e., classes $\{0, 3\}$ are depleted identically (0.9955) and classes $\{1, 2\}$ are enriched identically (1.0045). The anti-Perron eigenvector's mod-4 mass respects the Klein four-group $\mathbb{Z}/2 \times \mathbb{Z}/2$ structure: the $(+2 \bmod 4)$ involution is an exact symmetry of the mass distribution. This is **not** a parity effect — Analysis 2 confirms that even and odd residues carry exactly 50.000%/50.000% of the total mass. The asymmetry is purely *within* each parity sector.

#### 9.4.3 Partition-Independent Confirmation

Spatial autocorrelation of $|v_r|^2$ along the residue ring — a partition-independent observable — is essentially identical at both moduli:

| Lag | $m = 30030$ | $m = 15015$ |
|---|---|---|
| 1 | $+0.833$ | $+0.815$ |
| 2 | $+0.828$ | $+0.819$ |
| 5 | $+0.831$ | $+0.832$ |
| 10 | $+0.827$ | $+0.820$ |

The mass clustering structure is identical. Only the *direction* of the mod-4 bias reverses.

#### 9.4.4 Established Findings

1. **The sign flip is genuine.** Sector-normalized asymmetry reverses from $-0.003$ to $+0.005$; not a partition artifact.
2. **Klein four-group symmetry.** Anti-Perron mass at $m = 15015$ pairs classes $\{0 \leftrightarrow 3\}$ and $\{1 \leftrightarrow 2\}$ exactly.
3. **Parity is exactly conserved.** Even/odd residues carry 50.000%/50.000% — the asymmetry is intra-parity, not inter-parity.
4. **The Perron mode is flat.** Enrichment $\approx 1.0000$ at all classes, both moduli.
5. **Spatial autocorrelation preserved.** Mass clustering ($\rho \approx 0.82$) is independent of the sign flip.
6. Removing the factor of 2 from $m$ acts as a **topological switch** on the anti-Perron exhaust direction: the 2-torsion in the multiplicative group $\mathbb{Z}/m\mathbb{Z}^\times$ controls the sign of $\mathcal{A}$.

### 9.5 Relation to Plasma Physics: Falsified Predictions and Open Questions

The chirality balance hypothesis (Section 1.2) was originally framed in the language of tokamak physics: the primorial hierarchy as a model for plasma confinement regimes, the commutator $[D, P_\tau]$ as drift-wave turbulence. Several specific predictions tested in this paper have been falsified.

**Structural obstruction.** The fundamental obstacle to any direct plasma mapping is the non-locality of the multiplicative inversion permutation $P_\tau$: the map $r \to r^{-1} \bmod m$ sends a residue to its multiplicative inverse, which can be anywhere on the ring. Plasma turbulence cascades via local triad interactions ($\vec{k}_1 + \vec{k}_2 = \vec{k}_3$), where modes couple to their neighbors in wavenumber space. The operator $H(\alpha)$ computes modular inverses; the Navier-Stokes equations do not. This structural difference — non-local algebraic coupling vs. local convective coupling — constrains the class of physical systems that could realize $H(\alpha)$.

**Eigenvector predictions.** The flat-band condensation (Section 4) forces the Perron eigenvector to be uniform — $100\%$ DC power, $0.0\%$ spectral power at every tested period, chiral asymmetry $\mathcal{A} = -2.5 \times 10^{-6}$ (active parity compensation, not a flat-vector artifact; see Observation 9.3). Any prediction that requires the Perron eigenvector to encode spatial structure (coil winding angles, exhaust asymmetry, periodic modulation) is theorem-level impossible: flat-band condensation guarantees uniformity. The anti-Perron eigenvector does carry structure ($\mathcal{A}_{\text{anti}} = -0.003$, Klein four-group symmetry), but this mode represents the most *negative* eigenvalue, not the dominant transport mode.

**Torsion defect coincidence.** The 3-adic torsion defect $\Delta v_3 = 1$ at $m = 30030$ coincidentally matches the $T_i/T_e = 2$ ion-electron temperature ratio observed at the JET L-H transition (Andrew et al., 2003 [9]). However, Ryter et al. (2014) [10] performed the decisive test on ASDEX Upgrade: in pure ECRH plasmas (electron heating only, no NBI), $T_e > T_i$ at the L-H transition — the ratio is *inverted*. The edge ion heat flux $Q_{i,\text{edge}}$ governs the transition, but the $T_i/T_e = 2$ value is an NBI heating artifact, not a fundamental property of the L-H transition. The 3-adic coincidence is **falsified**.

The Active Transport Operator $H(\alpha)$ is a genuinely new object in algebraic spectral theory. Its mathematical properties (flat-band condensation, the sharp phase transition, torsion defects, Klein four-group symmetry) stand on their own. Whether a physical system exists whose transport operator *is* $H(\alpha)$ remains an open question — but that system, if it exists, would need to implement modular inversion as a physical operation.

---

## 10. Perturbation Response and the CRT Block Architecture

### 10.1 Selective Screening of Prime-Periodic Perturbations

The flat-band condensation (Section 4) naturally raises the question of robustness: does the flat band survive perturbation? We treat $H(\alpha)$ as a tight-binding Hamiltonian and study the time-averaged Loschmidt echo [13]

$$\langle L \rangle_\infty = \lim_{T \to \infty} \frac{1}{T} \int_0^T |\langle r{=}1 | e^{-iHt} | r{=}1 \rangle|^2 \, dt = \sum_k |c_k|^4$$

where $c_k = \langle k | r{=}1 \rangle$ are the overlaps of the site state $|r{=}1\rangle$ with the eigenstates of $H$. For the unperturbed operator at $m = 2310$ ($\varphi = 480$), the flat-band condensation concentrates spectral weight so that $\langle L \rangle_\infty \approx 384 / \varphi$, far above the ergodic baseline $1/\varphi$.

We perturb the Hamiltonian as $H_\varepsilon = H(\alpha) + \varepsilon V$, where $V$ is normalized to $\|V\|_F = 1$, and compare three classes of perturbation:

1. **Random Hermitian noise.** $V = (R + R^\dagger)/2$ with i.i.d. complex Gaussian entries.
2. **Prime-periodic diagonal potential.** $V_p = \text{diag}(\cos(2\pi r_i / p))$, normalized.
3. **Prime-periodic hopping modulation.** $(V_p^{\text{hop}})_{ij} = \cos(2\pi(r_i - r_j)/p)$, normalized.

#### 10.1.1 Random Noise Fragility

Under random Hermitian perturbation, the Loschmidt echo decays rapidly:

| $\varepsilon$ | $\varepsilon / \|A\|_F$ | Survival $\langle L \rangle_\varepsilon / \langle L \rangle_0$ |
|---|---|---|
| 0.01 | 0.9% | 81.3% |
| 0.03 | 2.6% | 35.6% |
| 0.10 | 8.7% | 11.3% |
| 0.30 | 26% | 4.0% |
| 1.00 | 87% | 1.4% |

At $\varepsilon = 0.03$ (less than $3\%$ of the operator norm), two-thirds of the caging is destroyed. The flat band is **fragile** against unstructured noise. This rules out any claim of topological protection in the conventional sense.

#### 10.1.2 Prime-Periodic Perturbations and Selective Screening

The situation changes qualitatively when the perturbation has arithmetic structure. At $m = 2310 = 2 \cdot 3 \cdot 5 \cdot 7 \cdot 11$ and $\alpha = \alpha_c \approx \sqrt{135/88}$, we apply the diagonal potential $V_p = \text{diag}(\cos(2\pi r / p))$ for primes $p$ both in and out of the factorization of $m$:

| Perturbation | Type | Survival at $\varepsilon = 0.03$ | Survival at $\varepsilon = 0.1$ | Survival at $\varepsilon = 1.0$ |
|---|---|---|---|---|
| Random Hermitian | unstructured | 35.6% | 11.3% | 1.4% |
| $V_{13}$ (diagonal) | missing prime | 93.8% | 81.3% | 66.8% |
| $V_{17}$ (diagonal) | missing prime | 96.8% | 90.5% | 74.4% |
| $V_{19}$ (diagonal) | missing prime | 93.8% | 85.0% | 68.7% |
| $V_{11}$ (diagonal) | in factorization | 96.8% | 93.3% | 93.6% |
| $V_7$ (diagonal) | in factorization | 96.2% | 91.7% | 90.7% |
| $V_5$ (diagonal) | in factorization | 95.8% | 91.7% | 90.6% |
| $V_3$ (diagonal) | in factorization | 100.0% | 100.0% | 100.0% |
| $V_2$ (diagonal) | in factorization | 100.0% | 100.0% | 100.0% |

Prime-periodic noise is $10$–$20\times$ less destructive than random noise at every tested $\varepsilon$. Within the prime-periodic class, the lattice **selectively screens its own primes**: in-factorization perturbations are less destructive than missing-prime perturbations. The aggregate comparison:

| $\varepsilon$ | Avg survival (IN primes) | Avg survival (OUT primes) | Ratio IN/OUT | Random |
|---|---|---|---|---|
| 0.03 | 97.8% | 94.8% | 1.031 | 35.6% |
| 0.10 | 95.3% | 85.6% | 1.114 | 11.3% |
| 0.30 | 95.3% | 73.8% | 1.291 | 4.0% |
| 1.00 | 95.0% | 70.0% | 1.357 | 1.4% |

**Observation 10.1** (Screening ratio amplification at $\alpha_c$). At $\alpha = 0$, the IN/OUT screening ratio at $\varepsilon = 1.0$ is $1.147$. At $\alpha = \alpha_c$, it rises to $1.357$. The critical point does not control flat-band caging (which is universal across $\alpha$), but it **amplifies the discrimination** between in-factorization and missing-prime perturbations.

#### 10.1.3 Perfect Screening of $p = 2$ and $p = 3$

The perturbations $V_2 = \text{diag}(\cos(\pi r))$ and $V_3 = \text{diag}(\cos(2\pi r / 3))$ produce **exactly $100\%$ survival** at all tested $\varepsilon \in [0.001, 1.0]$, at both $\alpha = 0$ and $\alpha = \alpha_c$.

The mechanism is arithmetic, not dynamical. Every coprime residue $r \in \mathcal{R}$ satisfies $\gcd(r, m) = 1$, so $r$ is odd (not divisible by 2) and not divisible by 3. Therefore:

- $\cos(\pi r) = \cos(\pi \cdot \text{odd}) = -1$ for all $r \in \mathcal{R}$. The perturbation $V_2$ is a constant on the Hilbert space — a scalar shift that commutes with $H$.
- $\cos(2\pi r / 3)$: since $3 \mid m$ and $\gcd(r, m) = 1$, we have $r \not\equiv 0 \pmod{3}$, so $r \bmod 3 \in \{1, 2\}$. By the CRT structure, these two classes are equally populated and $V_3$ has exactly two distinct values on $\mathcal{R}$. But since $V_3$ commutes with the block structure of the coprime sieve, it shifts eigenvalues without mixing eigenstates.

The coprime sieve removes the Fourier components at frequencies $2\pi/2$ and $2\pi/3$ from the Hilbert space by construction. This is a direct consequence of the Permanent Irreconcilability Theorem [1]: the number 3 does not merely create the topology of $H(\alpha)$ — it creates a noise firewall. Any perturbation whose spatial frequency is commensurate with a factor of $m$ is projected out of the relevant subspace.

#### 10.1.4 Hopping Perturbations

Off-diagonal (hopping) perturbations $V_p^{\text{hop}}$ with $(V_p^{\text{hop}})_{ij} = \cos(2\pi(r_i - r_j)/p)$ behave qualitatively differently from diagonal perturbations. At $\varepsilon = 1.0$ and $\alpha = \alpha_c$:

| Perturbation | Survival |
|---|---|
| $V_{11}^{\text{hop}}$ (in factorization) | 99.2% |
| $V_{13}^{\text{hop}}$ (missing prime) | 99.1% |

Hopping perturbations preserve $\sim 99\%$ of the Loschmidt echo regardless of whether $p$ divides $m$. The IN/OUT distinction vanishes. The mechanism is that hopping modulations couple to the *off-diagonal* structure of the flat band, which has a different symmetry class than diagonal potentials.

#### 10.1.5 First-Order Perturbation Theory

The structural analysis confirms the screening mechanism. The first-order Perron eigenvalue shift $\delta\lambda_P = \langle \psi_P | V | \psi_P \rangle$ under each perturbation:

| Perturbation | $\|P_{\text{flat}} V P_{\text{flat}}\| / \|V\|$ | $\delta\lambda_P$ |
|---|---|---|
| $V_2$ [IN] | 0.984 | $-0.046$ |
| $V_3$ [IN] | 0.984 | $-0.046$ |
| $V_5$ [IN] | 0.971 | $-0.019$ |
| $V_7$ [IN] | 0.970 | $-0.012$ |
| $V_{11}$ [IN] | 0.969 | $-0.007$ |
| $V_{13}$ [OUT] | 0.969 | $+0.000$ |
| $V_{17}$ [OUT] | 0.969 | $+0.000$ |
| $V_{19}$ [OUT] | 0.969 | $-0.000$ |
| Random | 0.969 | $-0.001$ |

All perturbations project $\sim 97\%$ of their weight onto the flat band. However, the first-order Perron shift reveals a clear dichotomy: **in-factorization primes shift the Perron eigenvalue** ($\delta\lambda_P$ between $-0.007$ and $-0.046$, scaling as $\sim 1/p$), while **missing primes leave it unchanged** ($|\delta\lambda_P| < 0.001$). The in-factorization perturbations couple to the Perron mode and shift its energy, but do not mix it with the flat band — this is the screening mechanism. Missing primes neither couple to the Perron mode nor are screened by the sieve, so they attack the flat-band degeneracy directly.

#### 10.1.6 Established Findings

**(5) Selective screening.** The flat band of $H(\alpha)$ is fragile against unstructured noise but selectively robust against prime-periodic perturbations. At $m = 2310$, $\alpha = \alpha_c$, $\varepsilon = 1.0$:

| Noise class | Survival | Mechanism |
|---|---|---|
| Random Hermitian | 1.4% | Destroys flat-band degeneracy |
| Missing primes ($p = 13, 17, 19$) | 70.0% (avg) | Incommensurate; partially couples to flat band |
| In-factorization primes ($p = 5, 7, 11$) | 91.6% (avg) | Commensurate; Perron shift without flat-band mixing |
| Core primes ($p = 2, 3$) | 100.0% | Coprime sieve eliminates these Fourier components entirely |
| Hopping modulations | 99.1% | Off-diagonal; different symmetry class |

**(6) Screening amplification at $\alpha_c$.** The IN/OUT ratio grows from $1.15$ at $\alpha = 0$ to $1.36$ at $\alpha = \alpha_c$ — the critical shear enhances the lattice's ability to discriminate between its own primes and foreign ones.

**(7) Coprime sieve as noise firewall.** The perfect screening of $p = 2$ and $p = 3$ is a theorem-level consequence of the coprime residue construction: the Hilbert space of $H(\alpha)$ is built on $\{r : \gcd(r, m) = 1\}$, which excludes all multiples of every prime factor of $m$. Perturbations at these frequencies are constant (or nearly so) on the surviving residues and therefore commute with $H$.

### 10.2 The Perron–Anti-Perron Gap: Free Precession, Not Rabi Oscillation

The initial measurements in this section suggested a "noise-immune Rabi oscillation" constituting a "primorial qubit." Falsification testing (`primorial_universality_test.py`, `crt_block_hamiltonian.py`, `block9_deep_dive.py`, `the_connection.py`) established that this specific mechanism is wrong: (a) the period-$3.51$ oscillation is **free precession** of the Perron–anti-Perron gap, not a driven rotation; (b) the maximum Rabi transfer is $0.41$, not $1.0$; (c) the CRT blocks form a 16-dimensional qudit, not the qubit basis (the correct qubit degree of freedom — $\sigma$-parity — is identified in Section 11); and (d) the "noise immunity" is explained by the CRT scalar identity ($V_5, V_7$ act as exact scalars within CRT blocks). Kill condition K3 was triggered: all primes show identical period $T = 2\pi/\Delta\lambda = 3.60$ — the signature of undriven precession. The original measurements are preserved below for reproducibility; the CRT block reinterpretation follows in Section 10.3, and the correct qubit structure in Section 11.

#### 10.2.1 The Flat-Band Obstruction

An initial test (`primorial_logic_gate.py`) applied $V_{13} = \text{diag}(\cos(2\pi r / 13))$ at amplitude $A = 0.5$ to the flat-band-projected identity residue state $|\psi_0\rangle = P_{\text{flat}} |r=1\rangle / \|P_{\text{flat}} |r=1\rangle\|$. The fidelity $F(t) = |\langle \psi_0 | e^{-iHt} | \psi_0 \rangle|^2$ dropped from $1.0$ to a minimum of $0.865$ and oscillated with a swing of only $4.6\%$ — not a useful gate rotation.

To diagnose the obstruction, we compute the restriction of $V_{13}$ to the flat-band subspace. Let $U_{\text{flat}}$ be the $\varphi \times n_{\text{flat}}$ matrix of flat-band eigenvectors of $H(\alpha_c)$. The restricted operator $V_{13}^{\text{FB}} = U_{\text{flat}}^\dagger V_{13} U_{\text{flat}}$ is a $473 \times 473$ Hermitian matrix with eigenvalue range $[-0.971, +1.000]$.

Decomposing $|\psi_0\rangle$ in the eigenbasis of $V_{13}^{\text{FB}}$ reveals the obstruction:

| Metric | Value |
|---|---|
| Largest single-eigenstate overlap | $11.5\%$ |
| Top 5 eigenstates total overlap | $36.5\%$ |
| Inverse participation ratio | $23.6 / 473$ |
| $\langle \psi_0 \lvert V_{13} \rvert \psi_0 \rangle$ | $+0.857$ |
| Std dev of $V_{13}$ eigenvalues in $\lvert\psi_0\rangle$ | $0.121$ |

The state $|\psi_0\rangle$ is spread across $\sim 24$ eigenstates, but these are **clustered on a single eigenvalue plateau** at $\lambda \approx +0.885$: the top 15 eigenstates by overlap all share this eigenvalue. The identity residue $r = 1$ satisfies $1 \equiv 1 \pmod{p}$ for all $p$, placing it at a fixed point of the CRT decomposition. Its flat-band projection inherits this arithmetic rigidity, landing on a $V_{13}$ eigenvalue plateau. The standard deviation of $V_{13}$ eigenvalues in $|\psi_0\rangle$ is only $0.121$ across the full range $[-0.971, +1.000]$ — the perturbation acts nearly as a scalar, producing a phase wobble rather than a rotation.

#### 10.2.2 The Perron–Anti-Perron Free Precession

The spectrum of $H(\alpha_c)$ at $m = 2310$ consists of four components: the Perron mode ($\lambda_P = +1.000$), the anti-Perron mode ($\lambda_{\text{AP}} = -0.745$), 5 dispersive bulk states ($0.01 < |\lambda| < 0.34$), and the degenerate flat band ($473$ states at $|\lambda| < 0.01$). Accounting: $1 + 1 + 5 + 473 = 480 = \varphi(2310)$. The Perron–anti-Perron gap is $\Delta\lambda = 1.745$.

The measured oscillation data (preserved for reproducibility):

| Initial state | $F_{\text{avg}}$ | $F_{\min}$ | $F_{\max}$ | Swing | Period |
|---|---|---|---|---|---|
| $P_{\text{flat}} \lvert r=1\rangle$ (flat-band projected) | $0.981$ | $0.955$ | $1.000$ | $0.046$ | — |
| $\lvert\psi_{\text{Rabi}}\rangle$ (optimal flat-band pair) | $0.494$ | $0.000$ | $1.000$ | $1.000$ | $\sim 6.4$ |
| $\lvert\psi_{\text{Perron}}\rangle$ (ground state) | $0.999$ | $0.998$ | $1.000$ | $0.002$ | — |
| $\lvert\psi_{\text{Cavity}}\rangle$ (Perron + anti-Perron) $/\sqrt{2}$ | $0.497$ | $0.000$ | $1.000$ | $1.000$ | $3.51$ |

**Interpretation.** The cavity state period $T = 3.51$ matches the free precession frequency $T_{\text{free}} = 2\pi / \Delta\lambda = 2\pi / 1.745 = 3.60$ to within the measurement window. The universality test (`primorial_universality_test.py`) confirmed that **all** primes — $V_5, V_7, V_{11}, V_{13}, V_{17}, V_{19}$ — produce the identical period $T \approx 3.60$, triggering kill condition K3. This is the hallmark of free precession: the oscillation frequency is set by the energy gap, not by the drive. A true Rabi oscillation would show a drive-dependent frequency $\Omega_R = \sqrt{\Delta^2 + g^2}$ where $g$ varies with the perturbation operator.

An earlier draft described this as "$V_{13}$-mediated coupling between the Perron and anti-Perron modes." The correct description: **the superposition state precesses freely in the energy eigenbasis; the perturbation $V_{13}$ does not drive the rotation, the gap $\Delta\lambda$ does.**

#### 10.2.3 Noise Immunity as CRT Scalar Identity

The original noise immunity result:

| Initial state | Clean swing | Noisy swing | $\max \lvert F_{\text{clean}} - F_{\text{noisy}}\rvert$ |
|---|---|---|---|
| $\lvert\psi_{\text{Cavity}}\rangle$ | $1.000$ | $1.000$ | $< 10^{-3}$ |
| $\lvert\psi_{\text{Rabi}}\rangle$ (flat-band) | $1.000$ | $1.000$ | $0.019$ |

**Interpretation.** The "noise immunity" under $V_5 + V_7$ perturbation is explained by a deeper result: on the coprime residues modulo $m = 2310$, the operators $V_2$ and $V_3$ are exact scalars ($V_2 = -1 \cdot I$, $V_3 = -\tfrac{1}{2} \cdot I$), and the remaining in-factorization operators $V_5, V_7, V_{11}$ are exact scalars **within each CRT block** (spread $= 3.2 \times 10^{-13}$, i.e., numerical zero at eigensolver tolerance). A scalar perturbation $\varepsilon \cdot c \cdot I$ adds a global phase and cannot change any observable. The "screening" is therefore not a dynamical effect but an **algebraic identity**: the CRT decomposition of $(\mathbb{Z}/m\mathbb{Z})^*$ forces in-factorization cosine operators to act trivially within the protected subspaces.

#### 10.2.4 What Remains Valid

Despite the falsification of the "Rabi oscillation" and "primorial qubit" interpretations, the following measurements from Sections 10.1–10.2 are confirmed:

1. **Flat-band obstruction** (Section 10.2.1): the $|r = 1\rangle$ state cannot be rotated. This is real and explained by CRT fixed-point arithmetic.
2. **Selective screening** (Section 10.1): IN-factorization primes are screened more strongly than OUT primes. The screening ratio $1.36\times$ at $\alpha_c$ is genuine, now understood as a consequence of CRT block structure.
3. **Perfect screening of $p = 2, 3$** (Section 10.1.3): $100\%$ survival is exact because these operators are scalars on the entire coprime Hilbert space.
4. **The three-component spectrum** (Perron, anti-Perron, flat band) is a verified structural fact.

What is falsified:
- "Perfect Rabi oscillation" → free precession ($T = 2\pi/\Delta\lambda$)
- "Primorial qubit" (Perron–anti-Perron Rabi gate) → **falsified**; the CRT blocks form a 16-dimensional qudit (max transfer $0.41$), but the true qubit degree of freedom is $\sigma$-parity (Section 11)
- "Noise-immune gate" → CRT scalar identity (algebraic, not dynamical)
- "Arithmetic bandpass filter" metaphor → partially correct but the mechanism is CRT block decomposition, not dynamical screening

#### 10.2.5 Falsification Table

| Original claim (Sections 10.1–10.2) | Status | Evidence |
|---|---|---|
| $V_{13}$ drives Rabi oscillation (period 3.51) | **Falsified** | All primes give identical period $\approx 3.60 = 2\pi/\Delta\lambda$; K3 kill triggered |
| Perron–anti-Perron = logical qubit | **Falsified** | Max Rabi transfer $= 0.41$; CRT blocks form 8-level qudit per $\sigma$-sector (16D pre-projection), not the qubit basis (true qubit is $\sigma$-parity, Section 11) |
| $V_5 + V_7$ noise immunity = dynamical screening | **Falsified** | $V_5, V_7$ are scalars within CRT blocks (spread $= 3.2 \times 10^{-13}$) |
| Flat-band = decoherence-free memory | **Partially valid** | Flat band exists; "decoherence-free" overstated — fragile to random noise (Section 10.1.1) |
| Selective screening ratio $1.36\times$ at $\alpha_c$ | **Valid** | Confirmed; mechanism is CRT block structure |
| $p = 2, 3$ perfect screening | **Valid** | Exact; $V_2 = -I$, $V_3 = -\tfrac{1}{2}I$ on coprimes (proven in `the_connection.py`) |

### 10.3 CRT Block Qudit Architecture

The falsification of the primorial qubit (Section 10.2) led directly to the correct structure: a **CRT block qudit architecture** in which the prime factorization of $m$ organizes the Hilbert space into arithmetically protected subspaces.

**A note on terminology.** The vocabulary of "qudits," "gates," "leakage," and "coherence" used below is the language of quantum information theory applied to the representation theory of $(\mathbb{Z}/m\mathbb{Z})^\times$. The algebraic content is that the CRT decomposition block-diagonalizes the in-factorization operators exactly, while the missing-prime operators do not commute with this decomposition and therefore mix the blocks. What the quantum-mechanical language adds is a *dynamical* interpretation: the time-evolution operator $e^{-iHt}$ propagates states through the block structure, and the measured containment ($97.1\%$ at $t = 50$, $92.6\%$ at $t = 200$) quantifies the rate at which the non-commuting terms cause inter-block transitions under unitary evolution. The "leakage" is not a static algebraic projection — it is measured as amplitude loss under $e^{-iHt}$, which is the physically meaningful quantity if $H(\alpha)$ is ever realized as the Hamiltonian of a physical system. We retain the quantum vocabulary because it provides the natural framework for these dynamical measurements, while emphasizing that the *protection* mechanism is purely algebraic (the CRT scalar identities).

#### 10.3.1 The Cosine Sieve and the $p = 3$ Boundary

Define the cosine perturbation operator $V_p = \text{diag}(\cos(2\pi r / p))$ on the coprime residues modulo $m = 2310$. The key observation (`the_connection.py`):

$$V_2 = -1 \cdot I, \qquad V_3 = -\tfrac{1}{2} \cdot I$$

on the coprime Hilbert space. These are **exact scalar operators** — not approximate. The proof: for $p = 2$, every coprime residue $r$ is odd, so $\cos(2\pi r / 2) = \cos(\pi r) = -1$. For $p = 3$, every coprime residue satisfies $r \not\equiv 0 \pmod{3}$, so $r \bmod 3 \in \{1, 2\}$ and $\cos(2\pi / 3) = \cos(4\pi / 3) = -1/2$.

For $p \geq 5$, the coprime residues sample $p - 1$ of the $p$ roots of $\cos(2\pi r / p)$ (excluding $r \equiv 0 \pmod{p}$), producing a **non-degenerate spectrum**. This is the *$p = 3$ boundary*: the cosine operator is constant on coprimes if and only if $p \leq 3$.

The hierarchy of cosine eigenvalue classes on coprime residues of $m = 2310$:

| Prime $p$ | Distinct $V_p$ values on coprimes | Eigenvalue structure | Role |
|---|---|---|---|
| 2 | 1 ($-1$) | Exact scalar $-I$ | Deleted by sieve |
| 3 | 1 ($-1/2$) | Exact scalar $-\tfrac{1}{2}I$ | Deleted by sieve |
| 5 | 4 | $\cos(2\pi k/5)$, $k = 1,2,3,4$ | Block label |
| 7 | 6 | $\cos(2\pi k/7)$, $k = 1,\ldots,6$ | Block label |
| 11 | 10 | $\cos(2\pi k/11)$, $k = 1,\ldots,10$ | Block label |
| 13 | 12 | Non-trivial within blocks | Gate operator |
| 17 | 16 | Non-trivial within blocks | Gate operator |
| 19 | 18 | Non-trivial within blocks | Gate operator |

#### 10.3.2 CRT Block Decomposition

The joint eigenspaces of $(V_5, V_7, V_{11})$ decompose the $\varphi(m) = 480$ coprime residues into $30$ blocks of dimension $16$ each. This is a direct manifestation of the Chinese Remainder Theorem: the ring $\mathbb{Z}/2310\mathbb{Z}$ factors as $\mathbb{Z}/2\mathbb{Z} \times \mathbb{Z}/3\mathbb{Z} \times \mathbb{Z}/5\mathbb{Z} \times \mathbb{Z}/7\mathbb{Z} \times \mathbb{Z}/11\mathbb{Z}$, and after removing the $p = 2, 3$ components (killed by the coprime sieve), the residue classes modulo $5 \times 7 \times 11 = 385$ provide the block labels:

The block count arises from the distinct joint eigenvalue tuples of $(V_5, V_7, V_{11})$. The cosine operator $V_p = \text{diag}(\cos(2\pi r/p))$ on coprime residues has $(p-1)/2$ distinct eigenvalues (since $\cos(2\pi k/p) = \cos(2\pi(p-k)/p)$ pairs residues). For the blocking primes: $V_5$ has $2$ distinct values, $V_7$ has $3$, and $V_{11}$ has $5$, yielding $2 \times 3 \times 5 = 30$ joint eigenspaces. Each block has dimension $\varphi(2310)/30 = 480/30 = 16$, as confirmed computationally by `crt_block_hamiltonian.py` (all 30 blocks have exactly dimension 16).

#### 10.3.3 Block-Internal Protection

Within each CRT block, the in-factorization operators $V_5, V_7, V_{11}$ act as exact scalars:

| Operator | Spread within any block (max over 30 blocks) | Interpretation |
|---|---|---|
| $V_5$ | $3.2 \times 10^{-13}$ | Numerical zero (eigensolver tolerance) — exact scalar by construction |
| $V_7$ | $3.2 \times 10^{-13}$ | Numerical zero (eigensolver tolerance) — exact scalar by construction |
| $V_{11}$ | $3.2 \times 10^{-13}$ | Numerical zero (eigensolver tolerance) — exact scalar by construction |

This protection is **exact by construction** — it follows from the definition of the CRT blocks as joint eigenspaces. It is not an approximation, not a dynamical effect, and not fragile. Any perturbation with the arithmetic structure of an in-factorization cosine operator is guaranteed to act as a scalar within each block.

The block-internal Hamiltonian $H_b = U_b^\dagger H(\alpha_c) U_b$ is a $16 \times 16$ Hermitian matrix. Noise measured as $\|H_b - \bar{\lambda}_b I\| / \|H_b\|$ (deviation from scalar):

| Block | Noise | Best gate ($V_p$, lowest noise) |
|---|---|---|
| Block 9 | $0.010$ | $V_{19}$ ($0.007$) |
| Block 20 | $0.0084$ | $V_{19}$ ($0.006$) |
| Worst block | $0.068$ | — |
| Average | $0.034$ | — |

Block 20 achieves noise immunity $0.0084$ — the Hamiltonian is $99.2\%$ scalar within this block.

#### 10.3.4 Gate Operators

The missing-prime operators $V_{13}, V_{17}, V_{19}$ are **not** scalar within CRT blocks and therefore act as non-trivial gates. Within a full 16D block (before $\sigma$-projection) the gate structure is:

| Gate | Block 9 noise | Block 20 noise | Max Rabi transfer |
|---|---|---|---|
| $V_{13}$ | $0.008$ | $0.008$ | $0.41$ |
| $V_{17}$ | $0.014$ | $0.012$ | $0.38$ |
| $V_{19}$ | $0.007$ | $0.006$ | $0.41$ |

The maximum Rabi transfer of $0.41$ (measured as the largest population transfer between any two eigenstates within a block under a gate perturbation) confirms that the CRT blocks are not the qubit basis. Because the CRT blocking operators commute exactly with $P_\sigma$ (since $\cos(2\pi(m-r)/p) = \cos(2\pi r/p)$ for $p \mid m$), each 16D block bisects into two 8D halves under $\sigma$-parity. Within the active $\sigma$-even register, the fine structure is an **8-level qudit** per qubit state — the physically relevant dimensionality. Full population inversion between any pair of qudit levels is not achieved; the gate operators distribute amplitude across multiple levels simultaneously. The true 2-level structure — the arithmetic qubit — lives at the $\sigma$-parity level (Section 11), with the CRT blocks providing Zeeman-like fine structure within each qubit level.

#### 10.3.5 Block Leakage

The CRT blocks are not invariant subspaces of $H(\alpha_c)$ itself — the Hamiltonian mixes blocks. The leakage rate, measured as the probability of finding a state initialized in Block 9 outside that block after time evolution:

| Time | Containment (probability remaining in block) |
|---|---|
| $t = 50$ | $97.1\%$ |
| $t = 100$ | $94.8\%$ |
| $t = 200$ | $92.6\%$ |

At $t = 200$, $92.6\%$ of the probability remains within the initial block. The leakage is slow but not zero, giving an estimated **coherence window of $\sim 60$ gate operations** (at the measured gate period) before the block population drops below $90\%$.

#### 10.3.6 The CRT Block Architecture (superseded by Section 11)

The CRT decomposition and the cosine sieve analysis yield the following internal architecture:

| Layer | Mechanism | What it does |
|---|---|---|
| **1. Coprime sieve** | $V_2 = -I$, $V_3 = -\tfrac{1}{2}I$ | Deletes $p = 2, 3$ from the Hilbert space — these primes cannot perturb anything |
| **2. CRT blocks** | Joint eigenspaces of $(V_5, V_7, V_{11})$ | 30 blocks × 16 states; in-factorization operators act as exact scalars within each block |
| **3. Flat band** | 473 states with $\lvert\lambda\rvert < 0.01$ | Provides a slow clock for gate timing; $\sim 60$ gates before decoherence from block leakage |
| **4. Gate operators** | $V_{13}, V_{17}, V_{19}$ (missing primes) | Non-trivial 16D rotations within each block; noise $< 0.015$ |

Each layer was identified by a separate analytical step: (1) the coprime sieve follows immediately from the definition of $\mathcal{R}$; (2) the CRT block decomposition is forced by the Chinese Remainder Theorem applied to the joint eigenspaces of in-factorization cosine operators; (3) the flat band is established by Theorem 4.1; and (4) the gate operators are identified by exclusion — the missing primes are precisely those whose cosine operators are non-scalar within blocks. The $p = 3$ boundary (Section 10.3.1) is the structural fulcrum: it separates the primes that the coprime sieve deletes from those that label the CRT blocks.

**Superseded.** The above table describes the internal fine structure of each qubit level, but it is not the full protection hierarchy. Section 11.5 gives the revised four-layer architecture in which $\sigma$-parity conservation (Theorem 11.1) provides the top-level algebraic isolation, the coprime sieve provides number-theoretic mode deletion, the flat band provides spectral averaging, and collective encoding across 240 dimensions provides dimensional dilution. The CRT blocks are the Zeeman fine structure within each $\tau$-parity subspace — spectral sublevels, not the computational basis.

---

## 11. Exact $\sigma$-Parity Conservation and the Arithmetic Qubit

The Klein-4 orbit analysis (Observation 4.6), the orbit Hamiltonian projection, and the falsification of the CRT qudit basis converged on a single question: what *is* the conserved quantum number of $H(\alpha)$? The answer is the simplest symmetry of the integers: the additive palindromic mirror $P_\sigma: r \mapsto m - r$.

**Terminology.** We use *arithmetic qubit* to describe the mathematical two-level structure of $\tau$-parity within each $\sigma$-sector — a structural property of the operator $H(\alpha)$. This is not a claim that a physical device has been built or that the system admits state preparation and measurement in the sense of quantum computing. The term denotes a proven algebraic two-level system with exact conservation laws and four layers of number-theoretic protection, not an engineered quantum register.

### 11.1 $\sigma$-Parity Conservation Theorem

**Theorem 11.1** ($\sigma$-Parity Conservation). *The active transport Hamiltonian*

$$H(\alpha) = \frac{D_{\text{sym}}}{\lambda_P} - i\alpha\,\frac{[|A|_{\text{ew}}, P_\tau]}{\lambda_P}$$

*commutes with the palindromic permutation $P_\sigma: r \mapsto m - r$ for all $\alpha \geq 0$:*

$$[H(\alpha), P_\sigma] = 0.$$

*Consequently, every eigenstate of $H(\alpha)$ has exact $\sigma$-parity $\pm 1$, and the Hilbert space decomposes as $\mathcal{H} = \mathcal{H}_+ \oplus \mathcal{H}_-$ with $\dim \mathcal{H}_\pm = \varphi(m)/2$.*

*Proof.* Two ingredients suffice.

(i) $[D_{\text{sym}}, P_\sigma] = 0$. The palindromic distance satisfies $d(m{-}r_i,\, m{-}r_j) = d(r_i, r_j)$ for all $r_i, r_j$, because $d(a,b) = \min(|a-b|,\, m - |a-b|)$ and $|(m-r_i)-(m-r_j)| = |r_i - r_j|$. Therefore $P_\sigma$ permutes the rows and columns of $D_{\text{sym}}$ without changing any entry.

(ii) $[[|A|_{\text{ew}}, P_\tau], P_\sigma] = 0$. The element-wise absolute value $|A|_{\text{ew}}$ inherits the palindromic symmetry from $D_{\text{sym}}$ (via Identity 2.5), so $[|A|_{\text{ew}}, P_\sigma] = 0$. The multiplicative inverse commutes with the palindromic mirror: $\sigma(\tau(r)) = m - r^{-1} = \tau(\sigma(r))$, because $(m - r)^{-1} \equiv m - r^{-1} \pmod{m}$ for all $r$ coprime to $m$. Therefore $[P_\tau, P_\sigma] = 0$. Since both $|A|_{\text{ew}}$ and $P_\tau$ commute with $P_\sigma$, so does their commutator $[|A|_{\text{ew}}, P_\tau]$.

Both parts of $H(\alpha)$ commute with $P_\sigma$; hence $H(\alpha)$ does. $\square$

**Verification.** Frobenius residual $\|[H(\alpha), P_\sigma]\|_F = 0.00 \times 10^{0}$ at both $m = 1155$ (odd, $\varphi = 480$) and $m = 2310$ (even, $\varphi = 480$) for all tested $\alpha$ values. The identity is algebraic, not numerical.

| $m$ | $\alpha$ | $\|[H(\alpha), P_\sigma]\|_F$ | $\sigma$-even | $\sigma$-odd | Mixed |
|---|---|---|---|---|---|
| $1155$ | $0.00$ | $0.00 \times 10^{0}$ | 240 | 240 | 0 |
| $1155$ | $0.50$ | $0.00 \times 10^{0}$ | 240 | 240 | 0 |
| $1155$ | $1.00$ | $0.00 \times 10^{0}$ | 240 | 240 | 0 |
| $1155$ | $1.24$ | $0.00 \times 10^{0}$ | 240 | 240 | 0 |
| $1155$ | $2.50$ | $0.00 \times 10^{0}$ | 240 | 240 | 0 |
| $2310$ | $0.00$ | $0.00 \times 10^{0}$ | 240 | 240 | 0 |
| $2310$ | $0.50$ | $0.00 \times 10^{0}$ | 240 | 240 | 0 |
| $2310$ | $1.00$ | $0.00 \times 10^{0}$ | 240 | 240 | 0 |
| $2310$ | $1.24$ | $0.00 \times 10^{0}$ | 240 | 240 | 0 |
| $2310$ | $2.50$ | $0.00 \times 10^{0}$ | 240 | 240 | 0 |

Every eigenstate at every $\alpha$ has exact $\sigma$-parity at both odd and even moduli. Zero mixed states.

### 11.2 Sector Structure

The $\sigma$-even sector ($\mathcal{H}_+$, dim 240) and $\sigma$-odd sector ($\mathcal{H}_-$, dim 240) have sharply distinct spectral properties at $\alpha = \alpha_c$:

| Property | $\sigma$-even ($\mathcal{H}_+$) | $\sigma$-odd ($\mathcal{H}_-$) |
|---|---|---|
| Dimension | 240 | 240 |
| Max eigenvalue ($m = 1155$) | $1.000004$ (Perron) | $0.336894$ |
| Max eigenvalue ($m = 2310$) | $1.000003$ (Perron) | $0.337587$ |
| Min eigenvalue ($m = 1155$) | $-0.744$ (anti-Perron) | $-0.742$ |
| Min eigenvalue ($m = 2310$) | $-0.745$ (anti-Perron) | $-0.742$ |
| Flat-band modes | 236 / 240 | 237 / 240 |
| Contains Perron? | Yes | No |
| Contains anti-Perron? | Yes | No |

The Perron and anti-Perron eigenstates — the two thermodynamically dominant modes — both live in the $\sigma$-even sector. The $\sigma$-odd sector is an excited manifold with a maximum eigenvalue only one-third of the Perron value. This asymmetry is a consequence of the palindromic distance matrix $D_{\text{sym}}$ having its Perron eigenvector (which converges to uniform, Observation 4.4) in the symmetric subspace of $P_\sigma$. Results are consistent across torsion structures: the odd modulus ($m = 1155$) and even modulus ($m = 2310$) yield the same sector dimensions, the same Perron/anti-Perron assignment, and match to three decimal places in all eigenvalues.

### 11.3 $\tau$-Parity Rabi Structure

Within the $\sigma$-even sector, the Klein four-group irreducible representations decompose the 240 states into $\tau$-parity subspaces:

- $(+,+)$: $\sigma$-even, $\tau$-even — 124 states (additive standing wave)
- $(+,-)$: $\sigma$-even, $\tau$-odd — 116 states (multiplicative current)

The dimension asymmetry $124 - 116 = 8$ equals the number of size-2 Klein-4 orbits. These are precisely the residues satisfying $r \equiv r^{-1} \pmod{m}$ (i.e., $r^2 \equiv 1$), which are fixed by $\tau$ and therefore contribute only to the $\tau$-even subspace, breaking the generic size-4 orbit structure. A proof that this equality is structural (not coincidental) is not given here; see Problem 20.

**Observation 11.2** ($\tau$-Parity as Dynamical Rabi Variable). While $\sigma$-parity is exactly conserved, $\tau$-parity is *dynamical* in the Rabi sense — it is not conserved, and the degree of $\tau$-mixing varies continuously with the parameter $\alpha$ (not with time; $H(\alpha)$ is a static parametric family). The $\alpha$-sweep of Klein-4 irrep weights within the $\sigma$-even sector:

| $\alpha$ | $(+,+)$ share | Cross-irrep share | Interpretation |
|---|---|---|---|
| $0.00$ | $77.7\%$ | $12.4\%$ | Additive ($\tau$-even) dominated |
| $1.24$ ($\alpha_c$) | $43.7\%$ | $50.7\%$ | **Resonant balance** — drive matches detuning |
| $2.50$ | $18.9\%$ | $78.7\%$ | Multiplicative ($\tau$-odd) dominated |

The phase transition at $\alpha_c$ is a Rabi crossing [11]: the point where the off-diagonal coupling (the $\tau$-mixing commutator $\alpha \cdot [|A|_{\text{ew}}, P_\tau]$) exactly balances the diagonal detuning ($D_{\text{sym}}$ eigenvalue splitting between $\tau$-subspaces), producing maximum $\tau$-mixing. The $50.7\%$ cross-irrep share is a Frobenius norm partition of the *operator*, not a state-vector superposition coefficient \u2014 it measures the fraction of the Hamiltonian's power that couples $\tau$-even to $\tau$-odd. Below $\alpha_c$, the system is predominantly additive geometry. Above $\alpha_c$, multiplicative topology dominates.

**$\tau$-parity smearing.** Only 1 out of 480 eigenstates (the Perron vector, with $\langle P_\tau \rangle = 1.000$) has clean $\tau$-parity. All other eigenstates have $\langle P_\tau \rangle$ ranging from $-0.75$ to $+0.53$ ($\sigma$-odd) and $-0.75$ to $+1.00$ ($\sigma$-even). The commutator $[|A|_{\text{ew}}, P_\tau]$ — the Rabi drive — generates continuous rotations between $\tau$-eigenspaces.

**Commutator strength.** Let $P_{\tau,+}$ and $P_{\tau,-}$ denote the $\tau$-even and $\tau$-odd projectors within the $\sigma$-even sector. The squared Frobenius norms of the $\tau$-off-diagonal blocks satisfy

$$\frac{\|P_{\tau,+}\, [\lvert A\rvert_{\text{ew}}, P_\tau]\, P_{\tau,-}\|_F^2}{\|P_{\tau,+}\, D_{\text{sym}}\, P_{\tau,-}\|_F^2} = 4.00$$

at both $m = 1155$ and $m = 2310$ (values: $0.166470 / 0.041617 = 4.00$ and $0.166908 / 0.041727 = 4.00$). This integer ratio suggests an algebraic origin (Problem 17).

### 11.4 The Quantum Number Hierarchy

The complete architecture is a three-level hierarchy discovered through successive falsification:

$$\underbrace{\sigma\text{-parity}}_{\text{which register}} \;\supset\; \underbrace{\tau\text{-parity}}_{\text{which qubit state}} \;\perp\; \underbrace{\text{CRT block}}_{\text{which sublevel}}$$

where $\supset$ denotes nesting (each $\sigma$-sector contains both $\tau$-subspaces) and $\perp$ denotes transversality: the CRT blocking operators $V_5, V_7, V_{11}$ do not commute with $P_\tau$ ($|A|_{\text{ew}}$ is $97.1\%$ off-diagonal in the CRT basis), so $\tau$-subspaces and CRT blocks cross independently within each $\sigma$-sector.

| Level | Quantum number | Status | Role | Dimension |
|---|---|---|---|---|
| Register | $\sigma$-parity ($\pm 1$) | **Exact** (algebraic) | Selects $\mathcal{H}_+$ or $\mathcal{H}_-$ | $2 \times 240$ |
| Qubit state | $\tau$-parity | **Dynamical** (Rabi) | Rotated by $\alpha \cdot [|A|_{\text{ew}}, P_\tau]$ | $124\,(\tau\text{-even}) \oplus 116\,(\tau\text{-odd})$ per sector |
| Sublevel | CRT block index | **Spectral** (fine structure) | Labels fine structure within each $\sigma$-sector (transverse to $\tau$-parity) | $30$ blocks $\times$ $8$D per $\sigma$-sector ($16$D pre-projection; $V_p$ commutes with $P_\sigma$ for $p \mid m$) |

The system comprises two $\sigma$-sectors: the **ground register** ($\mathcal{H}_+$, containing Perron eigenvalue $1.000$ and anti-Perron eigenvalue $-0.744$) and the **excited register** ($\mathcal{H}_-$, max eigenvalue $0.337$, no Perron, no anti-Perron). Each sector supports a $\tau$-parity qubit with the same Rabi structure, but the sectors are not spectral copies — the ground register carries the thermodynamically dominant modes. Within each sector:
- The computational basis is $\{|\tau\text{-even}\rangle,\, |\tau\text{-odd}\rangle\}$
- The Rabi drive is $\alpha \cdot [|A|_{\text{ew}}, P_\tau] / \lambda_P$
- The detuning is the $D_{\text{sym}}$ eigenvalue splitting between $\tau$-subspaces
- The resonance condition $\alpha = \alpha_c \approx \sqrt{135/88}$ marks the 50/50 crossover

### 11.5 Four-Layer Arithmetic Error Protection

The arithmetic qubit is protected by four independent mechanisms, each from a distinct mathematical source:

| Layer | Protects | Mechanism | Strength (quantified at $m = 2310$) |
|---|---|---|---|
| **1. $\sigma$-conservation** | Register selection | $[H(\alpha), P_\sigma] = 0$ — algebraic identity | Exact for palindromic-symmetric noise; linear breaking otherwise |
| **2. Coprime sieve** | Fourier mode space | Number-theoretic deletion of all modes divisible by $p \leq p_k$ | $Q \sim 10^7$ (transition width $\sim 10^{-7}$) |
| **3. Flat band** | Phase coherence | 473/480 eigenvalues degenerate near zero | Robust to $\varepsilon \lesssim 0.1$ ($\sim 6\%$ of $\|H\|_F$) |
| **4. Collective encoding** | Single-site errors | Qubit state spread across 240 dimensions | Perron overlap $99.8\%$ at $\varepsilon = \|H\|_F$; flat band $464/473$ at $\varepsilon = 10$ |

**Layer 1: Algebraic isolation (conditional on palindromic symmetry).** The $\sigma$-parity conservation is an algebraic identity (Theorem 11.1). Under any perturbation $E$ that respects the palindromic involution ($P_\sigma E P_\sigma = E$), the commutator $[H(\alpha) + \varepsilon E,\, P_\sigma] = 0$ remains identically zero at all $\varepsilon$ — verified computationally at $m = 2310$ through $\varepsilon = 10$ ($6.5\times \|H\|_F$). Generic perturbations that break the $r \mapsto m-r$ symmetry produce $\|[H', P_\sigma]\|_F \propto \varepsilon$ (linear degradation), destroying eigenvector $\sigma$-purity rapidly: at $\varepsilon = 0.001$ ($0.07\%$ of $\|H\|_F$), only $35/480$ eigenstates retain $|\langle P_\sigma \rangle| > 0.99$. The protection is therefore *exact and infinite* for palindromic-symmetric noise, but *conditional* — generic noise breaks it at first order.

**Layer 2: Number-theoretic deletion.** The Hilbert space is $(\mathbb{Z}/m\mathbb{Z})^\times$ — only the $\varphi(m)$ residues coprime to $m$. For $m = 30030$, this deletes $77\%$ of all residues. Every Fourier mode divisible by $2, 3, 5, 7, 11$, or $13$ is absent from the state space. A perturbation would need to excite one of these deleted modes — they are not suppressed, they are *structurally absent*.

**Layer 3: Spectral averaging (robust to $\sim 10\%$ of $\|H\|_F$).** Within the $\sigma$-even sector, 236 of 240 eigenvalues sit in the degenerate flat band near zero. A random perturbation must coherently shift these $\sim 236$ near-degenerate eigenvalues to cause phase decoherence. Quantified at $m = 2310$: the flat band ($473/480$ states at $|\lambda| < 0.01$) survives intact through $\varepsilon = 0.03$ ($2\%$ of $\|H\|_F$), degrades to $463/480$ at $\varepsilon = 0.1$ ($6.5\%$), and collapses sharply at $\varepsilon = 0.3$ ($20\%$; only $216/480$ survive). Single-site perturbations leave the flat band nearly untouched ($464/480$ at $\varepsilon = 10$) — the collective degeneracy is robust to local errors.

**Layer 4: Dimensional dilution (quantified).** The qubit is not a single site — it is a collective mode across 240 dimensions. Single-site perturbations at $m=2310$ shift the collective Perron mode by $\Delta\lambda \approx \varepsilon/\varphi$: at $\varepsilon = 0.1$, the shift is $+0.00023$ with $99.997\%$ overlap to the clean Perron state, matching the first-order prediction $\varepsilon/\varphi = 0.00021$. Beyond $\varepsilon \sim 0.3$, the localized defect state hybridizes with the collective mode; at $\varepsilon \geq 3$, the defect separates and the collective mode recovers (overlap $> 99.5\%$). The flat band remains intact ($464/473$ states) even at $\varepsilon = 10$. The Perron state overlap $|\langle v_0 | v_\varepsilon \rangle|^2$ remains $99.8\%$ at $\varepsilon = 1.0$ and $98.1\%$ at $\varepsilon = 3.0$ ($2\times \|H\|_F$). The max/min amplitude ratio grows from $1.007$ (clean) to $1.48$ at $\varepsilon = 1.0$, confirming gradual degradation rather than catastrophic failure.

**Noise model specification.** The four-layer protection hierarchy is effective against (a) any perturbation preserving palindromic symmetry (Layer 1, exact), (b) perturbations with amplitude $\varepsilon \lesssim 0.1\|H\|_F$ (Layer 3, flat-band threshold), and (c) local single-site errors which are diluted by the collective encoding factor $\sim \varphi/2$ (Layer 4). Generic symmetry-breaking noise at amplitude $\varepsilon$ produces $\sigma$-parity violation $\propto \varepsilon$. No fault-tolerance threshold theorem is claimed — the protection is *arithmetic*, not engineered.

### 11.6 Discovery Narrative

The path to the arithmetic qubit followed a successive-falsification approach — each wrong prediction eliminated a candidate basis until the correct one emerged:

| Prediction | Test | Result |
|---|---|---|
| CRT blocks are the qudit basis | Eigenstates respect blocks? | **Falsified** — $7.4\%$ leakage at $t = 200$ |
| Klein-4 orbits are the basis | OPR $> 90\%$? | **Falsified** — $41\%$ max, $0/473$ modes trapped |
| Orbits respect CRT blocks | Orbits in single block? | **Falsified** — $0/124$ at odd $m$ |
| The system is frustrated | Orbit Hamiltonian diagonal? | **Falsified** — $93.7\%$ off-diagonal |
| $\sigma$ exact, $\tau$ Rabi | $[H(\alpha), P_\sigma] = 0$? | **Confirmed** — $0.00 \times 10^{0}$ at all $\alpha$ |

The orbit Hamiltonian projection (capturing $44.1\%$ of the norm, with $O(2)$ as a $99.2\%$ resonance and $93.7\%$ CRT off-diagonal power) initially suggested fundamental frustration. The resolution came from recognizing that the correct block-diagonalization is not by CRT labels or Klein-4 orbits, but by the exact $\sigma$-symmetry — the simplest involution of the system, hiding in plain sight.

---

## 12. Open Problems

### Spectral Constants and Scaling Laws

**Problem 1** (Flat-band constant). Prove analytically that $\text{Var}(\lambda) \cdot \varphi(m) \to V_\infty$ and determine the value of $V_\infty \approx 1.780$.

**Problem 2** (Anti-Perron eigenvalue). Identify the limiting anti-Perron eigenvalue $\lambda_-$ as the root of a polynomial equation, or prove it is transcendental. At $m = 30030$ ($\alpha = \alpha_c$), $\lambda_- = -0.7441353$ with gap $\Delta\lambda = \lambda_P - \lambda_- = 1.7441 \pm 3.7 \times 10^{-4}$, stable across 203 $\alpha$-values spanning the full crash window $[\sqrt{3/2},\, \alpha_c]$. The value depends on $\alpha$; at $\alpha = \sqrt{2/3}$, $\lambda_- \approx -0.591$.

**Problem 3** (Critical shear closed form). Confirm or refute Conjecture 6.1: $\alpha_c = \sqrt{135/88}$. The $m = 30030$ micro-sweep gives $+0.015\%$ error — convergence continues but is sharply decelerating (local exponent $\beta \approx 2.0 \to 0.24$ across successive primorial pairs), ruling out a stable power law. Two scenarios remain open: (a) $\alpha_c = \sqrt{135/88}$ exactly, with logarithmic corrections slowing convergence; or (b) the true limit is a nearby transcendental number, and $\sqrt{135/88}$ is a Padé-type approximant to the local geometry. Distinguishing these requires $m = 510510$ ($\varphi = 92160$). If confirmed algebraic, explain why the critical ratio involves $3^3 \cdot 5 / (2^3 \cdot 11)$ — why these specific prime powers?

**Problem 4** (IPR growth). Determine whether the IPR ratio grows as $\Theta(\ln \varphi)$ or follows a different scaling law. A proof of logarithmic growth would establish $H(\alpha)$ as exhibiting *critical delocalization* (neither Anderson-localized [6] nor fully delocalized).

### Phase Transition Universality and Torsion

**Problem 5** (Extension beyond primorials). The non-primorial control $m = 385 = 5 \cdot 7 \cdot 11$ ($\varphi = 240$) exhibits the phase transition (Section 6.6), confirming that it is not primorial-specific. Characterize the class of squarefree moduli that admit the transition. Does every squarefree modulus with $\varphi(m) > 48$ exhibit it?

**Problem 6** (Torsion width scaling). Sylow torsion smears the transition without shifting $\alpha_c$ (Section 8.5). Does the width ratio (primorial vs non-primorial) scale with the number of missing small primes? The torsion–localization gap is negligible at $\varphi = 480$ ($2.46\times$ vs $2.47\times$) but reaches $7.9\%$ at $\varphi = 5760$ ($2.88\times$ vs $2.67\times$). Does this gap continue to widen, and does it saturate?

**Problem 7** (Washing failure mechanism). The mod-4 partition non-multiplicativity (Section 8.2) means torsion defects survive primorial expansion. Does the defect magnitude $\Delta v_q$ grow, stabilize, or oscillate along the primorial hierarchy? The computed values $\Delta v_5(2310) = 4$ and $\Delta v_3(510510) = \Delta v_5(510510) = 16$ suggest growth, but the data is sparse.

**Problem 8** ($m = 510510$ spectral computation). Compute $H(\alpha)$ at $m = 510510 = 2 \cdot 3 \cdot 5 \cdot 7 \cdot 11 \cdot 13 \cdot 17$ ($\varphi = 92160$, matrices $92160 \times 92160$). This modulus has chirality balance $I = 0$ and non-zero torsion defects ($\Delta v_3 = \Delta v_5 = 16$). The computation would require approximately $64$ GB of RAM and distributed eigensolvers.

### Perturbation Response and CRT Block Architecture

**Problem 9** (Screening ratio scaling). At $m = 2310$ ($\varphi = 480$), the IN/OUT selective screening ratio at $\alpha_c$ is $1.36$ (Section 10.1). Does this ratio grow with $m$? At $m = 30030$ ($\varphi = 5760$), does IN/OUT reach $2\times$ or higher? Conversely, at $m = 210$ ($\varphi = 48$), is the effect weaker or absent? A monotone increase would establish selective screening as an asymptotic property of the primorial hierarchy. The perfect screening of $p = 2, 3$ (Section 10.1.3) is a sieve identity and therefore holds at all $m$; the open question is whether the *partial* screening of larger in-factorization primes ($p = 5, 7, 11, \ldots$) strengthens with $\varphi$.

**Problem 10** (CRT qudit gate universality). Within each $\sigma$-sector, the CRT blocks are 8-dimensional (Section 10.3.4). Do $V_{13}, V_{17}, V_{19}$ generate the full $SU(8)$ (or a useful subgroup) on each $\sigma$-projected block? If not, what subgroup do they generate? The measured gate noise values ($0.007$–$0.014$) and max Rabi transfer ($0.41$) on the best blocks (Block 9, Block 20) suggest non-trivial but incomplete rotations.

**Problem 11** (Block leakage suppression). Can the Hamiltonian be modified to make CRT blocks exactly invariant? One approach: project $H(\alpha_c)$ into the block-diagonal subspace and use the projected Hamiltonian. What is the energy cost of the off-diagonal (inter-block) terms?

**Problem 12** (Scaling to $m = 30030$). At $m = 30030 = 2 \cdot 3 \cdot 5 \cdot 7 \cdot 11 \cdot 13$, the prime $p = 13$ moves from "gate operator" to "block label." The CRT blocks become joint eigenspaces of $(V_5, V_7, V_{11}, V_{13})$, of higher count but lower dimension. Does the protection strengthen or does the reduced block dimension limit gate expressivity? Does the architecture persist, and does the screening ratio increase?

### Physical Realization

**Problem 13** (Landau-Zener dynamics and physical realization). The narrow transition width $\delta \sim 10^{-7}$ at $m = 30030$ implies that if $\alpha$ is a time-dependent parameter, even slow sweeps produce complete non-adiabatic population transfer via the Landau-Zener mechanism (Section 9.3). Does the anti-Perron asymmetry $\mathcal{A}_{\text{anti}} \approx -0.003$ determine the spatial direction of energy release in such a transition? Establishing this requires identifying a physical system whose transport operator is $H(\alpha)$. The gap stability ($\Delta\lambda = 1.7441 \pm 3.7 \times 10^{-4}$ over 203 $\alpha$-values) guarantees a fixed energy scale for the transition.

**Problem 14** (Physical implementation). The wiring topology of $D_{\text{sym}}$ on coprime residues is global — every coprime residue couples to every other. Is there a local approximation that preserves the CRT block structure?

### Chirality Identity Chain and Isospectral Non-Conjugacy

**Problem 15** (Isospectral non-conjugacy and the folding miracle). The element-wise operation $A \to |A|$ (where $|A|_{ij} = |A_{ij}|$) is non-linear: it maps the antisymmetric directed distance matrix to its symmetric Hadamard absolute value. Despite this non-linearity, $|A|/m$ has *exactly identical eigenvalues* for matched pairs (same $\varphi$, different 2-torsion) — verified to $< 10^{-13}$ at $\varphi = 8, 48, 480$ (Section 2.6). Yet the eigenvectors differ: the Frobenius non-conjugacy error under the coprime bijection $\sigma$ grows with $\varphi$, from $1.5$ to $9.8$, and the counterfactual $P_\tau$ swap (Section 4.5) confirms that these eigenvector differences alone determine the torsion–localization gap. Two sub-questions: (a) Prove algebraically why the folding $A \to |A|$ preserves the singular value spectrum across matched pairs. Is this a consequence of a hidden symmetry in the cotangent eigenphase structure of $A$, or does it require a deeper invariant? (b) Characterize the non-conjugacy: is it detectable by a finite set of matrix moments $\operatorname{tr}((|A|/m)^k)$, or does it require the full eigenvector structure?

**Problem 16** (Kurtosis scaling). At $\varphi = 480$, the eigenvector kurtosis ratio between odd and even moduli is $117 / 6.9 \approx 17\times$ (Section 4.5, Experiment 3). Does this ratio grow with $\varphi$? If so, at what rate? The coprime gap parity constraint (odd $m$ allows odd gaps; even $m$ forces even gaps) is the identified mechanism; a proof that odd-gap moduli produce heavier-tailed eigenvector distributions would explain the torsion crossover analytically.

### $\sigma$-Parity and the Arithmetic Qubit

**Problem 17** (Commutator-to-$D_{\text{sym}}$ ratio). The ratio of commutator power to $D_{\text{sym}}$ power in the $\tau$-off-diagonal channels is exactly $4.00$ at $m = 1155$ (both moduli of the coupling matrices). Is this ratio universal across all moduli? If so, derive it algebraically from the structure of $|A|_{\text{ew}}$ and $P_\tau$. The integer value suggests a combinatorial identity connecting the palindromic distance to the multiplicative commutator.

**Problem 18** ($\sigma$-odd sector structure). The $\sigma$-odd sector has maximum eigenvalue $0.337$ at $m = 1155$ ($\alpha = \alpha_c$), compared to $1.000$ in the $\sigma$-even sector. Characterize the $\sigma$-odd spectrum: does the maximum eigenvalue converge to a universal constant as $\varphi \to \infty$? Why does the $\sigma$-odd sector contain neither the Perron nor the anti-Perron eigenstate?

**Problem 19** ($\sigma$-parity at larger moduli). Theorem 11.1 ($[H(\alpha), P_\sigma] = 0$) and the $\tau$-Rabi crossing are verified at $m = 1155$ and $m = 2310$ (Section 11.1). Extend to $m = 30030$ and $m = 510510$: does the $(+,+)$ share and cross-irrep share converge, and does the 50/50 crossover at $\alpha_c$ persist? Universality of the arithmetic qubit structure requires confirmation at moduli where CRT blocks are higher-count and lower-dimension.

**Problem 20** (Klein-4 irrep dimension asymptotics). At $\varphi = 480$, the dimension split is $124 + 116$ ($\sigma$-even, $\tau$-even vs $\tau$-odd), with the gap of $8$ equal to the number of size-2 Klein-4 orbits. Does this gap grow, stabilize, or vanish relative to $\varphi$ as $m \to \infty$? If it vanishes, the detuning between $\tau$-subspaces approaches zero, and the Rabi oscillation becomes symmetric. **Hint:** By the Chinese Remainder Theorem, the number of solutions to $r^2 \equiv 1 \pmod{m}$ for $m = 2 \cdot 3 \cdots p_k$ is exactly $2^{k-1}$, so the gap is $2^{k-2}$. Since $\varphi(m)/2$ grows super-exponentially while $2^{k-2}$ grows only exponentially in $k$, the detuning asymmetry provably vanishes: the Rabi oscillation becomes symmetric in the primorial limit.

---

## 13. Appendix A: Computational Details

### 13.1 Implementation

All computations were performed in Python 3.13 using NumPy with scipy-openblas for linear algebra. The initial scripts (`chiral_quench_test.py`, `chiral_phase_transition.py`) are single-file, single-threaded. Later scripts (`micro_sweep.py`, `overnight_batch.py`) use `concurrent.futures.ProcessPoolExecutor` with 8 workers and 4 OMP threads each for parallel eigendecomposition on a 32-core AMD consumer PC.

### 13.2 Matrix Construction

- **$D_{\text{sym}}$** is constructed via broadcasting: $D_{ij} = \min(|r_i - r_j|, m - |r_i - r_j|)$ using `np.minimum(diff, m - diff)` where `diff = np.abs(r[:, None] - r[None, :])`.
- **$P_\tau$** is represented as an index permutation array (not the full matrix) computed via `pow(r, -1, m)`.
- **$[D, P_\tau]$** is computed without full matrix multiplication: `D[:, perm] - D[perm, :]` permutes columns and rows directly ($\mathcal{O}(n^2)$, not $\mathcal{O}(n^3)$).

### 13.3 Eigendecomposition

The Hermitian matrix $H(\alpha) = A + i\alpha B$ (real symmetric $A$, real skew-symmetric $B$) is constructed as a complex Hermitian matrix and decomposed using `numpy.linalg.eigvalsh` (eigenvalues only) or `numpy.linalg.eigh` (eigenvalues and eigenvectors, for IPR computation). The LAPACK `zheevd` routine provides $\mathcal{O}(n^3)$ performance.

### 13.4 Timing

Single-threaded timing (original scripts, numpy.linalg.eigvalsh backed by LAPACK zheevd):

| $m$ | $\varphi(m)$ | $D_{\text{sym}}$ construction | Eigendecomposition (per $\alpha$) | Total (62 $\alpha$-values) |
|---|---|---|---|---|
| 30 | 8 | $< 0.01$s | $< 0.01$s | $< 0.1$s |
| 210 | 48 | $< 0.01$s | $< 0.01$s | $< 0.1$s |
| 2310 | 480 | $< 0.01$s | 0.13s | 8.0s |
| 30030 | 5760 | 0.2s | 16.3s | 1008s |

Parallel timing (overnight_batch.py, 8 workers × 4 OMP threads, scipy-openblas):

| $m$ | $\varphi(m)$ | Eigendecomposition (per $\alpha$) | micro_sweep 2003 points |
|---|---|---|---|
| 385 | 240 | $< 0.01$s | ~1 min |
| 2310 | 480 | 0.05s | ~5 min |
| 30030 | 5760 | 14.7s | ~8.2 h (est.) |

The fine-resolution sweep (`chiral_phase_transition.py`, single-threaded) required 222 eigendecompositions per modulus. At $m = 30030$, each $5760 \times 5760$ Hermitian eigendecomposition takes $\sim 16.2$s, yielding **3593 seconds (60 minutes)** for the sweep. The micro-sweep (`micro_sweep.py`, parallel) resolves 2003 points in the crash window; at $m = 2310$ this completes in $\sim 5$ minutes.

### 13.5 Numerical Stability

- The Perron eigenvalue is computed by power iteration (500 iterations, tolerance $10^{-14}$), providing the normalization $\lambda_P$ to machine precision before the full eigendecomposition.
- Skew-symmetry of $[D, P_\tau]$ is verified to be exact (residual $0.00$) at each modulus, confirming no floating-point contamination of the Hermitian structure.
- The IPR computation at $m = 30030$ requires a full eigenvector solve ($5760 \times 5760$ complex), taking 119 seconds.

---

## 14. Script Inventory

The scripts below are presented in the order the investigation unfolded — from initial discovery through falsification to the corrected architecture. Each script is self-contained, reads no external data files, and produces all results from first principles.

### 14.1 Spectral Discovery (Sections 2–7)

These scripts established that the active transport operator $H(\alpha)$ undergoes a sharp spectral phase transition and that the flat band condenses $99.7\%$ of all eigenvalues.

| Script | What it does and what it found | Runtime |
|---|---|---|
| `chiral_quench_test.py` | Constructs $H(\alpha)$ from scratch for four primorials ($m = 30$ to $30030$), sweeps 62 values of the coupling $\alpha$, and computes the full eigenspectrum at each. This is the script that first revealed the flat-band condensation and the spectral width convergence to $W = 1.5908$. | ~19 min (full) |
| `chiral_phase_transition.py` | Zooms into the phase transition with 222 $\alpha$-values per modulus. Hunts for the exact crash point where the gap position jumps, measures the transition width ($10^{-7}$ at $m = 30030$), and computes IPR profiles across the crash. Produced the derivative $-4{,}989{,}133$ that quantifies the transition sharpness. | ~62 min (full) |
| `micro_sweep.py` | Ultra-high-resolution probe: 2003 points in the crash window. Resolves the transition profile at a resolution where individual eigenvalue trajectories become visible. Uses `ProcessPoolExecutor` for parallel eigendecomposition. | ~5 min ($m = 2310$) |
| `overnight_batch.py` | Parallel batch runner (8 workers, 4 OMP threads each) that executed five independent jobs in a single overnight session: the $m = 385$ non-primorial test, full post-transition eigenvector extraction, the $m = 30030$ micro-sweep, the $m = 15015$ comparison, and large-scale eigenvector solves. | ~7 h |
| `cache_results.py` | Caching utilities — saves eigenvalue arrays to JSON/NPZ so that downstream scripts can reload results without recomputing. | — |
| `data_run_ipr_gap.py` | Focused data-collection script for paper revision. Task 1: computes IPR at $m = 385$ (non-primorial, $\varphi = 240$) at four $\alpha$ values including $\alpha_c(385)$. Task 2: dense gap-stability sweep at $m = 30030$ with 203 $\alpha$-values across $[\sqrt{3/2},\, \alpha_c]$, confirming $\Delta\lambda = 1.7441 \pm 3.7 \times 10^{-4}$. Uses `ProcessPoolExecutor` for parallel eigendecomposition. | < 1 s (IPR, from cached eigenvectors); ~68 min (gap, parallel eigendecomposition) |
| `chirality_analysis.py` | Proves Identity 2.5 ($D_{\text{sym}} = m/2 \cdot (J-I) - \lvert A\rvert$) and Identity 2.6 ($[D_{\text{sym}}, P_\tau] = -[\lvert A\rvert, P_\tau]$) computationally at five moduli. Verifies $\lvert A\rvert/m$ exact isospectrality for matched pairs. Classifies $A$ entries (integer at even $m$, half-integer at odd $m$). Tests V₂ chirality projection of all cotangent modes (zero at odd $m$). | < 1 min |
| `chirality_analysis_part2.py` | Measures $P_\tau$ parity mixing (0% at even $m$, 50–58% at odd $m$), commutator norms (ratio $\approx 1.000$ for matched pairs), exact $P_\tau$ conjugacy under $\sigma$, and imaginary-vs-real energy decomposition. | < 1 min |
| `chirality_analysis_part3.py` | The counterfactual $P_\tau$ swap experiment (Section 4.5): proves $D_{\text{sym}}$ alone determines IPR, $P_\tau$ is spectrally inert, $\lvert A\rvert/m$ is isospectral but not conjugate, and documents the kurtosis explosion (117 vs 6.9 at $\varphi = 480$). The definitive experiment for the torsion crossover mechanism. | < 2 min |

### 14.2 Torsion and Eigenvector Structure (Sections 8–9)

With the spectral landscape mapped, these scripts probed the arithmetic structure hidden in the eigenvalues and eigenvectors — torsion defects, ghost matrices, and the sign-flip that revealed the role of 2-torsion.

| Script | What it does and what it found | Runtime |
|---|---|---|
| `sylow_defect_predictor.py` | Predicts $p$-adic valuation defects $\Delta v_q$ using the Bragg condition $p \equiv 1 \pmod{q}$. A purely arithmetic computation that explains why the eigenvalue spacing at different primes deviates from the flat-band prediction. | < 1s |
| `verify_sylow_5adic.py` | Verifies the 5-adic torsion defect $\Delta v_5 = 4$ at $m = 2310$ by five independent methods (direct eigenvalue counting, CRT residue analysis, Bragg condition, explicit Sylow orbit construction, and power-iteration convergence rate). | < 1 min |
| `verify_510510_zero_defect.py` | Tests the "washing conjecture" — whether torsion defects vanish at $m = 510510$ where all small primes are included. Constructs a $92160 \times 92160$ matrix. Result: FALSIFIED — defects persist. This is the largest computation in the paper. | ~30 min |
| `additive_ghost_extractor.py` | Builds the ghost matrix $D_{\text{actual}} - D_{\text{ideal}}$ that captures the additive imprint of prime gaps on the distance matrix. Tests whether the ghost's asymmetry can predict chiral properties of the eigenvectors. | ~1 min |
| `actual_eigenvector_test.py` | Tests two early conjectures about the Perron eigenvector: (1) that it carries measurable chiral asymmetry, and (2) that it contains a 13-periodic Fourier component. Both FALSIFIED — the flat-band condensation forces the Perron eigenvector to be uniform, making both predictions theorem-level impossible. | < 1 min |
| `nozzle_test.py` | Compares the anti-Perron eigenvector at even modulus ($m = 30030$) vs. odd modulus ($m = 15015$), holding $\varphi = 5760$ constant. Discovered the chiral asymmetry sign flip: $\mathcal{A}_{\text{anti}} = -0.003$ at even $m$ reverses to $+0.002$ at odd $m$. This was the first evidence that 2-torsion controls the exhaust direction (Section 9.4). | < 1 min |
| `partition_4class.py` | Resolves whether the sign flip is genuine or a partition artifact. Decomposes the anti-Perron eigenvector mass into all four mod-4 residue classes at both moduli. Confirms the sign flip is real (sector-normalized asymmetry amplifies it from $+0.002$ to $+0.005$) and discovers the Klein four-group mirror symmetry: classes $\{0, 3\}$ are depleted identically, classes $\{1, 2\}$ enriched identically (Section 9.4). | < 1 min |

### 14.3 Perturbation Robustness (Section 10.1)

The flat band looked like a protected subspace. These scripts asked: protected against what?

| Script | What it does and what it found | Runtime |
|---|---|---|
| `quantum_prime_dynamics.py` | Computes the time-averaged Loschmidt echo $\langle L \rangle_\infty$ — the probability that a localized state returns to itself — at three values of $\alpha$ ($0$, $\alpha_c$, and $2.5$). Establishes that flat-band caging is universal across $\alpha$, not specific to the critical point. | < 2 min ($m = 2310$) |
| `perturbation_robustness.py` | Applies random Hermitian noise at five amplitudes ($\varepsilon = 0.01$ to $1.0$) and measures survival. Establishes fragility: at $\varepsilon = 0.03$ (less than $3\%$ of operator norm), two-thirds of the caging is already destroyed. Rules out topological protection. | < 5 min ($m = 2310$) |
| `prime_resonance_perturbation.py` | Replaces random noise with arithmetically structured perturbations — cosine potentials $V_p = \text{diag}(\cos(2\pi r / p))$ for primes both in and out of the factorization of $m$. Discovers selective screening: in-factorization primes are screened $10$–$20\times$ better than random noise, and the IN/OUT ratio reaches $1.36\times$ at $\alpha_c$. Also tests hopping perturbations (off-diagonal), which show no IN/OUT distinction. | < 5 min ($m = 2310$) |

### 14.4 The Qubit Hypothesis and Its Correction (Section 10.2)

The screening results suggested the flat band might support a noise-immune quantum gate. These scripts tested that hypothesis and identified the correct mechanism.

| Script | What it does and what it found | Runtime |
|---|---|---|
| `primorial_logic_gate.py` | First attempt at a quantum gate: applies $V_{13}$ to the flat-band-projected identity residue $\lvert r = 1\rangle$ and tracks the fidelity oscillation. Result: swing of only $4.6\%$ — not a useful rotation. Diagnoses the flat-band obstruction: the CRT fixed-point arithmetic of $r = 1$ pins the state to a single eigenvalue plateau. | < 10 min ($m = 2310$) |
| `flatband_internal_structure.py` | Decomposes $V_{13}$ within the flat-band subspace to understand the obstruction. Discovers the eigenvalue plateau at $\lambda \approx 0.885$. Constructs an optimal Rabi state from a flat-band pair and a cavity state from Perron + anti-Perron superposition. The cavity state showed period $3.51$ and swing $1.0$ — initially interpreted as a perfect Rabi oscillation. | < 10 min ($m = 2310$) |
| `primorial_universality_test.py` | Pre-registered falsification test with explicit kill conditions. Tests whether $V_5, V_7, V_{11}, V_{13}, V_{17}, V_{19}$ produce different oscillation periods (as required for a driven Rabi gate). Result: **all primes produce the identical period** $T \approx 3.60 = 2\pi/\Delta\lambda$. Kill condition K3 triggered — the oscillation is free precession of the energy gap, not a driven rotation. The "primorial qubit" (Rabi mechanism) is dead. | < 10 min ($m = 2310$) |
| `corrected_rabi.py` | Confirmation test: measures the maximum population transfer between the two levels of the supposed Perron–anti-Perron qubit. Result: $0.41$, not $1.0$ — the CRT blocks form a 16-dimensional qudit, not the qubit basis (the true qubit is $\sigma$-parity, Section 11). | < 5 min ($m = 2310$) |
| `reconstructed_rabi.py` | Independent reconstruction of the Rabi analysis from scratch, starting from the raw Hamiltonian. Confirms the free precession interpretation with no reliance on cached intermediate results. | < 5 min ($m = 2310$) |
| `true_qubit.py` | Searches for any 2-level subspace within the CRT blocks that supports full population inversion. Result: NOT FOUND within CRT blocks — the 16-dimensional block structure distributes amplitude across too many levels. This failure redirected the search from CRT blocks to $\sigma$-parity, where the true qubit lives (Section 11). | < 5 min ($m = 2310$) |

### 14.5 CRT Block Architecture (Section 10.3)

The falsification of the primorial qubit (Perron–anti-Perron Rabi mechanism) revealed an intermediate layer: the Chinese Remainder Theorem organizes the Hilbert space into arithmetically protected blocks. These scripts characterize that architecture.

| Script | What it does and what it found | Runtime |
|---|---|---|
| `the_connection.py` | The key script. Proves that $V_2 = -I$ and $V_3 = -\frac{1}{2}I$ on the coprime Hilbert space (exact identities, not approximations). Establishes the $p = 3$ boundary: the cosine operator is constant on coprimes if and only if $p \leq 3$. This is the mechanism behind perfect screening of $p = 2, 3$ and the starting point of the CRT block decomposition. | < 5 min ($m = 2310$) |
| `crt_block_hamiltonian.py` | Constructs the joint eigenspaces of $(V_5, V_7, V_{11})$, yielding 30 blocks of dimension 16. Projects $H(\alpha_c)$ into each block and measures the within-block spread of $V_5$: $3.2 \times 10^{-13}$ (numerical zero, eigensolver tolerance). The in-factorization operators are exact scalars within each CRT block — this is the algebraic identity behind the "noise immunity." | < 10 min ($m = 2310$) |
| `block9_deep_dive.py` | Systematic survey of all 30 CRT blocks: measures block-internal noise, tests each missing-prime operator ($V_{13}, V_{17}, V_{19}$) as a gate, and records the maximum Rabi transfer within each block. Identifies Block 20 as the best-protected block (noise $= 0.0084$, i.e., $99.2\%$ scalar) and measures the global maximum Rabi transfer of $0.41$. | < 15 min ($m = 2310$) |
| `flatband_qubit_search.py` | Intermediate script on the path from screening to CRT. Discovers the three-tier screening hierarchy (core primes → in-factorization → missing primes) that motivated the search for block structure. | < 10 min ($m = 2310$) |
| `find_connection.py` | Exploratory search for the algebraic link between the screening measurements (Section 10.1) and the CRT decomposition (Section 10.3). Tested several candidate connections before `the_connection.py` found the correct one. | < 5 min ($m = 2310$) |
| `invariant_subspace.py` | Tests whether CRT blocks are exactly invariant under $H(\alpha_c)$. They are not — the Hamiltonian mixes blocks, producing slow leakage. Measures containment: $97.1\%$ at $t = 50$, $92.6\%$ at $t = 200$, giving a coherence window of $\sim 60$ gate operations (Section 10.3.5). | < 5 min ($m = 2310$) |
| `flatband_cavity_resist.py` | Tests how long a state initialized in the flat-band cavity resists leakage into the dispersive modes. Intermediate result that informed the leakage analysis. | < 5 min ($m = 2310$) |

### 14.6 CRT-Chirality Bridge (Sections 2.5, 4.5)

With the chirality identity chain established (Section 2.5) and the CRT block architecture characterized (Section 10.3), these scripts tested whether the two structures are connected — specifically, whether the element-wise absolute value $|A|$ respects CRT block boundaries.

| Script | What it does and what it found | Runtime |
|---|---|---|
| `crt_chirality_bridge.py` | Tests five predictions connecting the chirality identity chain to CRT qudit architecture: (P1) $\|A\|$ block-diagonalizes in CRT basis, (P2) flat band = $N{-}1$ passive blocks, (P3) odd $m$ promotes $V_2$ to block label ($60 \times 8$D), (P4) commutator $[\|A\|, P_\tau]$ is a sparse routing bus, (P5) kurtosis hoards into active block. **Critical bug:** used matrix absolute value $\|A\| = (A^\top A)^{1/2}$ instead of element-wise $\|A\|_{ij} = \|A_{ij}\|$, producing 110% Identity 2.5 verification error. Results P1, P2, P4 invalidated. P3 confirmed ($60 \times 8$D at odd $m$). P5 falsified (kurtosis spread across all blocks, BPR $= 0.926$). Superseded by `flatband_hotspot_hunter.py`. | < 2 min |
| `flatband_hotspot_hunter.py` | Corrected CRT-chirality bridge using element-wise $\|A\|_{ij} = \|A_{ij}\|$. **Path A:** Verifies Identity 2.5 at numerical zero ($0.00$) for both $m = 2310$ and $m = 1155$; confirms Identity 2.6 at numerical zero. Tests $\|A\|_{\text{ew}}$ block-diagonality in CRT basis: **97.1% off-diagonal** at $m = 2310$, **98.7%** at $m = 1155$ — additive and multiplicative structures are genuinely orthogonal. **Path B (Flat-band autopsy):** Extracts all 473 flat-band eigenvectors, ranks by IPR, measures Block Participation Ratio. Kurtosis explosion lives entirely in the flat band (115.1 at odd $m$ vs 4.1 at even $m$). Most localized mode at $m = 1155$ spikes at the mirror pair $r = 2$ and $r = m{-}2$ (10.66% each) — boundary localization forced by $\tau$-symmetry and coprime gap parity. BPR of extreme-IPR modes $= 0.63{-}0.67$ of $N_{\text{blocks}}$ (intermediate). $\|A\|_{\text{ew}}$ intra-block spectral CV $= 3{-}4\%$ (not isospectral across blocks). | < 30 s |
| `klein4_orbit_opr.py` | Partitions coprime residues into Klein four-group orbits under $\langle \sigma, \tau \rangle$ ($\sigma(r) = m{-}r$, $\tau(r) = r^{-1}$) and measures Orbit Participation Ratio (OPR) for flat-band eigenvectors. **Odd $m = 1155$:** 124 orbits (116 size-4, 8 size-2). $O(2) = \{2, 577, 578, 1153\}$ exists; top mode ($\text{IPR} \times \varphi = 21.8$) concentrates 41.3% mass in $O(2)$, OPR $= 5.53/124$ orbits. 4/10 top-localized modes dominated by $O(2)$. Perron eigenvector fully delocalized (OPR $= 122/124$). $\tau$-parity smeared: $\langle v \| P_\tau \| v \rangle \in [0.005, 0.49]$ — confirms $[\|A\|, P_\tau]$ rotates between $\tau$-eigenspaces. **Even $m = 2310$:** $O(2)$ deleted ($\gcd(2,2310) \neq 1$). Top mode OPR $= 26.1/124$ (5$\times$ more spread), peak orbit mass 10.9%. $\tau$-parity closer to $\pm 1$ (up to 0.82). **Conclusion:** flat band is a coupled orbit network with $O(2)$ as leaky boundary trap, not decoupled Klein-4 rings. No mode has $> 50\%$ mass in a single orbit. | < 30 s |

### 14.7 Exact $\sigma$-Parity and the Arithmetic Qubit (Section 11)

With the CRT qudit basis falsified and the Klein-4 orbit analysis showing a leaky trap rather than a clean basis, these scripts searched for the true conserved quantum number — and found it in the simplest symmetry of the system.

| Script | What it does and what it found | Runtime |
|---|---|---|
| `orbit_hamiltonian.py` | Projects $H(\alpha)$ into Klein-4 orbit space (124 orbits at $\varphi = 480$) to test the "arithmetic crystal" hypothesis. **Results:** Norm capture $44.1\%$ — orbits are a lossy basis. $O(2)$ is $99.2\%$ of the resonance. $0/124$ orbits fit in a single CRT block at odd $m$. Orbit Hamiltonian is $93.7\%$ off-diagonal by CRT pattern — initially interpreted as fundamental frustration. Klein-4 irrep decomposition: $(+,+)$ sector $= 44.1\%$, cross-irrep $= 49.2\%$. Commutator is a pure $\tau$-flipper with $0\%$ diagonal power. **This frustration was later resolved by the $\sigma$-exact discovery in `rabi_rabbit.py`.** | < 30 s |
| `rabi_rabbit.py` | **THE CENTRAL DISCOVERY SCRIPT.** Tests whether $[H(\alpha), P_\sigma] = 0$ exactly (not numerically) at both $m = 1155$ (odd) and $m = 2310$ (even). **Results ($m = 1155$):** (1) $\|[H(\alpha), P_\sigma]\|_F = 0.00 \times 10^{0}$ at all tested $\alpha$ values — algebraic identity confirmed. (2) Exact 240/240 split into $\sigma$-even and $\sigma$-odd eigenstates; zero mixed. (3) Perron ($1.000004$) and anti-Perron ($-0.744$) both in $\sigma$-even sector; $\sigma$-odd max eigenvalue only $0.337$. (4) $\tau$-parity smeared: only $1/480$ eigenstates has clean $\tau$-parity. (5) At $\alpha_c$: cross-irrep $= 50.7\%$ — the Rabi crossing. (6) Commutator-to-$D_{\text{sym}}$ $\tau$-off-diagonal ratio exactly $4.00$. **Results ($m = 2310$):** All results replicated — $\|[H(\alpha), P_\sigma]\|_F = 0.00 \times 10^{0}$, 240/240 split, zero mixed, Perron ($1.000003$) and anti-Perron ($-0.745$) in $\sigma$-even, $\sigma$-odd max $0.338$, cross-irrep $50.8\%$ at $\alpha_c$, ratio $4.00$. **Conclusion:** The qubit is $\tau$-parity within a $\sigma$-sector. Verified at both odd and even moduli. | < 60 s |
| `error_protection_test.py` | **Quantified noise robustness for Section 11.5.** Five tests at $m = 2310$ ($\varphi = 480$), seeded RNG (`np.random.default_rng(42)`). **Test 1 (Random Hermitian noise):** $\sigma$-breaking $\propto \varepsilon$ (linear). Flat band (counted at $|\lambda| < 0.01$; baseline $465/473$ sector states pass this threshold): $465 \to 463$ at $\varepsilon = 0.1$, $216$ at $\varepsilon = 0.3$, $67$ at $\varepsilon = 1.0$. **Test 2 ($\sigma$-symmetric noise):** $\|[H', P_\sigma]\|_F = 0.00$ at all $\varepsilon$ through $10$. Flat-band degradation matches Test 1. **Test 3 (Single-site perturbation):** Flat band $464/473$ at $\varepsilon = 10$; Perron shift $\Delta\lambda = +0.046$ at $\varepsilon = 1.0$. Confirms collective encoding. **Test 4 (Eigenvector $\sigma$-purity):** Shatters fast under generic noise — $480 \to 35$ ($|\sigma| > 0.99$) at $\varepsilon = 0.001$. **Test 5 (Perron state stability):** Overlap $99.998\%$ at $\varepsilon = 0.1$, $99.8\%$ at $\varepsilon = 1.0$, $98.1\%$ at $\varepsilon = 3.0$. | < 30 s |

### 14.8 Running the Scripts

**Modes (original sweep scripts):**
- `python <script>.py` — Full run: $m \in \{30, 210, 2310, 30030\}$
- `python <script>.py --medium` — Medium: $m \in \{30, 210, 2310\}$
- `python <script>.py --quick` — Quick: $m \in \{30, 210\}$

**Dependencies:** Python 3.10+, NumPy, Matplotlib. Parallel scripts additionally use `concurrent.futures` (stdlib) and benefit from scipy-openblas thread control via `OMP_NUM_THREADS`.

**Reproducibility:** A consumer PC with 8 GB RAM can run all scripts except `verify_510510_zero_defect.py` (which constructs a $92160 \times 92160$ matrix, requiring ~128 GB). The $m = 30030$ modulus requires $\sim 500$ MB for the $5760 \times 5760$ complex matrix. All computations are deterministic (no randomness), except `perturbation_robustness.py`, the random Hermitian baseline in `prime_resonance_perturbation.py`, and `error_protection_test.py` (all seeded RNG, `np.random.default_rng(42)`).

**Code availability.** All scripts described in this section are archived at [https://github.com/fancyland-llc/active-transport-prime-gas](https://github.com/fancyland-llc/active-transport-prime-gas).
---

## 15. Summary

| Result | Status | Evidence |
|---|---|---|
| Proposition 2.1: $[D_{\text{sym}}, P_\tau]$ is skew-symmetric | **Proved** | Algebraic proof; residual exactly 0.00 at all 4 primorial levels |
| Proposition 2.3: $H(\alpha)$ is Hermitian with real spectrum | **Proved** | Follows from skew-symmetry of $[D, P_\tau]$ |
| $1/\sqrt{2}$ commutator ratio: $\|[D,P_\tau]\|_F / \|D\|_F \to 1/\sqrt{2}$ | **Verified** | $0.528 \to 0.698 \to 0.706 \to 0.707$; derived from $\sqrt{2/3} / (2/\sqrt{3})$ |
| Alignment $A \to 3/4$: palindromic matrix retains $3/4$ of structure under inversion conjugation | **Verified** | $0.860 \to 0.757 \to 0.751 \to 0.750$ across 4 levels |
| Flat-band condensation: $\text{Var}(\lambda) \cdot \varphi \to 1.780$ | **Verified** | $1.624 \to 1.768 \to 1.777 \to 1.780$; flat band carries 99.7% of spectral weight at $m = 30030$ |
| Perron gap universality: $\Delta_P \to 0.8146$ | **Verified** | $0.857 \to 0.814 \to 0.8146 \to 0.8146$; lock-in by $m = 2310$ |
| Spectral width convergence: $W \to 1.5908$ | **Verified** | $1.634 \to 1.603 \to 1.592 \to 1.591$ |
| Gap ratio convergence: $\Delta_P / W \to 0.512$ | **Verified** | Close to but not exactly $1/2$ |
| Gap position phase transition at $\alpha_c$ | **Verified** | Width $10^{-7}$, derivative $-4{,}989{,}133$ at $m = 30030$; sharpens $38\times$ from $m = 2310$ |
| Non-primorial universality: $m = 385$ exhibits the crash | **Verified** | $\alpha_c(385) = 1.24227$, error $+0.30\%$ from $\sqrt{135/88}$ |
| Mod-4 partition non-multiplicativity ($4 \nmid m$) | **Verified** | Explicit counterexamples at $m = 30, 210, 2310$ |
| 5-adic torsion defect $\Delta v_5(2310) = +4$ (Bragg condition) | **Verified** | Exact arithmetic, confirmed by 5 independent methods |
| Observation 9.2: Perron–anti-Perron gap $\Delta\lambda = 1.7441 \pm 3.7 \times 10^{-4}$ is constant across crash | **Verified** | 203 sampled $\alpha$ values across $[\sqrt{3/2},\, \alpha_c]$; no drift, narrowing, or oscillation |
| Observation 9.3: Perron symmetry is a theorem-level consequence of flat-band condensation | **Verified** | Population imbalance $-2/5760$ predicts $\mathcal{A} \approx -3.5 \times 10^{-4}$; measured $-2.5 \times 10^{-6}$ |
| Anti-Perron chiral asymmetry at $m = 30030$: $\mathcal{A}_{\text{anti}} = -0.00296$ | **Observed** | Stable across crash window ($\sigma < 10^{-10}$); even-block dominant |
| Anti-Perron asymmetry sign flip at $m = 15015$: $\mathcal{A}_{\text{anti}} = +0.00227$ | **Verified** | 4-class analysis confirms genuine sign flip; sector-normalized $\mathcal{A} = +0.005$; Klein four-group mirror symmetry; spatial autocorrelation identical |
| Klein four-group symmetry in anti-Perron eigenvector at odd $m$ | **Verified** | $\mathrm{enrichment}(c) = \mathrm{enrichment}(c+2 \bmod 4)$ exactly; $\{0,3\}$ depleted, $\{1,2\}$ enriched; parity exactly conserved |
| Eigenvalue gap $\Delta\lambda \approx 1.7441$: dimension-determined | **Verified** | $\Delta\lambda(30030) = 1.74414$, $\Delta\lambda(15015) = 1.74416$; $0.0012\%$ difference; not torsion-sensitive |
| $T_i/T_e = 2$ torsion coincidence (Section 9.5) | **Falsified** | NBI heating artifact; ECRH-only L-H transitions on ASDEX Upgrade show $T_e > T_i$ (Ryter 2014 [10]) |
| IPR ratio $\approx 0.38 \ln(\varphi)$: logarithmic eigenvector localization | **Observed** | $1.43 \to 1.92 \to 2.10 \to 2.46 \to 2.47 \to 2.67 \to 2.88$ across $\varphi = 8 \to 5760$; 7 data points (4 primorial + 3 non-primorial); torsion–localization gap grows from $< 1\%$ ($\varphi = 480$) to $8\%$ ($\varphi = 5760$) |
| Conjecture 6.1: $\alpha_c = \sqrt{135/88}$ | **Conjectured** | Error $+2.56\% \to +0.30\% \to +0.027\% \to +0.015\%$; convergence sharply decelerating (local $\beta: 2.0 \to 0.24$), possibly logarithmic |
| Sylow torsion affects transition (Section 8.5) | **Resolved** | Same $\alpha_c$ ($\delta = +0.00018\%$) but $80\times$ wider transition; torsion smears, doesn't shift |
| Landau-Zener connection: near-complete non-adiabatic transfer at $\delta \sim 10^{-7}$ | **Mathematical observation** | $P_{LZ} \to 1$ for any sweep rate; physical realization requires identifying $\alpha$ with a time-dependent parameter |
| $\sqrt{2/3}$ quench hypothesis | **Falsified** | All 4 indicators pass smoothly through $\alpha = 0.8165$; no transition |
| Chirality balance hypothesis: $I = 0 \Rightarrow$ stronger localization | **Falsified** | $m = 30$ ($I = 0$) is the MOST dispersive; localization scales with $\varphi$, not chirality |
| $\sqrt{3/2}$ conjecture: $\kappa \cdot \alpha_c = 1$ | **Falsified** | $\alpha_c = 1.2389$ is $+1.16\%$ above $\sqrt{3/2}$; product $\kappa \cdot \alpha_c = \sqrt{45/44} \neq 1$ |
| Washing conjecture: primorial expansion clears all torsion defects | **Falsified** | $m = 510510$: $\Delta v_3 = \Delta v_5 = +16$ persist; non-multiplicativity of mod-4 partition is root cause |
| Perron eigenvector chiral asymmetry (Observation 4.4) | **Falsified** | $\mathcal{A} = -2.5 \times 10^{-6}$; $50.000\%/50.000\%$ mass split; flat-band forces uniformity |
| 13-periodic Fourier structure in Perron eigenvector (Observation 4.4) | **Falsified** | $100\%$ DC power; $0.0\%$ at every tested period; eigenvector is completely flat |
| MASER / stimulated topological emission mechanism | **Falsified** | $H(\alpha)$ is a static parametric family with no time axis; Rabi oscillations require a time-dependent Hamiltonian |
| Flat-band robustness to random noise (Section 10.1.1) | **Falsified** | At $\varepsilon = 0.03$ ($2.6\%$ of $\|A\|_F$), survival drops to $35.6\%$; fragile to unstructured perturbations |
| Selective screening of prime-periodic perturbations (Section 10.1.2) | **Verified** | IN primes $95.0\%$ vs OUT primes $70.0\%$ at $\varepsilon = 1.0$, $\alpha_c$; ratio $1.36\times$ |
| Perfect screening of $p = 2, 3$ (Section 10.1.3) | **Verified** | $100\%$ survival at all $\varepsilon$; coprime sieve eliminates those Fourier components from Hilbert space |
| Screening ratio amplification at $\alpha_c$ (Observation 10.1) | **Verified** | IN/OUT ratio $1.15$ at $\alpha = 0$ vs $1.36$ at $\alpha = \alpha_c$; critical point enhances discrimination |
| Hopping perturbation immunity (Section 10.1.4) | **Verified** | $99.1\%$–$99.2\%$ survival at $\varepsilon = 1.0$; IN/OUT distinction vanishes for off-diagonal perturbations |
| Flat-band obstruction for $\lvert r{=}1\rangle$ under $V_{13}$ (Section 10.2.1) | **Verified** | Swing $= 0.046$; $V_{13}$ eigenvalue plateau at $\lambda \approx +0.885$, std dev $= 0.121$; CRT fixed point |
| Noise-immune Perron–anti-Perron Rabi oscillation (Section 10.2.2) | **Falsified** | Period $3.51 \approx 2\pi/\Delta\lambda = 3.60$ = free precession; K3 kill — all primes give same period |
| Rabi oscillation immune to in-factorization noise (Section 10.2.3) | **Explained** | $V_5, V_7$ are exact scalars within CRT blocks (spread $= 3.2 \times 10^{-13}$); algebraic identity, not dynamical screening |
| Spectral architecture as arithmetic bandpass filter (Section 10.2.4) | **Superseded** | Replaced by CRT block qudit architecture (Section 10.3) |
| $V_2 = -I$, $V_3 = -\tfrac{1}{2}I$ on coprime residues (Section 10.3.1) | **Proved** | Exact; cosine is constant on coprimes iff $p \leq 3$; the $p = 3$ boundary |
| CRT block decomposition: 30 blocks × 16D (Section 10.3.2) | **Proved** | Joint eigenspaces of $(V_5, V_7, V_{11})$; all blocks dimension 16 (8D per $\sigma$-sector since $V_p$ commutes with $P_\sigma$); exact by CRT |
| Block-internal noise immunity: $V_5$ spread $= 3.2 \times 10^{-13}$ (Section 10.3.3) | **Proved** | Numerical zero ($3.2 \times 10^{-13}$, eigensolver tolerance); protection is algebraic, not approximate |
| Gate operators $V_{13}, V_{17}, V_{19}$ non-trivial within blocks (Section 10.3.4) | **Verified** | Noise $0.007$–$0.014$; max Rabi transfer $0.41$; 8-level qudit per $\sigma$-sector (16D pre-projection) |
| Block leakage $7.4\%$ at $t = 200$ (Section 10.3.5) | **Measured** | $92.6\%$ contained; $\sim 60$ gate operations before $10\%$ leakage |
| Identity 2.5: $D_{\text{sym}} = m/2 \cdot (J-I) - \lvert A \rvert$ (element-wise) | **Proved** | Frobenius residual exactly $0.00$ at 7 moduli ($m = 30, 105, 210, 1155, 2310$ via `chirality_analysis.py`; confirmed at $m = 2310, 1155$ via `flatband_hotspot_hunter.py` with element-wise `np.abs(A)`) |
| Identity 2.6: $[D_{\text{sym}}, P_\tau] = -[\lvert A \rvert, P_\tau]$ | **Proved** | Frobenius residual exactly $0.00$ at 7 moduli |
| Corollary 2.7: $H(\alpha)$ canonical decomposition via $\lvert A \rvert$ | **Proved** | Follows algebraically from Identities 2.5–2.6 |
| $\lvert A \rvert/m$ exact isospectrality for matched pairs | **Verified** | Max eigenvalue difference $< 10^{-13}$ at $\varphi = 8, 48, 480$ |
| $\lvert A \rvert/m$ non-conjugacy under coprime bijection $\sigma$ | **Verified** | Frobenius error $1.5$ ($\varphi = 8$) to $9.8$ ($\varphi = 48$); same eigenvalues, different eigenvectors |
| $P_\tau$ exact conjugacy under $\sigma$ | **Verified** | $\lVert P_{\tau,\text{odd}} - \sigma P_{\tau,\text{even}} \sigma^\top \rVert_F = 0.0$ at all sizes |
| Counterfactual $P_\tau$ swap: $D_{\text{sym}}$ determines IPR | **Verified** | Swapping $P_\tau$ changes eigenvalues by $0.0$ and IPR by $0.0$; $D_{\text{sym}}$ alone controls localization |
| Kurtosis explosion at odd $m$: $117$ vs $6.9$ ($\varphi = 480$) | **Verified** | Odd coprime gaps create topological hotspots; max IPR $21.85/n$ (odd) vs $6.26/n$ (even) |
| Additive–multiplicative orthogonality: $\lvert A \rvert_{\text{ew}}$ vs CRT blocks | **Verified** | $\lvert A \rvert_{\text{ew}}$ is 97.1% off-diagonal in CRT basis ($m = 2310$), 98.7% ($m = 1155$); additive chirality and multiplicative CRT are independent decompositions |
| Flat-band kurtosis source: entirely from flat band | **Verified** | Flat-band-only kurtosis: 115.1 (odd $m = 1155$) vs 4.1 (even $m = 2310$); all top-10 localized states are flat-band eigenvectors |
| Boundary localization at odd $m$: mirror pair $r = 2$, $r = m{-}2$ | **Verified** | Most localized mode ($\text{IPR} \times \varphi = 21.8$) carries 10.66% each at $r = 2$ and $r = 1153$; $\tau$-symmetry forces mirror pairing |
| Block Participation Ratio: extreme-IPR modes partially CRT-spread | **Verified** | Avg BPR of top-20 localized flat-band states: $20.0/30$ ($m = 2310$), $38.1/60$ ($m = 1155$); ratio $\approx 0.63{-}0.67$ of $N_{\text{blocks}}$ |
| Klein-4 orbit decomposition: $\langle \sigma, \tau \rangle \cong \mathbb{Z}_2 \times \mathbb{Z}_2$ | **Verified** | 124 orbits at $\varphi = 480$; orbit sizes match group theory (116 size-4, 8 size-2); $O(2) = \{2, 577, 578, 1153\}$ confirmed via $\tau(2) = 578$, $\sigma(2) = 1153$ |
| $O(2)$ leaky trap at odd $m$: 41.3% mass concentration | **Verified** | Top flat-band mode ($\text{IPR} \times \varphi = 21.8$) puts 41.3% in $O(2)$, OPR $= 5.53/124$; $O(2)$ mass is $8\times$ background mean ($5.26\%$); 4/10 top modes dominated by $O(2)$ |
| $O(2)$ deletion at even $m$: coprime sieve erases boundary trap | **Verified** | $\gcd(2, 2310) = 2 \neq 1$ deletes $O(2)$; top mode OPR jumps to $26.1/124$ ($5\times$ delocalization); peak orbit mass drops to $10.9\%$; kurtosis 114 $\to$ 3.9 |
| Flat band is coupled orbit network, not decoupled Klein-4 rings | **Falsified (decoupled)** | 0/473 modes have $> 50\%$ mass in single orbit; 0/473 have $> 90\%$; modes spread across 5–26 orbits depending on $O(2)$ existence |
| $\tau$-parity non-diagonalization: $\langle v \lvert P_\tau \rvert v \rangle \neq \pm 1$ | **Verified** | Top-10 modes at odd $m$: $\langle v \lvert P_\tau \rvert v \rangle \in [0.005, 0.49]$; commutator $[\lvert A \rvert, P_\tau]$ rotates between $\tau$-eigenspaces. At even $m$ (no $O(2)$ funnel): values up to $0.82$ — $\tau$-parity partially restored |
| Theorem 11.1: $[H(\alpha), P_\sigma] = 0$ for all $\alpha$ ($\sigma$-parity conservation) | **Proved** | Algebraic identity; $\|[H(\alpha), P_\sigma]\|_F = 0.00 \times 10^{0}$ at all tested $\alpha$ and moduli; 240/240 $\sigma$-even/$\sigma$-odd split with zero mixed states |
| $\sigma$-even sector contains Perron and anti-Perron | **Verified** | Perron eigenvalue $1.000004$ and anti-Perron $-0.744$ both in $\mathcal{H}_+$; $\sigma$-odd max eigenvalue only $0.337$ |
| $\tau$-parity is dynamical Rabi variable | **Verified** | Only 1/480 eigenstates (Perron) has clean $\tau$-parity; all others $\langle P_\tau \rangle \in [-0.75, +1.00]$; commutator generates continuous rotations between $\tau$-eigenspaces |
| Phase transition = $\tau$-mixing 50/50 crossover at $\alpha_c$ | **Verified** | $(+,+)$ share: $77.7\% \to 43.7\% \to 18.9\%$ as $\alpha: 0 \to \alpha_c \to 2.5$; cross-irrep share $= 50.7\%$ at $\alpha_c$ |
| Commutator-to-$D_{\text{sym}}$ $\tau$-off-diagonal ratio $= 4.00$ | **Verified** | Integer ratio at both moduli; suggests algebraic origin (Problem 17) |
| Klein-4 irrep dimension split: $124 + 116$ ($\sigma$-even) | **Verified** | Gap $124 - 116 = 8$ = number of size-2 Klein-4 orbits; detuning between $\tau$-subspaces |
| Orbit Hamiltonian: frustrated crystal interpretation | **Superseded** | $44.1\%$ norm capture, $93.7\%$ CRT off-diagonal; resolved by $\sigma$-exact block-diagonalization (Theorem 11.1) |
| CRT blocks are qudit basis | **Superseded** | CRT blocks are Zeeman sublevels within each $\tau$-parity level, not the computational basis; true qubit is $\tau$-parity within $\sigma$-sector |
| Four-layer arithmetic error protection | **Verified** | Layer 1: $\sigma$-conservation (exact for palindromic-symmetric noise; linear breaking otherwise). Layer 2: coprime sieve ($Q \sim 10^7$). Layer 3: flat-band robust to $\varepsilon \lesssim 0.1$ ($6\%$ of $\|H\|_F$); collapses at $\varepsilon = 0.3$. Layer 4: single-site Perron overlap $99.8\%$ at $\varepsilon = \|H\|_F$; flat band $464/473$ at $\varepsilon = 10$. See `error_protection_test.py`. |

---

## 16. References

[1] A. P. Matos, "The Universal Two-Prime Formula and Spectral Structure of Palindromic Distance Matrices," preprint (2026). DOI: [10.5281/zenodo.19210625](https://doi.org/10.5281/zenodo.19210625).

[2] A. P. Matos, "Primorial Thermodynamics: Gauss Sums, the IPB98 Scaling Law, and the Hierarchy Break at $m = 30030$," preprint (2026). DOI: [10.5281/zenodo.19188924](https://doi.org/10.5281/zenodo.19188924).

[3] A. P. Matos, "Spectral Isotropy and the Exact Temperature of the Prime Gas," preprint (2026). DOI: [10.5281/zenodo.19156532](https://doi.org/10.5281/zenodo.19156532).

[4] G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed. Oxford University Press (2008).

[5] R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed. Cambridge University Press (2012).

[6] P. W. Anderson, "Absence of Diffusion in Certain Random Lattices," *Phys. Rev.* **109** (1958), 1492–1505.

[7] T. M. Apostol, *Introduction to Analytic Number Theory*. Springer (1976).

[8] D. Leykam, A. Andreanov, and S. Flach, "Artificial flat band systems: from lattice models to experiments," *Adv. Phys.: X* **3** (2018), 1473052.

[9] Y. Andrew, N. C. Hawkes, et al. "Edge Ion Parameters at the L-H Transition on JET," *Plasma Phys. Control. Fusion* **46** (2004), A87–A94.

[10] F. Ryter, L. Barrera Orte, B. Kurzan, R. M. McDermott, G. Tardini, E. Viezzer, M. Bernert, R. Fischer, and the ASDEX Upgrade Team, "Experimental evidence for the key role of the ion heat channel in the physics of the L-H transition," *Nucl. Fusion* **54** (2014), 083003.

[11] I. I. Rabi, "Space Quantization in a Gyrating Magnetic Field," *Phys. Rev.* **51** (1937), 652–654.

[12] C. Zener, "Non-Adiabatic Crossing of Energy Levels," *Proc. R. Soc. Lond. A* **137** (1932), 696–702.

[13] A. Peres, "Stability of quantum motion in chaotic and regular systems," *Phys. Rev. A* **30** (1984), 1610–1615.

---

## Acknowledgments

This work was conducted using the Lattice OS axiomatic research protocol, in which the author proposed hypotheses and two large language models — Claude Opus 4 (Anthropic) and Gemini 3.1 Pro (Google DeepMind) — served as adversarial reviewers. All claims were resolved by computational execution, not by model assertion. The LLM agents contributed to three aspects of the work: building and executing the $35{+}$ computational scripts that test each claim, running independent parallel sessions that converged on the CRT block architecture from different starting points, and helping falsify the author's own hypotheses (the primorial qubit, the CRT qudit basis, the washing conjecture, the chirality balance) so that corrected results could be pursued within the same working session. The final discovery — that $\sigma$-parity is exactly conserved and $\tau$-parity provides the 2-level degree of freedom mathematically isomorphic to a qubit — emerged from four successive wrong predictions (CRT blocks, Klein-4 orbits, frustrated crystal, free precession), each of which eliminated a wrong basis until the clean structure emerged. Whether a physical system exists whose Hamiltonian realizes $H(\alpha)$ remains an open question (Section 12).

---

<div align="center">

_Fancyland LLC — Lattice OS research infrastructure._

_The rabbit has been caught._

</div>
