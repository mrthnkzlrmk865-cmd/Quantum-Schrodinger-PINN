# Literature review and research gap

This review separates **established knowledge** (from the cited literature)
from what this project does. Every reference below was checked for
authors, venue, year and DOI or arXiv identifier. Where a numerical
prefactor depends on conventions, this is stated, and the project only
uses the order of magnitude.

## 1. Anderson localization and transport in disordered media

* **Anderson (1958)** showed that a quantum particle in a sufficiently
  disordered lattice can stop diffusing: its eigenstates become
  exponentially localized. This happens without any classical trapping, by
  interference of multiply scattered waves.
* **Abrahams, Anderson, Licciardello & Ramakrishnan (1979)**, scaling theory:
  in one and two dimensions all states are localized for arbitrarily weak
  uncorrelated disorder (orthogonal symmetry class). A mobility edge
  exists only in 3D.
* **Lee & Ramakrishnan (1985)**, review of disordered electronic systems.
  In weakly disordered 2D systems the localization length is exponentially
  large, ξ ≈ ℓ exp(π k ℓ / 2), where ℓ is the mean free path and k the wave
  number. So 2D localization is unobservable in finite samples unless
  kℓ is of order 1.
* **Kramer & MacKinnon (1993)** and **Evers & Mirlin (2008)** review
  localization theory, numerical methods (transfer matrix, level
  statistics, multifractality) and Anderson transitions.
* **Beenakker (1997)**, random-matrix theory of quantum transport in
  quasi-1D wires (DMPK equation). In a wire with N open channels the
  localization length is of order N ℓ; the exact prefactor depends on the
  symmetry class and on how ξ is defined. Transport is ballistic for
  L ≪ ℓ, diffusive (Ohmic, ⟨T⟩ ∝ ℓ/L) for ℓ ≪ L ≪ Nℓ, and localized
  (ln T ∝ −L/ξ) for L ≫ Nℓ. In the localized regime T is broadly
  (log-normally) distributed, so ⟨ln T⟩ and ⟨T⟩ behave differently.
* **Akkermans & Montambaux (2007)**, textbook: Born approximation for the
  elastic and transport mean free paths, and diffusion and weak
  localization of waves.
* **Lagendijk, van Tiggelen & Wiersma (2009)**, overview of 50 years of
  Anderson localization across electrons, light, sound and matter waves.

### Correlated disorder and matter waves

* **Kuhn, Sigwarth, Miniatura, Delande & Müller (2007)**: coherent
  matter-wave transport in correlated (speckle) potentials in 2D and 3D.
  The correlation length enters through kζ; for smooth disorder
  (kζ ≫ 1) scattering is forward-peaked and ℓ_tr ≫ ℓ_s.
* Experiments: **Billy et al. (2008)** and **Roati et al. (2008)** observed
  Anderson localization of Bose–Einstein condensates in 1D disorder
  (speckle and quasi-periodic). **Schwartz, Bartal, Fishman & Segev (2007)**
  observed transverse localization of light in 2D disordered photonic
  lattices. **Sanchez-Palencia & Lewenstein (2010)** review disordered
  quantum gases.
* **Shinozuka & Deodatis (1991)**: the spectral-representation (random
  phase) method for simulating stationary Gaussian random fields. This
  project's disorder generator is a version of it.

## 2. Numerical methods for the time-dependent Schrödinger equation

* **Feit, Fleck & Steiger (1982)** and **Kosloff & Kosloff (1983)**:
  split-operator / Fourier spectral propagation of the TDSE.
* **Strang (1968)**: symmetric operator splitting, second order in Δt.
  **Jahnke & Lubich (2000)** give error bounds for exponential operator
  splittings.
* **Muga, Palao, Navarro & Egusquiza (2004)** review complex absorbing
  potentials (CAPs) for open-boundary wave propagation. CAPs reflect
  long-wavelength components more strongly, which is visible in this
  project's validation V5.

## 3. Physics-informed neural networks

* **Raissi, Perdikaris & Karniadakis (2019)**: the PINN framework. It
  already includes a nonlinear Schrödinger example (1D, periodic).
* **Karniadakis et al. (2021)**: review of physics-informed machine learning.
* **Rahaman et al. (2019)**: spectral bias of neural networks.
  **Tancik et al. (2020)**: Fourier features. **Wang, Wang & Perdikaris (2021)**:
  Fourier features for multiscale PDEs in PINNs.
* Failure modes and training pathologies: **Wang, Teng & Perdikaris (2021)**
  on gradient imbalance between loss terms; **Krishnapriyan et al. (2021)**,
  who show PINNs can fail on mildly hard convection and reaction problems
  because of optimization, not network capacity; **Wang, Sankaran &
  Perdikaris (2024)** on causality violation in time-dependent problems
  and causal loss weighting.
* Collocation sampling: **Wu, Zhu, Tan, Kartha & Lu (2023)**, a systematic
  study of uniform vs residual-based adaptive (RAD, RAR-D) sampling.
* PINNs for the linear TDSE: **Shah, Stiller, Hoffmann & Cangi (2022,
  NeurIPS ML4PS workshop)** study PINN solvers for the time-dependent
  Schrödinger equation on harmonic-oscillator problems and their
  generalization.

## 4. Research gap (honest assessment)

**Established and not claimed here:** Anderson localization, the scaling
theory, the Born mean free paths of correlated disorder, quasi-1D
transport regimes, split-step spectral solvers, CAPs, PINNs, Fourier
features, adaptive sampling and causal weighting.

**What this project actually contributes** is a careful, reproducible
*computational study*, not a new theory or algorithm:

1. A validated time-dependent simulation of wave-packet transport through
   a finite slab of **Gaussian-correlated** disorder in a **2D hard-wall
   waveguide** with few open channels (N = 6), with ensemble statistics
   over many realizations.
2. A quantitative comparison of the measured trends (vs W, ζ, k, L)
   against the **perturbative Born / diffusion estimates**, identifying
   where the perturbative picture holds and where it fails.
3. A feasibility analysis of whether localization can be distinguished from
   finite-size and classical-trapping effects at this system size.
4. An assessment of PINN accuracy on the **same disordered potentials**
   (identical realization for PINN and reference), including documented
   failure modes, i.e. where PINNs are reliable and where they are not.

We did not find, and do not claim the absence of, prior studies that
combine exactly these elements. The literature search covered the topics
above but was not exhaustive. The PINN part in particular is an
*evaluation* of existing methods on a new test problem, not a new method.

## References

1. P. W. Anderson, "Absence of diffusion in certain random lattices," *Phys. Rev.* **109**, 1492 (1958). doi:10.1103/PhysRev.109.1492
2. E. Abrahams, P. W. Anderson, D. C. Licciardello, T. V. Ramakrishnan, "Scaling theory of localization: absence of quantum diffusion in two dimensions," *Phys. Rev. Lett.* **42**, 673 (1979). doi:10.1103/PhysRevLett.42.673
3. P. A. Lee, T. V. Ramakrishnan, "Disordered electronic systems," *Rev. Mod. Phys.* **57**, 287 (1985). doi:10.1103/RevModPhys.57.287
4. B. Kramer, A. MacKinnon, "Localization: theory and experiment," *Rep. Prog. Phys.* **56**, 1469 (1993). doi:10.1088/0034-4885/56/12/001
5. F. Evers, A. D. Mirlin, "Anderson transitions," *Rev. Mod. Phys.* **80**, 1355 (2008). doi:10.1103/RevModPhys.80.1355
6. C. W. J. Beenakker, "Random-matrix theory of quantum transport," *Rev. Mod. Phys.* **69**, 731 (1997). doi:10.1103/RevModPhys.69.731
7. E. Akkermans, G. Montambaux, *Mesoscopic Physics of Electrons and Photons*, Cambridge University Press (2007).
8. A. Lagendijk, B. van Tiggelen, D. S. Wiersma, "Fifty years of Anderson localization," *Physics Today* **62**(8), 24 (2009). doi:10.1063/1.3206091
9. R. C. Kuhn, O. Sigwarth, C. Miniatura, D. Delande, C. A. Müller, "Coherent matter wave transport in speckle potentials," *New J. Phys.* **9**, 161 (2007). doi:10.1088/1367-2630/9/6/161
10. J. Billy et al., "Direct observation of Anderson localization of matter waves in a controlled disorder," *Nature* **453**, 891 (2008). doi:10.1038/nature07000
11. G. Roati et al., "Anderson localization of a non-interacting Bose–Einstein condensate," *Nature* **453**, 895 (2008). doi:10.1038/nature07071
12. T. Schwartz, G. Bartal, S. Fishman, M. Segev, "Transport and Anderson localization in disordered two-dimensional photonic lattices," *Nature* **446**, 52 (2007). doi:10.1038/nature05623
13. L. Sanchez-Palencia, M. Lewenstein, "Disordered quantum gases under control," *Nature Physics* **6**, 87 (2010). doi:10.1038/nphys1507
14. M. Shinozuka, G. Deodatis, "Simulation of stochastic processes by spectral representation," *Appl. Mech. Rev.* **44**(4), 191 (1991). doi:10.1115/1.3119501
15. M. D. Feit, J. A. Fleck Jr., A. Steiger, "Solution of the Schrödinger equation by a spectral method," *J. Comput. Phys.* **47**, 412 (1982). doi:10.1016/0021-9991(82)90091-2
16. D. Kosloff, R. Kosloff, "A Fourier method solution for the time dependent Schrödinger equation as a tool in molecular dynamics," *J. Comput. Phys.* **52**, 35 (1983). doi:10.1016/0021-9991(83)90015-3
17. G. Strang, "On the construction and comparison of difference schemes," *SIAM J. Numer. Anal.* **5**, 506 (1968). doi:10.1137/0705041
18. T. Jahnke, C. Lubich, "Error bounds for exponential operator splittings," *BIT Numer. Math.* **40**, 735 (2000). doi:10.1023/A:1022396519656
19. J. G. Muga, J. P. Palao, B. Navarro, I. L. Egusquiza, "Complex absorbing potentials," *Phys. Rep.* **395**, 357 (2004). doi:10.1016/j.physrep.2004.03.002
20. M. Raissi, P. Perdikaris, G. E. Karniadakis, "Physics-informed neural networks: A deep learning framework for solving forward and inverse problems involving nonlinear partial differential equations," *J. Comput. Phys.* **378**, 686 (2019). doi:10.1016/j.jcp.2018.10.045
21. G. E. Karniadakis, I. G. Kevrekidis, L. Lu, P. Perdikaris, S. Wang, L. Yang, "Physics-informed machine learning," *Nat. Rev. Phys.* **3**, 422 (2021). doi:10.1038/s42254-021-00314-5
22. N. Rahaman et al., "On the spectral bias of neural networks," *ICML* (2019). arXiv:1806.08734
23. M. Tancik et al., "Fourier features let networks learn high frequency functions in low dimensional domains," *NeurIPS* 33 (2020). arXiv:2006.10739
24. S. Wang, H. Wang, P. Perdikaris, "On the eigenvector bias of Fourier feature networks: From regression to solving multi-scale PDEs with physics-informed neural networks," *Comput. Methods Appl. Mech. Eng.* **384**, 113938 (2021). doi:10.1016/j.cma.2021.113938
25. S. Wang, Y. Teng, P. Perdikaris, "Understanding and mitigating gradient flow pathologies in physics-informed neural networks," *SIAM J. Sci. Comput.* **43**, A3055 (2021). doi:10.1137/20M1318043
26. A. S. Krishnapriyan, A. Gholami, S. Zhe, R. M. Kirby, M. W. Mahoney, "Characterizing possible failure modes in physics-informed neural networks," *NeurIPS* 34 (2021). arXiv:2109.01050
27. S. Wang, S. Sankaran, P. Perdikaris, "Respecting causality for training physics-informed neural networks," *Comput. Methods Appl. Mech. Eng.* **421**, 116813 (2024). doi:10.1016/j.cma.2024.116813
28. C. Wu, M. Zhu, Q. Tan, Y. Kartha, L. Lu, "A comprehensive study of non-adaptive and residual-based adaptive sampling for physics-informed neural networks," *Comput. Methods Appl. Mech. Eng.* **403**, 115671 (2023). doi:10.1016/j.cma.2022.115671
29. K. Shah, P. Stiller, N. Hoffmann, A. Cangi, "Physics-informed neural networks as solvers for the time-dependent Schrödinger equation," NeurIPS 2022 Workshop on Machine Learning and the Physical Sciences. arXiv:2210.12522
30. D. C. Liu, J. Nocedal, "On the limited memory BFGS method for large scale optimization," *Math. Program.* **45**, 503 (1989). doi:10.1007/BF01589116
31. B. Efron, R. J. Tibshirani, *An Introduction to the Bootstrap*, Chapman & Hall (1993).

### Removed from the earlier bibliography

* Rebentrost et al. (2009), "Environment-assisted quantum transport," and
  Plenio & Huelga (2008), "Dephasing-assisted transport." Both are real,
  correct references, but they concern **dephasing / open-system**
  dynamics, which this project does not model. They were the motivation
  of the earlier (ENAQT) narrative and no longer support the work.
