# Summary

The full write-up is the paper: [paper/paper.pdf](paper/paper.pdf) (25 pages; source `paper/paper.md`).
This file is the one-page version.

**Subject.** QED in 1+1 dimensions (Schwinger model): how a gauge theory with only massless
fields gets a mass gap, when it confines, and what happens at the θ = π phase transition.
It is a controlled laboratory; it does not address the 4D Yang–Mills mass-gap problem.

**Chain of reasoning.** axial anomaly → gap M = g/√π → bosonized Hamiltonian → charge is
screened at m = 0 and fractional charge confined at m ≠ 0 → at θ = π charge conjugation
breaks above a critical mass (Ising class) → on an open chain the background field pins both
walls to the same vacuum, so the critical spectrum is the Ising boundary-CFT identity tower
(charge 0) or ε tower (charge −1).

**Tests.** Two independent codes on the lattice Hamiltonian, agreeing to 1e-11 where they
overlap: exact diagonalisation (≤ 24 sites) and DMRG with an exact matrix-product Coulomb
term (≤ 192 sites). No fitted physics parameters.

| Prediction | Theory | Lattice |
|---|---|---|
| Gap at m = 0 (continuum) | 0.56419 | 0.56415 – 0.56421 |
| dM/dm (continuum) | 1.781 | 1.775 |
| String tension at m = 0 | 0 | < 5×10⁻⁵ g² |
| String tension, charge ¼: linear, quadratic coefficient | 0.1599, −0.358 | 0.158 – 0.161, −0.33 to −0.36 |
| Level ratios, charge 0 (64 sites) | 3/2, 2, 2, 5/2 | 1.498, 1.993, 1.995, 2.487 |
| Level ratios, charge −1 (64 sites) | 2, 3, 4, 4 | 1.996, 2.986, 3.969, 3.986 |
| Ratio of first two gaps at the crossing, N → ∞ | 3/2 | 1.500(2) |
| Central charge | 1/2 | 0.500(7) Casimir; 0.49(1) entanglement |
| β/ν | 1/8 | 0.124(3) |
| 1/ν | 1 | 1.00(1) |
| Speed of light, continuum | 1 | 1.001(4) |
| Critical mass, continuum (no prediction) | — | 0.3333(2); literature 0.333561(4) |

**Findings beyond confirmation.**
- Finite-size corrections at the critical point are non-monotonic: the gap ratio is flat to
  four digits for 16–24 sites at a value 0.9% off, then moves to the correct limit.
- On the staggered lattice the even–odd alternation of the entanglement entropy is
  proportional to the order parameter, so single-parity fits give central charges near 0.63
  and 0.37; the parity average gives 0.50.
- Classical and Gaussian estimates of the critical mass fail (0.158; 0.265 with the wrong
  order of transition).

**Remaining limitations.** All uncertainties are systematic (spread between extrapolation
choices), not statistical. The critical mass is limited by the volumes reachable at the two
finest lattice spacings. dM/dm is 0.3% below theory, unexplained beyond the two-point
extrapolation in m. No continuum extrapolation of the string-tension coefficients is quoted.

**Prior work.** The boundary-CFT identification and a far more precise critical mass were
published by Arguello Cruz, Tarnopolsky and Xin (arXiv:2412.01902). This work rederives and
verifies them with independent codes, and adds the central-charge measurements, the
entanglement observation, the second-order string-tension test and the Gaussian analysis.
