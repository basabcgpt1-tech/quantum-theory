# Summary

The full write-up is the paper: [paper/paper.pdf](paper/paper.pdf) (source `paper/paper.md`).
This file is the one-page version.

**Subject.** QED in 1+1 dimensions (Schwinger model): how a gauge theory with only massless
fields gets a mass gap, when it confines, and what happens at the θ = π phase transition.
It is a controlled laboratory; it does not address the 4D Yang–Mills mass-gap problem.

**Chain of reasoning.** axial anomaly → gap M = g/√π → bosonized Hamiltonian → charge is
screened at m = 0 and fractional charge confined at m ≠ 0 → at θ = π charge conjugation
breaks above a critical mass (Ising class) → on an open chain the background field pins both
walls to the same vacuum, so the critical spectrum is the Ising boundary-CFT identity tower
(charge 0) or ε tower (charge −1).

**Tests** (exact diagonalisation, up to 24 sites, no fitted physics parameters):

| Prediction | Theory | Lattice |
|---|---|---|
| Gap at m = 0 | 0.56419 | 0.5636 – 0.5641 |
| dM/dm | 1.781 | 1.753, 1.798 |
| String tension at m = 0 | 0 | < 10⁻⁴ |
| String tension, charge ¼: linear, quadratic coefficient | 0.1599, −0.358 | 0.1588, −0.328 |
| Level ratios, charge 0 | 3/2, 2, 2, 5/2 | 1.497, 1.981, 1.985, 2.445 |
| Level ratios, charge −1 | 2, 3, 4, 4 | 1.965, 2.884, 3.739, 3.846 |
| Central charge (Casimir) | 1/2 | 0.50(4) |
| β/ν | 1/8 | 0.126 (using the extrapolation length from the gap) |
| ν | 1 | 1.0 ± 0.1 |

**What failed or is open.**
- Classical and Gaussian estimates of the critical mass (0.158; 0.265 with the wrong order).
- Central charge from entanglement entropy does not converge at these sizes (0.26 – 0.43).
- The level ratio at the crossing point is 1.486, 0.9% below 3/2, and stable; unexplained
  at these sizes.
- No continuum value of the critical mass: lattice values are 0.392, 0.358, 0.345 at
  ga = 2, 1.41, 1; the literature continuum value is 0.333561(4).

**Prior work.** The boundary-CFT identification was published by Arguello Cruz, Tarnopolsky
and Xin (arXiv:2412.01902). It was rederived here independently and then found in the
literature; this work is a verification by a different method, plus the Casimir central
charge, the second-order string tension test, and the Gaussian analysis.
