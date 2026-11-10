# Mass gap, confinement and the θ = π critical point of 1+1D QED

A working theory with eight quantitative predictions, each tested against an
independent first-principles lattice calculation (`qed2.py`, exact
diagonalisation, no fitted physics parameters).

## 1. The problem, and what this is not

The central unsolved problem of quantum field theory is non-perturbative: why a
gauge theory whose Lagrangian contains only massless fields has a **mass gap**
and **confines** charge. In four dimensions (Yang–Mills, QCD) that is open, and
nothing here solves it.

What is here is the one gauge theory where the mechanism can be followed end to
end — QED in 1+1 dimensions, the Schwinger model — pushed past its exactly
solvable corner into the regime that has no closed-form solution: finite fermion
mass and vacuum angle θ = π, where the theory has a quantum phase transition.

    L = −¼ F² + ψ̄ (i∂̸ − g A̸ − m) ψ + (g θ / 4π) ε^{μν} F_{μν}

g has dimension of mass; the only parameters are m/g and θ.

## 2. The theory

**Step 1 — the anomaly generates the gap.** The axial current is anomalous,
∂·j₅ = (g/π) E. In 1+1D j₅^μ = ε^{μν} j_ν, and Gauss's law is ∂₁E = g j⁰.
Together these close into (□ + g²/π) E = 0: the only gauge-invariant excitation
of the massless theory is a free boson of mass M = g/√π.

**Step 2 — bosonize.** With ψ̄γ^μψ = ε^{μν}∂_νφ/√π the full theory is a single
scalar:

    H = ½π² + ½(∂ₓφ)² + (g²/2π)(φ + θ/(2√π))² − c m g cos(2√π φ),   c = e^γ/(2π^{3/2})

θ = 2πα is exactly the background electric field made by external charges ±αg
at the ends of space. Shifting φ moves θ into the cosine:
V(φ) = (g²/2π) φ² − c m g cos(2√π φ − θ).

**Step 3 — screening versus confinement.** At m = 0 the shift removes θ
entirely: the vacuum energy does not depend on α, so the force between external
charges vanishes at large distance. At m ≠ 0 the cosine pins φ and the energy
density depends on α: a linear potential, except for integer α where a
fermion pair breaks the string.

**Step 4 — θ = π.** Here V(φ) = (g²/2π) φ² + c m g cos(2√π φ) is symmetric under
φ → −φ (charge conjugation). For small m the quadratic term wins and the vacuum
is unique. For large m the cosine wins and there are two vacua with electric
field ±g/2: charge conjugation breaks spontaneously. One real scalar with a Z₂
symmetry and a double well in 1+1D means the transition, if continuous, is in
the **2D Ising** universality class (c = ½).

**Step 5 — the boundary conditions (the sharp part).** Put the theory on an
interval with the background field g/2 entering at one end and leaving at the
other. Then the ordered vacuum "+g/2" costs nothing at the walls, while "−g/2"
needs a unit charge at each wall, costing ~2m. Both walls therefore pin the
Ising order parameter to the *same* sign. In boundary CFT, the Ising model with
equal fixed boundary conditions has the simplest possible spectrum — the
Virasoro tower of the identity alone:

    Z = χ₀(q) = 1 + q² + q³ + 2q⁴ + 2q⁵ + …
    E_n − E₀ = π v h_n / (L + ℓ),   h = 2, 3, 4, 4, 5, …
    E₀ = e·L + b − π c v / (24 (L + ℓ))

Everything about the spectrum at criticality is then fixed with no free
parameter. (ℓ is a non-universal extrapolation length; it cancels from ratios.)

### Predictions

| | Prediction | Origin |
|---|---|---|
| P1 | Gap at m = 0: M/g = 1/√π = 0.56419 | anomaly |
| P2 | dM/dm at m = 0 equals e^γ = 1.781 | cosine, 1st order |
| P3 | String tension at m = 0 is zero for every α | shift symmetry |
| P4 | σ/g² = c (m/g)(1 − cos 2πα) + b(α)(m/g)², with c = 0.15993, b(¼) = −0.358, b(½) = 0 | cosine, 1st + 2nd order |
| P5 | At θ = π there is a critical mass m_c with a gapless spectrum | Z₂ breaking |
| P6 | Level ratios at m_c: (E_n − E₀)/(E₁ − E₀) = 1, 3/2, 2, 2, 5/2 | identity tower |
| P7 | Casimir energy gives c = ½ | Ising |
| P8 | Speed of light v → 1 as a → 0 | Lorentz invariance |

The second-order coefficient in P4 uses the constant E₊ = −8.9139 from mass
perturbation theory (Adam 1997), **quoted from memory**; the overall sign is
fixed by requiring the second-order vacuum energy to be negative. The lattice
fit of b below is independent of it.

Two cheap analytic estimates of m_c from Step 4, to see how far "simple" theory
gets:
- classical (V''(0) = 0): m_c/g = 1/(4π²c) = **0.158**, second order;
- Gaussian variational (self-consistent normal ordering): **0.265**, but
  *first* order — the wrong kind of transition.

## 3. The test

Lattice Hamiltonian with staggered fermions on N sites and open boundaries;
Gauss's law is solved exactly so only N spins remain. Zero-charge sector,
Lanczos with full reorthogonalisation, up to N = 22 (705 432 states; N = 24 exceeds
the memory available for the Lanczos basis). x = 1/(ga)², physical length Lg = N/√x. The lattice mass is
m − g²a/8 (Dempsey et al. 2022), which removes O(a) errors. The solver is checked
against an independently built dense Hamiltonian (`test_qed2.py`).

### Away from θ = π

| | Prediction | Lattice |
|---|---|---|
| P1 | 0.56419 | 0.5636 (extrapolated in a²) / 0.5641 (in a), from x = 9, 12, 16 |
| P2 | 1.781 | 1.753 (x = 2), 1.798 (x = 4); 1.874 at x = 9, where Lg ≤ 6.7 is too small |
| P3 | 0 | |σ|/g² ≤ 1×10⁻⁴ at x = 2, 3 (α = ¼); < 10⁻⁵ at α = ½ |
| P4, α = ¼, linear coefficient | 0.1599 | 0.1588 (x = 3), 0.1661 (x = 2) |
| P4, α = ¼, quadratic coefficient | −0.358 | −0.328 (x = 3), −0.394 (x = 2) |
| P4, α = ½, linear coefficient | 0.3199 | 0.3268 (x = 3), 0.3247 (x = 2) |
| P4, α = ½, quadratic coefficient | 0 | −0.025 (x = 3), −0.034 (x = 2) |

Point by point at x = 3, lattice/theory for σ is 0.996 – 1.014 at α = ¼ and
1.013 – 1.019 at α = ½ for m/g ≤ 0.1.

### At θ = π

The critical mass m_c(N) is taken where R = (E₂−E₀)/(E₁−E₀) for sizes N−2 and N
cross; c and v are then evaluated there (`analyze_crit.py`).

| ga | m_c/g (N = 22) | R (N = 22) | R (N → ∞, in 1/N²) | c (N = 22) | v (N = 22) |
|---|---|---|---|---|---|
| 2.00 | 0.3955 | 1.4675 | 1.486 | 0.413 | 1.022 |
| 1.41 | 0.3588 | 1.4855 | 1.488 | 0.485 | 0.989 |
| 1.00 | 0.3361 | 1.5218 | 1.500 | 0.681 | 0.917 |

Full tower at ga = 1.41, m/g = 0.356 (`run_tower.py`):

| N | E₂ | E₃ | E₄ | E₅ |
|---|---|---|---|---|
| 12 | 1.4919 | 1.9434 | 1.9567 | 2.3422 |
| 16 | 1.4947 | 1.9683 | 1.9743 | 2.4125 |
| 20 | 1.4968 | 1.9814 | 1.9846 | 2.4454 |
| **P6** | **1.5** | **2** | **2** | **2.5** |

(columns are (E_n − E₀)/(E₁ − E₀)).

## 4. Verdict

**Confirmed**
- P1 to 0.1%, P2 to 1–2% where the volume is adequate, P3 to 10⁻⁴, P4 to 1–2% at first order and ~10% on the
  second-order coefficient. The 15–22% deficit in the fractional-charge string
  tension reported in the first pass of this work was the missing second-order
  term, not a lattice problem.
- P5: the gap closes as 1/L at a critical mass on every lattice spacing studied.
- P6: the spectrum is the identity tower. The ratios 3/2, 2, 2, 5/2 are
  reproduced to 0.2%, 0.9%, 0.8%, 2.2% at N = 20 and are still approaching, and
  the predicted two-fold degeneracy at h = 4 is there (splitting 0.3%, shrinking
  with N). The alternatives are far away: free boundaries would give
  E₂/E₁ = 3, opposite fixed boundaries 2.

**Supported, not tight**
- P6 via the crossing method gives R = 1.486 – 1.500 across the three lattice
  spacings: within 1% of 3/2, but at ga = 1.41, where the sequence has
  essentially converged, it sits 0.7% low. Not explained at these sizes.
- P7: c = 0.485 at ga = 1.41 (3% low, drifting slowly); the other two lattice
  spacings are still moving (0.41 rising, 0.68 falling) and bracket ½.
- P8: v = 0.99, 1.02, 0.92 — close to 1 but with no clear trend in a, because
  the three volumes are not comparable.

**Failed**
- Both analytic estimates of m_c. Classical gives 0.158; Gaussian gives 0.265
  and predicts a first-order transition, which the 1/L gap and conformal tower
  rule out.

**Not determined**
- The continuum value of m_c. The lattice values 0.396, 0.359, 0.336 (N = 22;
  N → ∞ estimates 0.392, 0.358, 0.345) are at coarse spacings ga = 2, 1.41, 1
  and are still moving with a. They head toward the literature value
  m_c/g ≈ 0.3335 (Byrnes et al. 2002, quoted from memory) but 22 sites cannot
  reach the continuum at the volume this transition needs (Lg ≳ 15); that
  requires tensor-network sizes.
- P2 at fine lattice spacing: the x = 9 value is 5% high and the volume there
  is too small to tell lattice error from finite-size error.

## 5. What is known and what is this work

Known: Steps 1–4 (Schwinger 1962; Coleman, Jackiw, Susskind 1975; Coleman 1976),
the Ising character of the θ = π transition and the value of m_c.

This work: Step 5 — the identification of the open-chain critical spectrum with
the identity Virasoro tower from the physical argument that the background
field pins both walls to the same vacuum — and its numerical test, P6–P7. I
have not searched the literature for prior statements of this, so treat the
novelty as unverified; the test itself stands either way.

Limits: the mechanism of Step 1 is special to 1+1D, where the gauge field has
no propagating degrees of freedom. This is a controlled laboratory for mass
generation and confinement, not a route to the 4D proof.

## 6. Reproduce

    gcc -O3 -march=native -fopenmp -shared -fPIC kernel.c -o kernel.so
    python3 test_qed2.py          # solver vs dense Hamiltonian
    python3 run_gap.py            # P1            (~6 min)
    python3 run_slope.py          # P2            (~3 min)
    python3 run_string.py         # P3, P4        (~10 min), then analyze_string.py
    python3 run_crit.py 8,10,12,14,16,18,20,22   # P5–P8 (~40 min), then analyze_crit.py
    python3 run_tower.py          # P6 tower      (~5 min)

Raw numbers are in `data/`.
