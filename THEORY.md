# Mass gap and confinement in a gauge theory: a tested working theory in 1+1D QED

## The problem

The central open problem of quantum field theory is non-perturbative: why does a
gauge theory with massless fields in its Lagrangian have a **mass gap** and
**confine** charge (the 4D Yang–Mills / QCD problem). Perturbation theory cannot
see either effect at any order.

This repo does **not** solve the 4D problem. It takes the one gauge theory where
the mechanism can be derived exactly — QED in 1+1 dimensions (the Schwinger
model) — states the mechanism as a theory with quantitative predictions, and
checks those predictions from first principles with an independent lattice
Hamiltonian calculation (`schwinger.py`, numpy only).

## The theory

L = −¼ F² + ψ̄(i∂̸ − gA̸ − m)ψ, with g of mass dimension 1.

1. **Anomaly ⇒ gap.** The axial current is anomalous: ∂·j₅ = (g/π) E. In 1+1D
   j₅^μ = ε^{μν} j_ν, so with Gauss's law ∂₁E = g j⁰ this closes into
   (□ + g²/π) E = 0. The only gauge-invariant excitation of the massless theory
   is a free boson of mass **M = g/√π**. The "photon" eats the fermion pair.
2. **Bosonized form.** With ψ̄γ^μψ = ε^{μν}∂_νφ/√π the full theory is
   L = ½(∂φ)² − (g²/2π)(φ + θ/(2√π))² + c·m·g·cos(2√π φ), c = e^γ/(2π^{3/2}).
   A background field θ = 2πα is what a pair of external charges ±αg creates.
3. **Predictions.**
   - P1: massless gap M/g = 1/√π = 0.5642.
   - P2: small fermion mass shifts it linearly, dM/dm = e^γ = 1.781.
   - P3: at m = 0 the shift φ → φ − const absorbs θ: string tension σ = 0 for
     every α (charge is **screened**, not confined).
   - P4: at m ≠ 0 the cosine blocks that shift: σ = c·m·g·(1 − cos 2πα), a linear
     potential for fractional charge (**confinement**), vanishing for integer α
     (string breaking).

So: the gap comes from the anomaly, and confinement of fractional charge is
controlled by the explicit chiral-symmetry breaking m.

## The test

Staggered fermions on N sites, open boundaries, Gauss's law solved so the gauge
field is eliminated exactly; zero-charge sector diagonalised by Lanczos
(up to N = 20, dimension 184 756). Lattice mass m_lat = m − g²a/8 is used,
which removes the O(a) error. Volume limit from gap² = M² + c/N² over
N = 14…20; x = 1/(ga)².

| quantity | prediction | lattice (x = 1 … 6) |
|---|---|---|
| P1: M/g at m = 0 | 0.5642 | 0.5646 – 0.5724 (0.1% – 1.5% high) |
| P2: dM/dm | 1.781 | 1.75 (x=1) – 1.80 (x=6), from m/g = 0, 0.1 |
| P3: σ/g² at m = 0, α = ¼, ½ | 0, 0 | −0.0006, 0.0000 |
| P4: σ/g², m/g = 0.05, α = ½ | 0.0160 | 0.0167 |
| P4: σ/g², m/g = 0.10, α = ½ | 0.0320 | 0.0333 |
| P4: σ/g², m/g = 0.05, α = ¼ | 0.0080 | 0.0068 |
| P4: σ/g², m/g = 0.10, α = ¼ | 0.0160 | 0.0124 |

String tension is measured at x = 4 from the bulk energy density difference
between N = 16 and N = 20. Raw output: `run_N20.log`, `results.json`.

## What holds and what does not

- P1–P3 and P4 at α = ½ are confirmed at the 1–5% level.
- P4 at α = ¼ is 15–22% low. Not resolved: it is a single lattice spacing with
  small volumes, and the formula is only first order in m/g.
- The gap **rises** with x instead of converging: at fixed N ≤ 20 the physical
  box N/√x shrinks as a → 0, so the large-x points carry finite-volume error.
  The a → 0 extrapolations printed by the script (0.570 linear in a², 0.597
  quadratic) are therefore worse than the best single points and should not be
  quoted as the continuum value. Fixing this needs larger N (tensor networks).
- The mechanism is special to 1+1D: the gauge field has no propagating degrees
  of freedom and the anomaly equation closes. In 4D Yang–Mills neither is true;
  this is a controlled laboratory for the problem, not a proof of it.

## Reproduce

    python3 schwinger.py 20   # ~2 min; `18` runs in ~20 s
