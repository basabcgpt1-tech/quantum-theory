---
title: "Mass generation, confinement and the θ = π Ising critical point of two-dimensional QED: derivations and an exact-diagonalisation test"
date: "10 November 2026"
lang: en
---

::: {.abstract}
**Abstract.** Quantum electrodynamics in 1+1 dimensions (the Schwinger model) is the simplest gauge theory that generates a mass gap from massless constituents, confines fractional charge, and — at vacuum angle $\theta=\pi$ — undergoes a quantum phase transition that has no closed-form solution. We give a self-contained derivation of the chain *anomaly → mass gap → bosonized Hamiltonian → screening versus confinement → $Z_2$ symmetry breaking → boundary conformal field theory*, and turn it into twelve quantitative predictions. Each is tested against an independent first-principles calculation: exact diagonalisation of the gauge-invariant Kogut–Susskind lattice Hamiltonian on open chains of up to 24 sites, with no fitted physics parameters. The massless gap is reproduced to $0.1\%$ ($M/g = 0.5636$–$0.5641$ against $1/\sqrt\pi = 0.5642$); the string tension of fractional charges agrees with mass perturbation theory to $1$–$2\%$ once the second-order term is included; and at the critical point the low-lying spectrum of the neutral sector is the Virasoro tower of the identity, with level ratios $1.497, 1.981, 1.985, 2.445$ against $\tfrac32, 2, 2, \tfrac52$, while the charge $-1$ sector gives $1.965, 2.884, 3.739, 3.846$ against $2, 3, 4, 4$. The Casimir energy gives central charge $c = 0.48$–$0.49$ on the best-converged lattice and the order-parameter exponent is $\beta/\nu = 0.126$ (Ising: $\tfrac12$, $\tfrac18$). We also record what fails: the classical and Gaussian-variational estimates of the critical mass ($0.158$ and $0.265$, the latter with the wrong order of transition), the entanglement-entropy estimate of $c$ at these sizes, and the continuum limit of the critical mass, which is out of reach of exact diagonalisation. The boundary-CFT identification tested here was established by Arguello Cruz, Tarnopolsky and Xin (2024); our contribution is an independent derivation and verification with a different method, and tests of the central charge and second-order string tension. Nothing here bears on the four-dimensional Yang–Mills mass-gap problem beyond providing a controlled laboratory for its mechanisms.
:::

# Introduction

## The problem

The defining difficulty of quantum field theory is that its most important phenomena are invisible in perturbation theory. Quantum chromodynamics has a Lagrangian containing only massless gluons and nearly massless quarks, yet its spectrum consists of massive, colour-neutral hadrons. Proving that four-dimensional Yang–Mills theory has a mass gap is an open problem; so is deriving confinement. No expansion in the gauge coupling can produce either effect at any finite order.

Progress on such questions has historically come from models in which the mechanism can be followed from the Lagrangian to the spectrum. The most instructive is QED in one space and one time dimension, introduced by Schwinger [@schwinger] precisely to show that gauge invariance does not require a massless vector particle. In two dimensions the gauge coupling $g$ has the dimension of mass, the gauge field has no propagating degrees of freedom, and with massless fermions the model is exactly solvable [@schwinger; @lowenstein]: the spectrum is a free boson of mass $g/\sqrt\pi$. With a fermion mass $m$ and a vacuum angle $\theta$ the model is no longer solvable, and it acquires a rich structure [@cjs; @coleman76]: confinement of fractional charge, string breaking, and at $\theta=\pi$ a line of first-order transitions ending at a critical point $m_c$ in the two-dimensional Ising universality class [@hamer82; @byrnes].

## Scope and claims

This paper has three aims.

1. **Derivation.** Sections 2–4 derive, from the Lagrangian, every formula that is later tested, including the finite-size spectrum at the critical point from boundary conformal field theory (BCFT).
2. **Test.** Sections 5–7 test those formulas against exact diagonalisation of the lattice Hamiltonian. The lattice calculation uses none of the continuum results as input.
3. **Honest accounting.** Section 8 separates what is confirmed, what is merely consistent, and what fails.

We state plainly what is and is not new. The physics of Sections 2–3 is classic [@schwinger; @cjs; @coleman76; @adam]. The identification of the open-chain critical spectrum with Ising boundary conditions, and its use to locate $m_c$ to five digits by tensor-network methods, is due to Arguello Cruz, Tarnopolsky and Xin [@act]; we arrived at the same identification independently by the physical argument of Section 4 before consulting that work, and our numerical results agree with theirs where they overlap. What this paper adds is (i) a verification by a different numerical method, with full degeneracy structure of two towers; (ii) an extraction of the central charge from the Casimir energy, which Ref. [@act] takes as input; (iii) a lattice test of the second-order mass correction to the string tension; and (iv) an explicit demonstration that the Gaussian variational method predicts the wrong order for the transition. We make no claim about four-dimensional gauge theory.

# The continuum theory

## Lagrangian and conventions

We use metric $(+,-)$, $\epsilon^{01}=+1$, $\gamma^5=\gamma^0\gamma^1$, and

$$\mathcal L = -\tfrac14 F_{\mu\nu}F^{\mu\nu} + \bar\psi\left(i\gamma^\mu\partial_\mu - g\gamma^\mu A_\mu - m\right)\psi - \frac{g\theta}{4\pi}\,\epsilon^{\mu\nu}F_{\mu\nu}. \qquad (2.1)$$

The electric field is $E\equiv F_{01}$, so $F^{\mu\nu}=-\epsilon^{\mu\nu}E$ and the last term is $(g\theta/2\pi)E$, a total derivative. The parameters are $m/g$ and $\theta$; the theory is super-renormalisable and $g$ sets the scale.

In the gauge $A_1=0$ the field $A_0$ is not dynamical. Its equation of motion is Gauss's law,

$$\partial_1 E = g\,j^0, \qquad j^\mu=\bar\psi\gamma^\mu\psi, \qquad (2.2)$$

whose solution is $E(x)=g\int^x j^0 + F$. The integration constant $F$ is a background field. Coleman [@coleman76] showed that it is the same thing as the vacuum angle, $F=g\theta/2\pi$, and that physics is periodic in $\theta$ with period $2\pi$: a background field larger than $g/2$ in magnitude is lowered by pair production. Equivalently, static external charges $\pm\alpha g$ placed at the ends of the system shift $\theta\to\theta+2\pi\alpha$ between them. The Hamiltonian is

$$H=\int dx\left[\bar\psi\left(-i\gamma^1\partial_1+m\right)\psi+\tfrac12E^2\right]. \qquad (2.3)$$

Eliminating $E$ gives a linear potential $-\tfrac{g^2}{4}\int\!\!\int j^0(x)\,|x-y|\,j^0(y)$ between charges. Classically, then, all charge is confined. Quantum mechanically the answer is subtler, and depends on $m$.

## The anomaly generates the mass gap

Set $m=0$. Classically both the vector current $j^\mu$ and the axial current $j_5^\mu=\bar\psi\gamma^\mu\gamma^5\psi$ are conserved. In two dimensions they are not independent: $\gamma^\mu\gamma^5=-\epsilon^{\mu\nu}\gamma_\nu$, so

$$j_5^\mu=-\epsilon^{\mu\nu}j_\nu. \qquad (2.4)$$

Quantum mechanically, gauge invariance forces the axial current to be anomalous,

$$\partial_\mu j_5^\mu=\frac{g}{\pi}E. \qquad (2.5)$$

(The overall sign depends on conventions for $\epsilon$ and $\gamma^5$; the sign relative to the coupling in Eq. (2.1) is fixed, and is the one that makes Eq. (2.7) non-tachyonic.) Maxwell's equation $\partial_\mu F^{\mu\nu}=g j^\nu$ reads $-\epsilon^{\mu\nu}\partial_\mu E=g j^\nu$. Taking the divergence of Eq. (2.4) and using $\epsilon^{\mu\nu}\epsilon_{\lambda\nu}=-\delta^\mu_\lambda$,

$$\partial_\mu j_5^\mu=-\epsilon^{\mu\nu}\partial_\mu j_\nu=-\frac1g\,\Box E. \qquad (2.6)$$

Equating Eqs. (2.5) and (2.6),

$$\left(\Box+\frac{g^2}{\pi}\right)E=0. \qquad (2.7)$$

The electric field is a free massive field. Since every gauge-invariant local operator of the massless model can be built from $E$, the entire spectrum is a free boson of mass

$$M_0=\frac{g}{\sqrt\pi}\approx0.56419\,g. \qquad (2.8)$$

No massless particle survives: the photon has absorbed the fermion–antifermion pair. Three features of this derivation are special to two dimensions and should be kept in mind when drawing analogies: the identity (2.4); the absence of transverse photons; and the fact that Eqs. (2.5)–(2.6) close on a single field.

## Bosonization

The same result follows, together with much more, from abelian bosonization [@coleman75; @mandelstam]. A Dirac fermion in 1+1 dimensions is equivalent to a real scalar $\phi$ with the dictionary

$$j^\mu=\frac{1}{\sqrt\pi}\epsilon^{\mu\nu}\partial_\nu\phi,\qquad \bar\psi\, i\gamma^\mu\partial_\mu\psi\to\tfrac12(\partial\phi)^2,\qquad \bar\psi\psi=-\frac{e^{\gamma}}{2\pi}\,\mu\,N_\mu\!\left[\cos 2\sqrt\pi\phi\right], \qquad (2.9)$$

where $\gamma=0.5772\ldots$ is Euler's constant and $N_\mu$ denotes normal ordering with respect to a free field of mass $\mu$. The vector current is identically conserved, and Gauss's law (2.2) integrates to

$$E=\frac{g}{\sqrt\pi}\phi+\frac{g\theta}{2\pi}=\frac{g}{\sqrt\pi}\left(\phi+\frac{\theta}{2\sqrt\pi}\right). \qquad (2.10)$$

Substituting into Eq. (2.3) and choosing $\mu=M_0$,

$$H=\int dx\;N_{M_0}\!\left[\tfrac12\Pi^2+\tfrac12(\partial_x\phi)^2+\frac{M_0^2}{2}\left(\phi+\frac{\theta}{2\sqrt\pi}\right)^2-c\,m g\cos\left(2\sqrt\pi\phi\right)\right], \qquad (2.11)$$

with $c=e^\gamma/2\pi^{3/2}\approx0.15993$.

Shifting $\phi\to\phi-\theta/2\sqrt\pi$ moves $\theta$ into the cosine:

$$V(\phi)=\frac{g^2}{2\pi}\phi^2-c\,m g\cos\left(2\sqrt\pi\phi-\theta\right). \qquad (2.12)$$

At $m=0$ this is a free field of mass $M_0$, reproducing Eq. (2.8), and $\theta$ has disappeared entirely.

## Screening, confinement and mass perturbation theory

The vacuum energy density $\mathcal E(\theta)$ controls the force between external charges $\pm\alpha g$ a distance $r$ apart: for large $r$ the energy cost is $\sigma r$ with string tension

$$\sigma(\alpha)=\mathcal E(\theta+2\pi\alpha)-\mathcal E(\theta). \qquad (2.13)$$

**Massless case.** $\mathcal E$ is independent of $\theta$, so $\sigma=0$ for every $\alpha$. Classical confinement is completely undone: the vacuum polarises and screens any charge. This is the same physics as the mass gap, seen from another side — a massive photon mediates a short-range force.

**Massive case.** To first order in $m$ the energy shift is the expectation value of the perturbation in the free vacuum, and $\langle N_{M_0}\cos(2\sqrt\pi\phi-\theta)\rangle_0=\cos\theta$, so $\mathcal E^{(1)}=-c\,mg\cos\theta$. The second-order term was computed by Adam [@adam]. In his notation $\Sigma=e^\gamma M_0/2\pi=c\,g$ and

$$\mathcal E(\theta)=-m\Sigma\cos\theta-\frac{m^2\Sigma^2}{4M_0^2}\left(\mathcal A_+\cos2\theta+\mathcal A_-\right)+O(m^3), \qquad (2.14)$$

with $\mathcal A_+=-8.9139$, $\mathcal A_-=9.7384$, where $\mathcal A_\pm=M_0^2\int d^2x\,(e^{\mp2K_0(M_0|x|)}-1)$. (Ref. [@adam] quotes $\epsilon=\ln Z/V$, which is minus the energy density; the sign in Eq. (2.14) is such that the $\theta$-averaged second-order energy is negative, as it must be.) Hence, at $\theta=0$,

$$\frac{\sigma(\alpha)}{g^2}=c\,\frac mg\left(1-\cos2\pi\alpha\right)+b(\alpha)\left(\frac mg\right)^2+\ldots,\qquad b(\alpha)=\frac{e^{2\gamma}}{16\pi^2}\,\mathcal A_+\left(1-\cos4\pi\alpha\right). \qquad (2.15)$$

Numerically $b(\tfrac14)=-0.358$ and $b(\tfrac12)=0$. We recomputed $\mathcal A_+$ by direct quadrature and reproduce $-8.9139$ (`check_formulas.py`). Fractional charges are confined by a linear potential; integer charges ($\alpha\in\mathbb Z$) are not, because a dynamical fermion pair breaks the string. Expanding the cosine in Eq. (2.12) about its minimum gives the boson mass

$$M^2=M_0^2+4\pi c\,mg\cos\theta+O(m^2)\quad\Longrightarrow\quad \frac{dM}{dm}\Big|_{m=0,\theta=0}=e^\gamma\approx1.7811. \qquad (2.16)$$

# The point θ = π

## Symmetry breaking

At $\theta=\pi$ Eq. (2.12) becomes

$$V(\phi)=\frac{g^2}{2\pi}\phi^2+c\,mg\cos2\sqrt\pi\phi, \qquad (3.1)$$

which is even in $\phi$. The symmetry $\phi\to-\phi$ is charge conjugation, which maps $E\to-E$ and is a symmetry only when the background field is $0$ or $g/2$. For small $m$ the parabola dominates, the vacuum is unique and $\langle E\rangle=0$ in the shifted variable — that is, the fermions half-screen the background field. For large $m$ the cosine dominates and there are two degenerate vacua near $\phi=\pm\sqrt\pi/2$, with electric field $\pm g/2$: these are the classical states "background field unscreened" and "background field reversed by one pair". Charge conjugation is spontaneously broken. One real scalar field with a $Z_2$ symmetry interpolating between a single and a double well is the Landau–Ginzburg description of the two-dimensional Ising model, so if the transition is continuous it must be in the Ising universality class, with central charge $c=\tfrac12$ and exponents $\nu=1$, $\beta=\tfrac18$ [@coleman76; @hamer82; @byrnes]. The order parameter is the electric field.

## Two analytic estimates of the critical mass

**Classical.** The curvature of Eq. (3.1) at the origin vanishes at

$$\frac{m_c}{g}\Big|_{\mathrm{cl}}=\frac{1}{4\pi^2c}=\frac{1}{2\sqrt\pi\,e^{\gamma}}\approx0.1584. \qquad (3.2)$$

This is identical to the point where the first-order mass formula (2.16) gives $M^2=0$, as it must be, since both keep one power of the normal-ordered cosine.

**Gaussian variational.** The classical estimate ignores the fact that the cosine is normal-ordered with respect to mass $M_0$, while near the transition the fluctuating field is much lighter. A Gaussian trial state with mean $\phi_0$ and mass $\mu$ accounts for this [@stevenson]. Using Coleman's re-ordering identities [@coleman75],

$$\begin{aligned}N_{M}\!\left[\tfrac12\Pi^2+\tfrac12\phi'^2\right]&=N_\mu\!\left[\tfrac12\Pi^2+\tfrac12\phi'^2\right]+\frac{\mu^2-M^2}{8\pi},\qquad N_M[\phi^2]=N_\mu[\phi^2]+\frac{1}{4\pi}\ln\frac{M^2}{\mu^2},\\ N_M[\cos\beta\phi]&=\left(\frac{\mu}{M}\right)^{\beta^2/4\pi}N_\mu[\cos\beta\phi],\end{aligned} \qquad (3.3)$$

with $\beta^2=4\pi$, the energy density of the trial state in units of $M_0^2$ is

$$\frac{\mathcal E(\mu,\phi_0)}{M_0^2}=\frac{\mu^2-1}{8\pi}-\frac{\ln\mu^2}{8\pi}+\frac{\phi_0^2}{2}+k\,\mu\cos2\sqrt\pi\phi_0,\qquad k=\pi c\,\frac mg, \qquad (3.4)$$

where $\mu$ is measured in units of $M_0$. Stationarity in $\mu$ gives $\mu=\tfrac12\left(-a+\sqrt{a^2+4}\right)$ with $a=4\pi k\cos2\sqrt\pi\phi_0$. A continuous transition would require the curvature in $\phi_0$ at the origin to vanish, $1=4\pi k\mu$, simultaneously with the gap equation at $\phi_0=0$, $\mu-1/\mu+4\pi k=0$. Eliminating $k$ gives $1/\mu=0$: **there is no solution**. The Gaussian approximation therefore cannot produce a continuous transition here. Minimising Eq. (3.4) numerically (Fig. 3b) one finds instead that a pair of minima at $\phi_0\approx\pm0.63$ drops below the symmetric minimum at

$$\frac{m_c}{g}\Big|_{\mathrm{Gauss}}\approx0.265\quad\text{(first order)}. \qquad (3.5)$$

Both estimates are tested, and found wanting, in Section 7.

![Effective potentials at $\theta=\pi$. (a) Classical potential, Eq. (3.1), below, at and above the classical critical mass. (b) Gaussian effective potential, Eq. (3.4) minimised over $\mu$: the symmetric minimum never becomes unstable, and the transition is first order.](fig/fig3_potential.svg){width=100%}

# Finite-size spectrum at criticality from boundary CFT

The lattice calculation works on an interval of length $L$, so we need the prediction for the spectrum in a finite box with definite boundary conditions. For a conformal field theory on a strip with boundary conditions $a$, $b$, Cardy [@cardy86; @cardy89] showed that the Hamiltonian is

$$H_{ab}=e_\infty L+e_a+e_b+\frac{\pi v}{L}\left(L_0-\frac{c}{24}\right), \qquad (4.1)$$

with $e_\infty$, $e_{a,b}$ non-universal bulk and boundary energies, $v$ the velocity, and $L_0$ acting on a direct sum of irreducible Virasoro representations determined by the pair $(a,b)$. The partition function is a sum of characters, $Z_{ab}=\sum_i n^i_{ab}\chi_i(q)$, $q=e^{-\pi v/LT}$. The Casimir term $-\pi cv/24L$ is the Blöte–Cardy–Nightingale–Affleck result [@bcn; @affleck].

For the Ising model ($c=\tfrac12$) there are three conformal boundary conditions — spin fixed up ($+$), fixed down ($-$), and free ($f$) — and

$$Z_{++}=\chi_0,\qquad Z_{+-}=\chi_{1/2},\qquad Z_{ff}=\chi_0+\chi_{1/2},\qquad Z_{+f}=\chi_{1/16}, \qquad (4.2)$$

$$\chi_0=q^{-1/48}\left(1+q^2+q^3+2q^4+2q^5+3q^6+\ldots\right),\qquad \chi_{1/2}=q^{1/2-1/48}\left(1+q+q^2+q^3+2q^4+2q^5+\ldots\right). \qquad (4.3)$$

**Which boundary condition does the gauge theory realise?** Consider the Schwinger model on an interval with the background field $g/2$ entering at the left wall, and total charge $Q$ on the interval, so that the field leaving at the right wall is $g(\tfrac12+Q)$. In the ordered phase the two bulk vacua have field $+g/2$ and $-g/2$.

- $Q=0$: the field is $+g/2$ at both walls. The vacuum $+g/2$ matches both walls at no cost. The vacuum $-g/2$ requires a unit charge at each wall, costing $\approx2m$. Both walls therefore act as a boundary magnetic field of the same sign on the Ising order parameter. Under the renormalisation group a boundary field flows to the fixed boundary condition, so the critical theory has boundary conditions $(+,+)$.
- $Q=-1$: the field is $+g/2$ on the left and $-g/2$ on the right. The walls favour opposite vacua: $(+,-)$.
- $Q=+1$: the field leaving is $3g/2$, which is not a vacuum. The extra charge binds to the right wall at an $O(g)$ energy cost, after which the bulk again sees $(+,+)$.

From Eqs. (4.1)–(4.3), the predictions at $m=m_c$ are parameter-free:

$$Q=0:\quad \frac{E_n-E_0}{E_1-E_0}=1,\;\tfrac32,\;2,\;2,\;\tfrac52,\ldots\qquad\qquad Q=-1:\quad \frac{E_n-E_0}{E_1-E_0}=1,\;2,\;3,\;4,\;4,\ldots \qquad (4.4)$$

$$E_1-E_0=\frac{2\pi v}{L+\ell},\qquad E_0=e_\infty L+2e_b-\frac{\pi c\,v}{24(L+\ell)}\qquad(Q=0). \qquad (4.5)$$

Here $\ell$ is an extrapolation length. It arises because the leading irrelevant boundary operator for a fixed boundary is the stress tensor, whose effect at first order is to rescale all levels by a common factor $1-\ell/L$; ratios are therefore corrected only at relative order $L^{-2}$. Eliminating $v$ and $\ell$ between the two relations (4.5) gives

$$E_0=e_\infty L+2e_b-\frac{c}{48}\left(E_1-E_0\right), \qquad (4.6)$$

an estimator of $c$ that needs neither. By contrast, free boundaries would give $(E_2-E_0)/(E_1-E_0)=3$.

The identification (4.4) of the $Q=0$ and $Q=\pm1$ sectors with $(+,+)$ and $(+,-)$ was made, and used as a criterion to locate $m_c$, in Ref. [@act].

**Further Ising predictions.** With $(+,+)$ boundaries the order parameter profile at criticality is $\langle\sigma(x)\rangle\propto[(L/\pi)\sin(\pi x/L)]^{-1/8}$, so the mid-chain electric field scales as

$$\langle E\rangle_{\mathrm{mid}}\propto(L+\ell)^{-\beta/\nu},\qquad\beta/\nu=\tfrac18, \qquad (4.7)$$

and the entanglement entropy of a block of length $x$ adjacent to a wall is [@cc]

$$S(x)=\frac c6\ln\left[\frac{2L}{\pi}\sin\frac{\pi x}{L}\right]+\text{const}. \qquad (4.8)$$

Near the critical point any dimensionless ratio is a function of $(m-m_c)L^{1/\nu}$ with $\nu=1$.

## Summary of predictions

| | Prediction | Value | Source |
|---|---|---|---|
| P1 | massless gap $M_0/g$ | $0.56419$ | Eq. (2.8) |
| P2 | $dM/dm$ at $m=0$ | $1.7811$ | Eq. (2.16) |
| P3 | string tension at $m=0$ | $0$ for all $\alpha$ | Sec. 2.4 |
| P4 | $\sigma$, linear coefficient | $0.15993\,(1-\cos2\pi\alpha)$ | Eq. (2.15) |
| P5 | $\sigma$, quadratic coefficient | $-0.358$ ($\alpha=\tfrac14$), $0$ ($\alpha=\tfrac12$) | Eq. (2.15) |
| P6 | critical point exists at $\theta=\pi$, gap $\propto1/L$ | — | Sec. 3.1 |
| P7 | $Q=0$ level ratios | $\tfrac32,2,2,\tfrac52$ | Eq. (4.4) |
| P8 | $Q=-1$ level ratios | $2,3,4,4$ | Eq. (4.4) |
| P9 | central charge | $\tfrac12$ | Eq. (4.6) |
| P10 | $\beta/\nu$ | $\tfrac18$ | Eq. (4.7) |
| P11 | $\nu$ | $1$ | Sec. 4 |
| P12 | $v\to1$ in the continuum limit | $1$ | Lorentz invariance |

# Lattice formulation

## Hamiltonian

We use the Kogut–Susskind Hamiltonian [@ks; @bsk] with one-component staggered fermions $\chi_n$ on sites $n=1,\ldots,N$ (spacing $a$) and link variables $U_n=e^{i\vartheta_n}$ with conjugate integer-spaced electric fields $L_n$, $[\vartheta_n,L_m]=i\delta_{nm}$:

$$H=-\frac{i}{2a}\sum_{n=1}^{N-1}\left(\chi_n^\dagger U_n\chi_{n+1}-\text{h.c.}\right)+m_{\mathrm{lat}}\sum_{n=1}^N(-1)^n\chi_n^\dagger\chi_n+\frac{g^2a}{2}\sum_{n=1}^{N-1}L_n^2. \qquad (5.1)$$

Even sites carry the upper (particle) component and odd sites the lower (antiparticle) component, so an occupied even site has charge $+1$ and an empty odd site charge $-1$. Gauss's law is

$$L_n-L_{n-1}=\chi_n^\dagger\chi_n-\tfrac12\left[1-(-1)^n\right]. \qquad (5.2)$$

With open boundaries and $L_0=\varepsilon_0\equiv\theta/2\pi$ fixed, Eq. (5.2) determines every $L_n$ in terms of the fermion occupations, and the link phases can be removed by the gauge transformation $\chi_n\to\prod_{k<n}U_k^\dagger\chi_n$. The gauge field is thereby eliminated exactly; what remains is its Coulomb interaction. A Jordan–Wigner transformation $\chi_n=\prod_{k<n}(i\sigma^z_k)\,\sigma^-_n$ gives a spin chain. In units $W=2H/(ag^2)$,

$$W=x\sum_{n=1}^{N-1}\left(\sigma^+_n\sigma^-_{n+1}+\text{h.c.}\right)+\frac{\mu}{2}\sum_{n=1}^{N}\left[1+(-1)^n\sigma^z_n\right]+\sum_{n=1}^{N-1}L_n^2, \qquad (5.3)$$

$$L_n=\varepsilon_0+\frac12\sum_{k=1}^{n}\left[\sigma^z_k+(-1)^k\right],\qquad x=\frac{1}{g^2a^2},\qquad\mu=\frac{2m_{\mathrm{lat}}}{g^2a}=2\sqrt x\,\frac{m_{\mathrm{lat}}}{g}. \qquad (5.4)$$

Energies in units of $g$ are $W/2\sqrt x$; the physical length is $Lg=N/\sqrt x$; the total charge is $Q=\tfrac12\sum_n\sigma^z_n$ (for even $N$) and is conserved. The continuum limit is $x\to\infty$ at fixed $N/\sqrt x\to\infty$. The term $\sum L_n^2$ is a long-range (linear) interaction between all pairs of spins — the lattice form of the Coulomb potential of Section 2.1.

## Mass shift

Dempsey, Klebanov, Pufu and Zan [@dkpz] observed that at $m_{\mathrm{lat}}=0$ the lattice Hamiltonian is not at the point corresponding to the massless continuum theory: translation by one site, which is the lattice remnant of the discrete chiral symmetry, is a symmetry only if

$$m_{\mathrm{lat}}=m-\frac{g^2a}{8}\qquad\Longleftrightarrow\qquad\mu=2\sqrt x\,\frac mg-\frac14. \qquad (5.5)$$

With this identification lattice artefacts start at $O(a^2)$ rather than $O(a)$. We use Eq. (5.5) throughout; $m$ always denotes the physical mass.

# Numerical method

**Basis.** For even $N$ and charge $Q$ the basis consists of all bit strings with $N/2+Q$ up spins; the dimension is $\binom{N}{N/2+Q}$, which is $705\,432$ at $N=22$ and $2\,704\,156$ at $N=24$ for $Q=0$. The diagonal of $W$ (mass and electric energy) is computed once per size as three vectors, since $\sum_nL_n^2=\sum_nc_n^2+2\varepsilon_0\sum_nc_n+(N-1)\varepsilon_0^2$ with $c_n=L_n-\varepsilon_0$, so that changing $m$, $x$ or $\varepsilon_0$ costs nothing.

**Matrix–vector product.** The hopping term connects each basis state to at most $N-1$ others. Their indices are tabulated once (binary search in the sorted basis) and the product is evaluated by an OpenMP-parallel C kernel (`kernel.c`, 12 lines).

**Eigensolver.** Lanczos iteration with full (twice-applied) Gram–Schmidt reorthogonalisation against all previous vectors. Full reorthogonalisation removes ghost eigenvalues, which matters here because we need several low-lying levels and exact degeneracies are part of the prediction. Iteration stops when the lowest $k$ Ritz values change by less than $10^{-10}$ (relative) over ten steps. Memory for the Lanczos basis is the binding constraint and limits us to $N\le24$.

**Validation.** The code is checked (`test_qed2.py`) against an independent dense Hamiltonian built from Kronecker products of Pauli matrices in the full $2^N$-dimensional space and then restricted to $Q=0$; the four lowest levels agree to $10^{-9}$ for three parameter sets including $\varepsilon_0\ne0$.

**Infinite-volume gap away from criticality.** In a massive phase the first excited state is a single boson in a box, with $\Delta^2=M^2+\pi^2v^2/(L+\ell)^2$. We fit $\Delta^2=M^2+C/(N+n_0)^2$ to sizes $N=12,\ldots,22$ with $M$, $C$, $n_0$ free.

**String tension.** $\sigma$ is the difference of bulk energy densities with and without background field. The bulk density is taken from $[E_0(N_2)-E_0(N_1)]/(L_2-L_1)$, which cancels boundary energies exactly and leaves corrections of order $e^{-ML}$. Because $\sigma\equiv0$ at $m=0$, the measured value there is a direct diagnostic of the residual finite-volume error; we keep only volumes $Lg\ge11$, for which it is below $10^{-3}g^2$.

**Critical point.** The ratio $R=(E_2-E_0)/(E_1-E_0)$ is computed on a grid of six masses spaced by $0.01g$ and interpolated by a cubic. The pseudo-critical mass $m_c(N)$ is defined as the crossing of $R(m;N-2)$ and $R(m;N)$; near $m_c$ the ratio depends on $m$ and $N$ only through $(m-m_c)N^{1/\nu}$, so curves for different sizes cross at $m_c$ up to corrections from irrelevant operators. The central charge then follows from Eq. (4.6) by taking second differences in $N$, which remove the terms linear in $L$:

$$c=-48\,\frac{E_0(N)-2E_0(N-2)+E_0(N-4)}{\Delta(N)-2\Delta(N-2)+\Delta(N-4)},\qquad v=\left[2\pi\,\frac{1/\Delta(N)-1/\Delta(N-2)}{L(N)-L(N-2)}\right]^{-1}, \qquad (6.1)$$

with $\Delta=E_1-E_0$. Both are evaluated at $m_c(N)$.

All scripts and raw numbers are in the repository; Appendix B lists them.

# Results

## Mass gap of the massless model (P1, P2)

Figure 1(a) shows the gap at $m=0$ for three lattice spacings. Plotted as $\Delta^2$ against $1/L^2$ the data collapse onto a straight line with intercept $1/\pi$, as expected for a free boson in a box. Table 1 and Fig. 1(b) give the infinite-volume gap. For $x\ge2$ it approaches the exact value monotonically from below; a linear extrapolation of the three finest spacings gives

$$\frac{M_0}{g}=0.5636\;\text{(linear in }a^2)\qquad\text{or}\qquad0.5641\;\text{(linear in }a),\qquad\text{exact: }0.56419.$$

The two extrapolations bracket the residual uncertainty; either agrees with Eq. (2.8) to $0.1\%$. The value at $x=1$ lies slightly above the trend, an indication that at $ga=1$ higher-order lattice artefacts are not negligible.

| $x$ | 1 | 2 | 4 | 6 | 9 | 12 | 16 |
|---|---|---|---|---|---|---|---|
| $ga$ | 1.00 | 0.71 | 0.50 | 0.41 | 0.33 | 0.29 | 0.25 |
| $Lg$ at $N=22$ | 22.0 | 15.6 | 11.0 | 9.0 | 7.3 | 6.4 | 5.5 |
| $M_0/g$ | 0.56463 | 0.56147 | 0.56208 | 0.56260 | 0.56287 | 0.56303 | 0.56318 |

: Infinite-volume gap of the massless model.

![Mass gap at $m=0$. (a) Finite-volume gap squared against inverse length squared; the dashed line is the exact $M_0^2=g^2/\pi$. (b) Infinite-volume gap against lattice spacing squared.](fig/fig1_gap.svg){width=100%}

For the slope (P2) we repeat the fit at $m/g=0.02$ and $0.04$ and extrapolate the difference quotient to $m=0$: $dM/dm=1.753$ ($x=2$), $1.798$ ($x=4$) and $1.874$ ($x=9$), against $e^\gamma=1.781$. The first two agree to $1$–$2\%$. The third is $5\%$ high; at $x=9$ the largest volume is only $Lg=6.7$ and we cannot separate lattice from finite-volume error there, so we do not regard it as a meaningful discrepancy, but neither is it a confirmation.

## Screening and confinement (P3–P5)

Table 2 and Fig. 2 show the string tension. At $m=0$ the measured $\sigma$ is zero within the finite-volume diagnostic: $|\sigma|/g^2\le1\times10^{-4}$ for $\alpha=\tfrac14$ at $x=2,3$ and below $10^{-5}$ for $\alpha=\tfrac12$. Charge is screened.

For $m>0$ the string tension is non-zero and, for $\alpha=\tfrac12$, linear in $m$ with slope within $2\%$ of $2c$. For $\alpha=\tfrac14$ the first-order formula overestimates $\sigma$ by $12\%$ at $m/g=0.05$ and $27\%$ at $m/g=0.1$. The second-order term of Eq. (2.15) accounts for this: with it, lattice and theory agree to between $0.4\%$ and $1.4\%$ for $m/g\le0.1$ at $x=3$. A straight-line fit of $\sigma/m$ against $m$ gives an independent lattice estimate of both coefficients (last rows of Table 2).

| $m/g$ | $\alpha=\tfrac14$ lattice | theory, 1st | theory, 1st+2nd | $\alpha=\tfrac12$ lattice | theory |
|---|---|---|---|---|---|
| 0.025 | 0.00376 | 0.00400 | 0.00377 | 0.00815 | 0.00800 |
| 0.050 | 0.00713 | 0.00800 | 0.00710 | 0.01628 | 0.01599 |
| 0.075 | 0.01006 | 0.01199 | 0.00998 | 0.02437 | 0.02399 |
| 0.100 | 0.01259 | 0.01599 | 0.01241 | 0.03241 | 0.03199 |
| 0.150 | 0.01662 | 0.02399 | 0.01593 | 0.04818 | 0.04798 |
| linear coeff. | 0.1588 | 0.1599 | | 0.3268 | 0.3199 |
| quadratic coeff. | $-0.328$ | | $-0.358$ | $-0.025$ | 0 |

: String tension $\sigma/g^2$ at $x=3$ ($ga=0.58$), from $N=18$ and $22$ ($Lg=10.4$, $12.7$). Fitted coefficients use $m/g\le0.1$.

![String tension divided by $mg$. Dashed: first-order mass perturbation theory. Solid: including the second-order term of Eq. (2.15). The lattice data approach the solid line as the lattice spacing decreases.](fig/fig2_string.svg){width=100%}

At $x=2$ the fitted coefficients are $0.1661$ and $-0.394$ ($\alpha=\tfrac14$) and at $x=1$ they are $0.187$ and $-0.57$, so the approach to the continuum values is from one side and consistent with $O(a^2)$ artefacts, but three spacings at different volumes do not justify a formal extrapolation.

## The critical point at θ = π (P6, P7, P11, P12)

Figure 4(a) shows the level ratio $R$ against $m$ for $x=0.5$. The curves for successive sizes cross, and the crossing point converges rapidly. Table 3 lists the crossing mass, the ratio there, and $c$ and $v$ from Eq. (6.1), for three lattice spacings. We deliberately work at coarse spacing, $ga=2$, $1.41$, $1$: the universality class does not depend on $a$, and coarse lattices give the largest physical volume ($Lg=44$, $31$–$34$, $22$ respectively), which is what the scaling regime requires. At $x\ge2$ no crossing exists for $N\le18$.

| $N$ | $x{=}0.25$: $m_c$ | $R$ | $c$ | $v$ | $x{=}0.5$: $m_c$ | $R$ | $c$ | $v$ | $x{=}1$: $m_c$ | $R$ | $c$ | $v$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 12 | 0.4047 | 1.4282 | 0.317 | 1.132 | 0.3623 | 1.4755 | 0.447 | 1.011 | 0.3223 | 1.5523 | 0.781 | 0.859 |
| 14 | 0.4009 | 1.4420 | 0.349 | 1.090 | 0.3601 | 1.4812 | 0.479 | 0.996 | 0.3244 | 1.5487 | 0.822 | 0.863 |
| 16 | 0.3987 | 1.4513 | 0.370 | 1.064 | 0.3593 | 1.4836 | 0.489 | 0.990 | 0.3280 | 1.5417 | 0.800 | 0.875 |
| 18 | 0.3972 | 1.4581 | 0.386 | 1.046 | 0.3590 | 1.4846 | 0.490 | 0.989 | 0.3314 | 1.5343 | 0.760 | 0.890 |
| 20 | 0.3962 | 1.4633 | 0.400 | 1.033 | 0.3589 | 1.4851 | 0.488 | 0.989 | 0.3341 | 1.5276 | 0.718 | 0.904 |
| 22 | 0.3955 | 1.4675 | 0.413 | 1.022 | 0.3588 | 1.4855 | 0.485 | 0.989 | 0.3361 | 1.5218 | 0.681 | 0.917 |
| 24 | | | | | 0.3587 | 1.4859 | 0.483 | 0.989 | | | | |
| $\infty$, $1/N^2$ | 0.392 | 1.486 | 0.46 | 0.98 | 0.358 | 1.488 | 0.47 | 0.988 | 0.345 | 1.500 | 0.55 | 0.96 |
| $\infty$, $1/N$ | 0.387 | 1.511 | 0.52 | 0.91 | 0.358 | 1.490 | 0.46 | 0.988 | 0.358 | 1.469 | 0.37 | 1.03 |

: Pseudo-critical mass $m_c/g$ (crossing of $R$ for sizes $N-2$ and $N$) and the ratio $R$, central charge $c$ and velocity $v$ evaluated there. The last two rows are linear extrapolations of the four largest sizes in $1/N^2$ and in $1/N$; Section 4 argues that $1/N^2$ is the appropriate form for $R$.

![(a) Level ratio $R$ against mass at $ga=1.41$ for increasing chain length (light to dark). The dashed line is the BCFT value $\tfrac32$. (b) Value of $R$ at the crossing point against $1/N^2$ for three lattice spacings.](fig/fig4_crossing.svg){width=100%}

**Existence of the transition (P6).** At the crossing the gap falls as $1/(L+\ell)$: the velocity estimator of Eq. (6.1), which is the inverse slope of $1/\Delta$ against $L$, is constant to three digits over $N=16$–$22$ at $x=0.5$. Away from the crossing, on either side, the gap saturates.

**Ratio (P7).** At the three spacings $R$ approaches $\tfrac32$ from below, from below, and from above (Fig. 4b). The $1/N^2$ extrapolations are $1.486$, $1.488$ and $1.500$. At $x=0.5$ the sequence has essentially converged and sits $0.8\%$ below $\tfrac32$; we return to this in Section 8.

**Velocity (P12).** $v=0.989$ at $x=0.5$, and $0.98$, $0.96$ ($1/N^2$ extrapolations) at $x=0.25$, $1$. The velocity is not protected on the lattice and is $1$ only in the continuum limit; that it is within a few per cent of $1$ at $ga\sim1$–$2$ is a consequence of the improved mass (5.5) and is not a test of the continuum prediction.

**Exponent $\nu$ (P11).** The slope $dR/dm$ at the crossing should grow as $N^{1/\nu}$. The effective exponent from successive sizes at $N=22$ is $1/\nu=0.87$, $0.89$, $1.06$ at $x=0.25$, $0.5$, $1$. Replacing $N$ by $N+n_0$ with the extrapolation length $n_0=3.8$ determined below raises the $x=0.5$ value to $1.05$. We conclude $\nu=1.0\pm0.1$.

## The conformal towers (P7, P8)

The crossing analysis uses only two excited levels. A much more stringent test is the full low-lying spectrum, including degeneracies. Table 4 and Fig. 5 give the first five excitation energies in units of the first, at $x=0.5$ and fixed $m/g=0.356$.

| $N$ | $Q=0$: $n=2$ | 3 | 4 | 5 | $Q=-1$: $n=2$ | 3 | 4 | 5 |
|---|---|---|---|---|---|---|---|---|
| 12 | 1.4919 | 1.9434 | 1.9567 | 2.3422 | 1.9271 | 2.7315 | 3.3616 | 3.7106 |
| 14 | 1.4934 | 1.9580 | 1.9671 | 2.3848 | 1.9430 | 2.7956 | 3.5176 | 3.7637 |
| 16 | 1.4947 | 1.9683 | 1.9743 | 2.4125 | 1.9531 | 2.8367 | 3.6196 | 3.8006 |
| 18 | 1.4958 | 1.9758 | 1.9799 | 2.4316 | 1.9598 | 2.8644 | 3.6892 | 3.8268 |
| 20 | 1.4968 | 1.9814 | 1.9846 | 2.4454 | 1.9646 | 2.8839 | 3.7386 | 3.8457 |
| **BCFT** | $\mathbf{3/2}$ | $\mathbf2$ | $\mathbf2$ | $\mathbf{5/2}$ | $\mathbf2$ | $\mathbf3$ | $\mathbf4$ | $\mathbf4$ |

: Excitation energies $(E_n-E_0)/(E_1-E_0)$ within the $Q=0$ and $Q=-1$ sectors at $x=0.5$, $m/g=0.356$, $\theta=\pi$.

![Excitation energies in units of the lowest, against $1/N^2$. Dashed lines are the BCFT values. (a) Neutral sector: the identity tower, with the two-fold degenerate level at $h=4$. (b) Charge $-1$: the $\varepsilon$ tower, with the two-fold degenerate level at $h=\tfrac92$.](fig/fig5_towers.svg){width=100%}

In the neutral sector the four ratios are within $0.2\%$, $0.9\%$, $0.8\%$ and $2.2\%$ of the identity-tower values at $N=20$ and all are still moving toward them. The two states predicted to be degenerate at $h=4$ ($L_{-4}|0\rangle$ and $L_{-2}^2|0\rangle$) appear as a pair split by $0.16\%$ of their energy, with the splitting decreasing monotonically (from $0.68\%$ at $N=12$). In the charge $-1$ sector the ratios are within $1.8\%$, $3.9\%$, $6.5\%$ and $3.9\%$ of $2,3,4,4$, converge more slowly, and again the degenerate pair is forming (splitting $2.8\%$ at $N=20$, down from $9.9\%$ at $N=12$). The $Q=+1$ sector reproduces the identity pattern ($1.485$, $1.948$, $1.964$, $2.378$ at $N=20$) at an energy offset that grows linearly with $N$ in units of the gap, i.e. is a constant in physical units, as expected for a charge bound to the wall.

The alternative boundary conditions are excluded by a wide margin: free boundaries would put the second neutral level at ratio $3$, not $\tfrac32$.

## Central charge (P9)

The Casimir estimator (6.1) gives the values in Table 3 and Fig. 6(a). At $x=0.5$ it has a broad plateau, $c=0.489$, $0.490$, $0.488$, $0.485$, $0.483$ for $N=16$–$24$, within $2$–$4\%$ of $\tfrac12$ and drifting slowly downward. The estimator is, however, very sensitive to the mass at which it is evaluated: at fixed $m/g=0.357$, only $0.0017$ below the crossing, it gives $c=0.542$–$0.545$, constant over the same sizes, so that $dc/dm\approx-30/g$. An uncertainty of $\pm0.001g$ in the critical mass therefore translates into $\pm0.03$ in $c$, and our result is $c=0.50(4)$. At the other two spacings it is still drifting — upward through $0.41$ at $x=0.25$, downward through $0.68$ at $x=1$ — so that the three sequences bracket $\tfrac12$ but only one of them can be said to have converged. This is the first point at which the limitations of 24 sites are clearly visible: the estimator is a ratio of second differences and amplifies corrections to scaling.

![(a) Central charge from the Casimir energy, Eq. (6.1), and (b) velocity, against $1/N$, at the crossing point.](fig/fig6_c_v.svg){width=100%}

## Order parameter and entanglement (P10)

At $x=0.5$, $m/g=0.3585$ we computed the ground-state wavefunction for $N=8,\ldots,22$. The electric field alternates strongly between even and odd links (a staggering artefact: $0.464$ and $0.152$ at the centre for $N=22$), so we average the two central links. A fit of $1/\Delta$ against $N$ gives the extrapolation length $n_0=\ell/a=3.78$ sites. With this value — determined from the spectrum, not from the order parameter — a log–log fit of the mid-chain field against $N+n_0$ over $N=14$–$22$ gives

$$\beta/\nu=0.1263\qquad(\text{Ising: }0.125).$$

Without the shift ($n_0=0$) the same fit gives $0.104$, so the agreement relies on the BCFT form (4.7); it is a consistency check of the whole picture rather than an unbiased measurement of the exponent (Fig. 7a).

The entanglement entropy (Fig. 7b) is dominated at these sizes by an oscillation between even and odd cuts of amplitude $\approx0.4$, larger than the total variation of the logarithm in Eq. (4.8). Fitting each parity separately gives $c=0.43$ (even) and $0.26$ (odd) at $N=22$; averaging adjacent cuts gives $0.36$. All three increase slowly with $N$. **We cannot extract $c$ from the entanglement entropy at $N\le22$**, and report this as a failure of the method at these sizes rather than of the prediction.

![(a) Mid-chain electric field against $N+n_0$ on logarithmic axes; the dashed line has slope $-\tfrac18$. (b) Entanglement entropy against the logarithm of the chord length at $N=22$; the dashed line has the slope $c/6$ for $c=\tfrac12$.](fig/fig7_order_ent.svg){width=100%}

## The critical mass

Table 3 gives $m_c/g\approx0.392$, $0.358$ and $0.345$ at $ga=2$, $1.41$ and $1$. These are properties of the lattice model at fixed spacing, they are well determined (the $x=0.5$ value is stable to four digits), and they decrease toward the continuum value $m_c/g=0.333561(4)$ of Ref. [@act] (earlier: $0.3335(2)$ [@byrnes], $0.325(20)$ [@hamer82]). We have not attempted a continuum extrapolation: three spacings with $ga\ge1$ are not in the regime where $m_c(a)=m_c+O(a^2)$ can be trusted, and finer lattices need $N\gtrsim30$ to reach the scaling volume.

Against this, the analytic estimates of Section 3.2 fare poorly. The classical value $0.158$ is a factor of two too small. The Gaussian value $0.265$ is closer but belongs to a first-order transition: the Gaussian state predicts a finite gap and a jump in $\langle E\rangle$ at the transition, whereas the lattice shows a gap closing as $1/L$ and a conformal spectrum. The failure is qualitative. Its origin is that a Gaussian state cannot represent the critical Ising fixed point, which is a free Majorana fermion and not a free boson.

# Discussion

## What is established

Table 5 collects the outcome.

| | Prediction | Theory | Lattice | Status |
|---|---|---|---|---|
| P1 | $M_0/g$ | 0.56419 | 0.5636–0.5641 | confirmed, 0.1% |
| P2 | $dM/dm$ | 1.781 | 1.753, 1.798 (1.874 at small volume) | confirmed, 1–2%, where volume is adequate |
| P3 | $\sigma(m=0)$ | 0 | $<10^{-4}g^2$ | confirmed |
| P4 | linear coeff. of $\sigma$, $\alpha=\tfrac14$; $\tfrac12$ | 0.1599; 0.3199 | 0.1588; 0.3268 | confirmed, 1–2% |
| P5 | quadratic coeff., $\alpha=\tfrac14$; $\tfrac12$ | $-0.358$; 0 | $-0.328$; $-0.025$ | confirmed, ~10% |
| P6 | gapless point at $\theta=\pi$ | yes | yes | confirmed |
| P7 | $Q=0$ ratios | 1.5, 2, 2, 2.5 | 1.497, 1.981, 1.985, 2.445 | confirmed, 0.2–2.2% |
| P8 | $Q=-1$ ratios | 2, 3, 4, 4 | 1.965, 2.884, 3.739, 3.846 | confirmed, 2–7%, converging |
| P9 | $c$ (Casimir) | 0.5 | 0.48–0.49 at $ga=1.41$; 0.41↑, 0.68↓ elsewhere | supported |
| P10 | $\beta/\nu$ | 0.125 | 0.126 (with $n_0$ from gap); 0.104 (raw) | consistent |
| P11 | $\nu$ | 1 | $1.0\pm0.1$ | consistent |
| P12 | $v$ | 1 (continuum) | 0.96–0.99 at $ga=1$–2 | not tested in the continuum |

: Summary of tests.

## Open points and sources of error

**The ratio at the crossing, $x=0.5$.** The crossing value of $R$ is $1.4857$–$1.4859$ for every size from $N=16$ to $24$ — a fixed point to four digits — but it sits $0.9\%$ below $\tfrac32$, while the same ratio at fixed $m/g=0.356$ is $1.4968$ and rising, and at $m/g=0.355$ it passes through $\tfrac32$ between $N=18$ and $20$ and keeps rising. Since $dR/dm\approx-4$ at $N=22$, a shift of the true critical mass by $-0.003g$ relative to the crossing would account for the difference. Crossing points converge to $m_c$ as $N^{-1/\nu-\omega}$ where $\omega$ is the leading correction exponent, so a small $\omega$, or two corrections of opposite sign that nearly cancel in the crossing but not in $R$, would produce exactly this pattern. Ref. [@act], working with chains of several hundred sites, fits the finite-size dependence of the same criterion with $N^{-3}$ and $N^{-4}$ terms of opposite sign and large coefficients, which supports the view that sizes $N\le24$ are not asymptotic and that the approach is not monotonic. We cannot resolve this with exact diagonalisation.

**Central charge.** Only one of three lattice spacings has converged. The estimator is also steeply mass-dependent (Section 7.5). The honest statement is $c=0.50(4)$ at $ga=1.41$ and "consistent with $\tfrac12$" elsewhere.

**Continuum limit.** Only the $\theta=0$ results (P1–P5) are at lattice spacings fine enough to speak of a continuum limit, and even there the string tension is limited to $ga\ge0.58$ by the volume requirement. All $\theta=\pi$ results are universal quantities measured at coarse spacing, which is legitimate for universal quantities and not for $m_c$.

**Systematic versus statistical error.** Exact diagonalisation has no statistical error; eigenvalues are converged to $10^{-10}$. Every uncertainty quoted is systematic — finite volume, finite lattice spacing, and the choice of extrapolation form — and we have indicated its size by quoting alternative extrapolations rather than a single error bar.

**Open boundaries and $\theta$.** With open boundaries the $2\pi$ periodicity in $\theta$ is lost in each fixed-charge sector and restored only by comparing sectors; this is why $Q=0$ and $Q=-1$ realise different boundary conditions at the same $\theta$.

## Relation to the four-dimensional problem

It is tempting to read Section 2.2 as a template: an anomaly converts a would-be massless gauge boson into a massive one. The template does not transfer. In four dimensions the axial current is not dual to the vector current, the gluon has propagating transverse modes, the anomaly equation does not close on a single field, and the gap is believed to arise from the dynamics of the non-abelian field itself, with no fermions required. What does transfer is methodological: the Hamiltonian lattice formulation with Gauss's law imposed exactly, the use of $\theta$ and external charges as probes of confinement, and the practice of testing a proposed mechanism through parameter-free finite-size predictions. The Schwinger model is a laboratory in which that practice can be carried to completion; it is not evidence about Yang–Mills theory.

# Conclusion

Starting from the Lagrangian of two-dimensional QED we derived the mass gap from the anomaly, the bosonized Hamiltonian, the screening of charge at $m=0$ and its confinement at $m\ne0$, the $Z_2$ transition at $\theta=\pi$, and — from the observation that the background field pins both walls of an open system to the same vacuum — the complete finite-size spectrum at criticality. Exact diagonalisation of the lattice gauge theory on up to 24 sites confirms the mass gap to $0.1\%$, the string tension through second order in the fermion mass to $1$–$2\%$, and the Ising boundary-CFT spectrum in two charge sectors including its degeneracies. The Casimir central charge and the exponents $\nu$ and $\beta$ are consistent with Ising values at the few-per-cent to ten-per-cent level. Simple analytic estimates of the critical mass fail, the Gaussian variational method qualitatively so. The continuum critical mass is beyond the method, and has in any case been determined to five digits by tensor networks [@act].

# Appendix A: Derivation of the Gaussian effective potential {.unnumbered}

Let $|\mu,\phi_0\rangle$ be the vacuum of a free field of mass $\mu$ displaced so that $\langle\phi\rangle=\phi_0$. For operators normal-ordered with respect to $\mu$, expectation values in this state are obtained by replacing $\phi\to\phi_0$. The Hamiltonian density (2.11) at $\theta=\pi$, after the shift, is $N_{M}[\tfrac12\Pi^2+\tfrac12\phi'^2+\tfrac12M^2\phi^2+\kappa\cos\beta\phi]$ with $\kappa=c\,mg$, $\beta=2\sqrt\pi$, $M=M_0$. The equal-time two-point functions of free fields of mass $M$ and $\mu$ differ by

$$\langle\phi^2\rangle_\mu-\langle\phi^2\rangle_M=\int\frac{dk}{4\pi}\left(\frac{1}{\sqrt{k^2+\mu^2}}-\frac{1}{\sqrt{k^2+M^2}}\right)=\frac{1}{4\pi}\ln\frac{M^2}{\mu^2},$$

$$\left\langle\tfrac12\Pi^2+\tfrac12\phi'^2\right\rangle_\mu-\langle\ldots\rangle_M=\int\frac{dk}{8\pi}\left(\frac{2k^2+\mu^2}{\sqrt{k^2+\mu^2}}-\frac{2k^2+M^2}{\sqrt{k^2+M^2}}\right)=\frac{\mu^2-M^2}{8\pi},$$

which are the first two identities of Eq. (3.3). For the exponential, $N_M[e^{i\beta\phi}]=e^{\frac{\beta^2}{2}\langle\phi^2\rangle_M}e^{i\beta\phi}$, so $N_M[e^{i\beta\phi}]=\exp\{\tfrac{\beta^2}{2}(\langle\phi^2\rangle_M-\langle\phi^2\rangle_\mu)\}N_\mu[e^{i\beta\phi}]=(\mu/M)^{\beta^2/4\pi}N_\mu[e^{i\beta\phi}]$. Collecting terms,

$$\mathcal E(\mu,\phi_0)=\frac{\mu^2-M^2}{8\pi}+\frac{M^2}{8\pi}\ln\frac{M^2}{\mu^2}+\frac{M^2}{2}\phi_0^2+\kappa\,\frac{\mu}{M}\cos2\sqrt\pi\phi_0,$$

which is Eq. (3.4) after dividing by $M^2$, with $k=\kappa/M^2=\pi c\,m/g$. At $\kappa=0$ the minimum is at $\mu=M$, $\phi_0=0$, as it should be. The exponent $\beta^2/4\pi=1$ is the scaling dimension of the fermion mass operator; it is this linear dependence on $\mu$ that prevents the curvature and gap conditions from being satisfied together.

# Appendix B: Reproducibility {.unnumbered}

| Script | Produces | Runtime (4 cores) |
|---|---|---|
| `qed2.py`, `kernel.c` | Hamiltonian and Lanczos solver | — |
| `test_qed2.py` | comparison with dense Hamiltonian | seconds |
| `run_gap.py` | Table 1 | 6 min |
| `run_slope.py` | $dM/dm$ | 3 min |
| `run_string.py`, `analyze_string.py` | Table 2 | 18 min |
| `run_crit.py`, `analyze_crit.py` | Table 3 | 45 min ($N\le22$); 1.5 h more for $N=24$ |
| `run_tower.py`, `run_sectors.py` | Table 4 | 12 min |
| `run_groundstate.py`, `analyze_extra.py` | Sec. 7.5, $\nu$ | 5 min |
| `check_formulas.py` | characters (4.3), $\mathcal A_+$, identities (3.3) | seconds |
| `make_figures.py` | all figures | seconds |
| `paper/build.py` | this document | seconds |

Raw results are stored as JSON in `data/`. The kernel is built with `gcc -O3 -march=native -fopenmp -shared -fPIC kernel.c -o kernel.so`.

# References {.unnumbered}

Bibliographic details of Refs. [@adam] and [@act] were checked against the arXiv records during this work, and the equations quoted from them against the text. The remaining entries are standard and are cited from memory; volume and page numbers should be verified before formal submission.

::: {#refs}
:::
