"""Independent numerical checks of formulas quoted in the paper."""
import numpy as np

# 1. Ising characters from the free-fermion product formulas (half-integer powers of q)
K = 14                                     # work in powers of q^(1/2) up to q^7
def prod(sign):
    p = np.zeros(2 * K + 1); p[0] = 1
    for n in range(1, 2 * K + 2, 2):       # exponents n/2, n odd
        new = p.copy(); new[n:] += sign * p[:-n] if n <= 2 * K else 0
        p = new
    return p
pp, pm = prod(+1), prod(-1)
chi0 = ((pp + pm) / 2)[0::2]; chi12 = ((pp - pm) / 2)[1::2]
print("chi_0   coefficients of q^0..q^6 :", chi0[:7].astype(int))
print("chi_1/2 coefficients of q^(1/2+n):", chi12[:7].astype(int))

# 2. Adam's constant A_+ = 2 pi int_0^inf r (exp(-2 K0(r)) - 1) dr
t = np.linspace(0, 12, 4001)
def K0(r):
    return np.trapezoid(np.exp(-np.outer(r, np.cosh(t))), t, axis=1)
r = np.concatenate([np.linspace(1e-6, 1, 4000), np.linspace(1, 40, 8000)[1:]])
A_plus = 2 * np.pi * np.trapezoid(r * (np.exp(-2 * K0(r)) - 1), r)
print(f"A_+ = {A_plus:.4f}   (paper: -8.9139)")
S = np.exp(np.euler_gamma) / (2 * np.pi)
print(f"b(1/4) = {S ** 2 / 4 * A_plus * 2:.4f}   (paper: -0.3582)")

# 3. Re-ordering identities used in the Gaussian effective potential
k = np.linspace(0, 4000, 4_000_001); mu, M = 0.6, 1.0
kin = np.trapezoid((2 * k**2 + mu**2) / np.sqrt(k**2 + mu**2) - (2 * k**2 + M**2) / np.sqrt(k**2 + M**2), k) * 2 / (8 * np.pi)
phi2 = np.trapezoid(1 / np.sqrt(k**2 + mu**2) - 1 / np.sqrt(k**2 + M**2), k) * 2 / (4 * np.pi)
print(f"kinetic: {kin:.5f} vs (mu^2-M^2)/8pi = {(mu**2 - M**2) / (8 * np.pi):.5f}")
print(f"phi^2  : {phi2:.5f} vs ln(M^2/mu^2)/4pi = {np.log(M**2 / mu**2) / (4 * np.pi):.5f}")

# 4. constants
c = np.exp(np.euler_gamma) / (2 * np.pi ** 1.5)
print(f"c = {c:.5f}; 1/sqrt(pi) = {1 / np.sqrt(np.pi):.5f}; e^gamma = {np.exp(np.euler_gamma):.4f}; "
      f"classical m_c = {1 / (4 * np.pi**2 * c):.4f}")
