"""Third-order mass perturbation theory for the vacuum energy of the Schwinger model.

With S = (Sigma/2) (e^{i theta} N[e^{+i b phi}] + e^{-i theta} N[e^{-i b phi}]), b = 2 sqrt(pi), and
<prod_j N[e^{i s_j b phi(x_j)}]>_0 = prod_{j<k} exp(-2 s_j s_k K0(mu0 r_jk)), the third cumulant gives

  eps_3(theta) = (m^3 Sigma^3 / (48 mu0^4)) [ 2 cos(3 theta) I3 + 6 cos(theta) I1 ],

  I3 = int d2x d2y [ g(x) g(y) g(x-y) - g(x) - g(y) - g(x-y) + 2 ],        g = exp(-2 K0)
  I1 = int d2x d2y [ g(x) h(y) h(x-y) - g(x) - h(y) - h(x-y) + 2 ],        h = exp(+2 K0)

(lengths in units of 1/mu0; Z = exp(V eps), energy density = -eps).  The second-order analogue,
E_+ = int d2x (g - 1) = -8.9139, is reproduced as a check of the machinery.
"""
import sys
import numpy as np

EG = np.euler_gamma
_t = np.linspace(0, 14, 2801)
_r = np.exp(np.linspace(np.log(1e-9), np.log(60.0), 6000))
_K = np.empty_like(_r)
for i0 in range(0, len(_r), 500):
    rr = _r[i0:i0 + 500]
    _K[i0:i0 + 500] = np.trapezoid(np.exp(-np.outer(rr, np.cosh(_t))), _t, axis=1)
small = _r < 1e-3                                   # series for small argument
_K[small] = -(np.log(_r[small] / 2) + EG) * (1 + _r[small] ** 2 / 4) + _r[small] ** 2 / 4


def K0(r):
    return np.interp(np.log(np.clip(r, 1e-9, 60.0)), np.log(_r), _K)


def gl(n):
    return np.polynomial.legendre.leggauss(n)


def integrals(nX=260, nR=320, nphi=160, Rmax=30.0, rmin=1e-6):
    """I3 and I1.  x = (X, 0); y in polar coordinates about the origin, restricted to the half-plane
    y.x < X^2/2 (|y| < |x - y|), which by the symmetry y <-> x - y is half of the integral."""
    uX, wX = gl(nX); uR, wR = gl(nR); uP, wP = gl(nphi)
    lX = 0.5 * (uX + 1) * (np.log(Rmax) - np.log(1e-4)) + np.log(1e-4); X = np.exp(lX)
    WX = 0.5 * (np.log(Rmax) - np.log(1e-4)) * wX * X                     # dX
    I3 = I1 = 0.0
    for X_, wx in zip(X, WX):
        gx = np.exp(-2 * K0(X_))
        # the angular range changes character at rho = X/2: integrate rho in two panels split there
        cuts = [rmin, X_ / 2, Rmax]
        rho = np.concatenate([np.exp(0.5 * (uR + 1) * (np.log(b_) - np.log(a_)) + np.log(a_)) for a_, b_ in zip(cuts[:-1], cuts[1:])])
        WR = np.concatenate([0.5 * (np.log(b_) - np.log(a_)) * wR for a_, b_ in zip(cuts[:-1], cuts[1:])]) * rho
        # allowed angles: cos(phi) < X/(2 rho)
        c = np.clip(X_ / (2 * rho), -1, 1); phi0 = np.where(X_ / (2 * rho) >= 1, 0.0, np.arccos(c))
        half = np.pi - phi0                                               # phi in [phi0, 2pi - phi0]
        phi = np.pi + np.outer(half, uP)                                  # (nR, nphi)
        wphi = np.outer(half, wP)
        d = np.sqrt(np.maximum(rho[:, None] ** 2 + X_ ** 2 - 2 * rho[:, None] * X_ * np.cos(phi), 1e-300))
        Ky = K0(rho)[:, None]; Kd = K0(d)
        gy, gd = np.exp(-2 * Ky), np.exp(-2 * Kd)
        hy, hd = np.exp(2 * Ky), np.exp(2 * Kd)
        f3 = gx * gy * gd - gx - gy - gd + 2
        # written to keep the cancellation at y -> 0 explicit: h(y) [g(x) h(x-y) - 1] - (h(x-y) - 1) - (g(x) - 1)
        f1 = hy * (gx * hd - 1.0) - (hd - 1.0) - (gx - 1.0)
        I3 += wx * 2 * np.pi * X_ * 2 * np.sum(WR * rho * np.sum(wphi * f3, axis=1))
        I1 += wx * 2 * np.pi * X_ * 2 * np.sum(WR * rho * np.sum(wphi * f1, axis=1))
    return I3, I1


if __name__ == "__main__":
    r = np.exp(np.linspace(np.log(1e-7), np.log(50), 20001))
    Ep = 2 * np.pi * np.trapezoid(r * r * (np.exp(-2 * K0(r)) - 1), np.log(r))
    print(f"check: E_+ = {Ep:.4f}  (Adam: -8.9139)")
    res = {}
    for n in ((100, 80, 64), (140, 110, 96), (200, 160, 128), (280, 220, 180)):
        I3, I1 = integrals(*n)
        res[n] = (I3, I1); print(f"grid {n}: I3 = {I3:.4f}   I1 = {I1:.4f}", flush=True)
    # the grid sequence converges roughly geometrically; extrapolate with the last three members
    ks = list(res)[-3:]
    def aitken(a0, a1, a2):
        d = (a2 - a1) - (a1 - a0)
        return a2 if abs(d) < 1e-12 else a2 - (a2 - a1) ** 2 / d
    I3 = aitken(*[res[k][0] for k in ks]); I1 = aitken(*[res[k][1] for k in ks])
    e3 = abs(I3 - res[ks[-1]][0]); e1 = abs(I1 - res[ks[-1]][1])
    print(f"extrapolated: I3 = {I3:.3f} +- {e3:.3f}   I1 = {I1:.3f} +- {e1:.3f}")
    S = np.exp(EG) / (2 * np.pi)                       # Sigma / mu0
    pref = S ** 3 * np.sqrt(np.pi) / 48                # m^3 Sigma^3/(48 mu0^4) in units g^2 (m/g)^3
    print(f"prefactor = {pref:.6e}")
    for al in (0.25, 0.5):
        c3 = -pref * (2 * I3 * (np.cos(6 * np.pi * al) - 1) + 6 * I1 * (np.cos(2 * np.pi * al) - 1))
        dc = pref * np.hypot(2 * e3 * (np.cos(6 * np.pi * al) - 1), 6 * e1 * (np.cos(2 * np.pi * al) - 1))
        print(f"alpha = {al}: predicted third-order coefficient c3 = {c3:+.4f} +- {dc:.4f}")
    import json
    json.dump({"I3": I3, "I1": I1, "I3_err": e3, "I1_err": e1, "prefactor": pref, "E_plus_check": Ep,
               "grids": {str(k): v for k, v in res.items()}}, open("data/third_order.json", "w"), indent=1)
