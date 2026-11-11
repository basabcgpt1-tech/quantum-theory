"""Error budget for the continuum critical mass and the other headline numbers.
Each source of systematic error is varied separately and the resulting shift recorded."""
import io, contextlib, json
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    import analyze_dmrg as A

cont = A.summary["cont"]
u = np.array([1 / c["x"] for c in cont]); y0 = np.array([c["mc"] for c in cont]); e0 = np.array([c["err"] for c in cont])


def fit(u, y, e, sel, deg):
    M = np.vander(u[sel], deg + 1) / e[sel, None]
    cov = np.linalg.inv(M.T @ M); p = cov @ M.T @ (y[sel] / e[sel])
    return p[-1], np.sqrt(cov[-1, -1])


central, stat = fit(u, y0, e0, u <= 0.5, 2)                     # x >= 2, quadratic in a^2
budget = {}
# 1. fit form and range
alts = [fit(u, y0, e0, u <= 1.0, 2)[0], fit(u, y0, e0, u <= 0.25, 2)[0], fit(u, y0, e0, u <= 0.125, 1)[0]]
budget["continuum fit form / range"] = max(abs(a - central) for a in alts)
# 2. propagated uncertainty of the individual points (finite size)
budget["finite size of individual points (propagated)"] = stat
# 3. finite size, coherent: all points shifted by their uncertainty in the direction of the last change
sgn = np.array([np.sign(c["m32"][-1] - c["m32"][-2]) for c in cont])
budget["finite size, coherent shift of all points"] = abs(fit(u, y0 + sgn * e0, e0, u <= 0.5, 2)[0] - central)
# 4. criterion: crossing masses instead of R = 3/2
yx = np.array([c["cross"] for c in cont]); ex = np.maximum(np.abs(yx - y0), 1e-5)
budget["criterion (crossing instead of R = 3/2)"] = abs(fit(u, yx, ex, u <= 0.5, 2)[0] - central)
# 4b. midpoint of the two criteria at every spacing, half their difference as uncertainty
ym = 0.5 * (yx + y0); em = np.maximum(0.5 * np.abs(yx - y0), 1e-5)
budget["criterion (midpoint of the two criteria)"] = abs(fit(u, ym, em, u <= 0.5, 2)[0] - central)
# 5. interpolation in m: linear between the two bracketing masses instead of the polynomial
ylin = []
for c in cont:
    x = c["x"]; N = c["sizes"][-1]
    ms = np.array(sorted(k[2] for k in A.D if k[0] == x and k[1] == N))
    Rv = np.array([(A.D[(x, N, m)]["E"][2] - A.D[(x, N, m)]["E"][0]) / (A.D[(x, N, m)]["E"][1] - A.D[(x, N, m)]["E"][0]) for m in ms])
    j = np.nonzero(np.sign(Rv[:-1] - 1.5) != np.sign(Rv[1:] - 1.5))[0][-1]
    ylin.append(ms[j] + (1.5 - Rv[j]) * (ms[j + 1] - ms[j]) / (Rv[j + 1] - Rv[j]))
ylin = np.array(ylin)
budget["interpolation in m (linear vs polynomial)"] = abs(fit(u, ylin, e0, u <= 0.5, 2)[0] - central)
# 6. weighting
budget["weighting (unweighted fit)"] = abs(fit(u, y0, np.ones_like(e0), u <= 0.5, 2)[0] - central)
# 7. leave one lattice spacing out
loo = []
for k in np.nonzero(u <= 0.5)[0]:
    sel = (u <= 0.5); sel[k] = False
    loo.append(fit(u, y0, e0, sel, 2 if sel.sum() > 2 else 1)[0])
budget["leave-one-spacing-out (max shift)"] = max(abs(a - central) for a in loo)
# 8. DMRG accuracy: energy error 3e-10 g on gaps of ~0.02 g, |dR/dm| ~ 10/g
budget["DMRG truncation and convergence"] = 3e-10 / 0.02 / 10
indep = ["continuum fit form / range", "finite size, coherent shift of all points", "criterion (midpoint of the two criteria)",
         "interpolation in m (linear vs polynomial)", "weighting (unweighted fit)", "leave-one-spacing-out (max shift)",
         "DMRG truncation and convergence"]
total = np.sqrt(sum(budget[k] ** 2 for k in indep))
conservative = np.sqrt(total ** 2 - budget["criterion (midpoint of the two criteria)"] ** 2
                       + budget["criterion (crossing instead of R = 3/2)"] ** 2)
print(f"critical mass, central fit (x >= 2, quadratic in a^2): {central:.5f}")
for k, v in budget.items():
    print(f"  {k:55s} {v:.1e}" + ("" if k in indep else "   (not added: alternative treatment of a source already counted)"))
print(f"  {'total (quadrature of the independent sources)':55s} {total:.1e}")
print(f"  {'conservative total (crossing criterion in full)':55s} {conservative:.1e}")
print(f"  reference 0.333561(4): difference {central - 0.333561:+.1e} = {abs(central - 0.333561) / total:.2f} sigma")
json.dump({"central": central, "budget": budget, "total": total, "conservative": conservative, "leave_one_out": loo, "m32_linear": ylin.tolist()},
          open("data/error_budget.json", "w"), indent=1, default=float)
