"""theta = 0 checks at fine lattice spacing from DMRG (data/dmrg0): massless gap, dM/dm, string tension."""
import glob, json
import numpy as np

D = {}
for f in glob.glob("data/dmrg0/*.json"):
    d = json.load(open(f)); D[(d["x"], d["N"], d["m"], d["eps0"], len(d["E"]))] = d["E"]
out = {}


def bulk_mass(x, m, Ns=(40, 56, 72)):
    """Solve gap^2 = M^2 + C/(N+n0)^2 exactly on three sizes."""
    g2 = np.array([(D[(x, N, m, 0.0, 2)][1] - D[(x, N, m, 0.0, 2)][0]) ** 2 for N in Ns]); N = np.array(Ns, float)
    f = lambda n0: ((g2[0] - g2[1]) / ((N[0] + n0) ** -2 - (N[1] + n0) ** -2)
                    - (g2[1] - g2[2]) / ((N[1] + n0) ** -2 - (N[2] + n0) ** -2))
    lo, hi = -5.0, 30.0
    grid = np.linspace(lo, hi, 3501); v = f(grid); j = np.nonzero(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    if len(j) == 0:      # no finite-size dependence resolvable: gap already converged
        return float(np.sqrt(g2[-1])), None
    n0 = grid[j[0]] - v[j[0]] * (grid[1] - grid[0]) / (v[j[0] + 1] - v[j[0]])
    C = (g2[1] - g2[2]) / ((N[1] + n0) ** -2 - (N[2] + n0) ** -2)
    return float(np.sqrt(g2[2] - C / (N[2] + n0) ** 2)), float(n0)


xs = sorted({k[0] for k in D if k[4] == 2})
print("massless gap and slope (DMRG, N = 40, 56, 72)")
rows = []
for x in xs:
    ms_x = sorted({k[2] for k in D if k[0] == x and k[4] == 2 and all((x, N, k[2], 0.0, 2) in D for N in (40, 56, 72))})
    M = np.array([bulk_mass(x, m)[0] for m in ms_x]); mm = np.array(ms_x)
    # slope at m = 0 from a polynomial in m: quadratic if three masses, cubic if more
    deg = 2 if len(mm) <= 3 else 3
    p = np.polyfit(mm, M, deg)
    two_pt = 2 * (M[mm == 0.02][0] - M[0]) / 0.02 - (M[mm == 0.04][0] - M[0]) / 0.04
    rows.append([x, M[0], p[-2], p[-3], two_pt, len(mm)])
    print(f"  x={x:5.1f}  Lg_max={72 / np.sqrt(x):5.1f}  M0/g={M[0]:.5f}  dM/dm={p[-2]:.4f} (degree-{deg} fit, {len(mm)} masses)"
          f"  curvature={p[-3]:+.3f}  [two-point estimate {two_pt:.4f}]")
out["gap"] = rows
r = np.array(rows)
p1 = np.polyfit(1 / r[-3:, 0], r[-3:, 1], 1)[1]; p2 = np.polyfit(1 / r[:, 0], r[:, 1], 2)[2]
out["M0/g"] = [p1, p2, 1 / np.sqrt(np.pi)]
print(f"  continuum M0/g: {p1:.5f} (linear in a^2, three finest)  {p2:.5f} (quadratic, all)   exact {1 / np.sqrt(np.pi):.5f}")
sel = r[:, 5] >= 5
if sel.sum() >= 2:
    q = np.polyfit(1 / r[sel, 0], r[sel, 2], 1); q2 = np.polyfit(r[sel, 0] ** -0.5, r[sel, 2], 1)
    out["dM/dm"] = [q[1], q2[1], np.exp(np.euler_gamma)]
    print(f"  continuum dM/dm: {q[1]:.4f} (linear in a^2)  {q2[1]:.4f} (linear in a)   exact {np.exp(np.euler_gamma):.4f}")
    qc = np.polyfit(1 / r[sel, 0], r[sel, 3], 1)
    out["curvature"] = qc[1]
    print(f"  continuum d2M/dm2 / 2: {qc[1]:+.3f}")

# slope at matched physical volume: the finite-size fit is biased at small Lg, so compare
# lattice spacings at equal Lg of the largest chain and extrapolate in a^2
print("\ndM/dm from triples of sizes (cubic fit in m, six masses)")
ms6 = [0.0, 0.005, 0.01, 0.02, 0.03, 0.04]; vol = {}
for x in (9.0, 16.0):
    for Ns in ((40, 56, 72), (56, 72, 96), (72, 96, 128)):
        if all((x, N, m, 0.0, 2) in D for N in Ns for m in ms6):
            M = np.array([bulk_mass(x, m, Ns)[0] for m in ms6]); p = np.polyfit(ms6, M, 3)
            vol[(x, Ns[-1])] = [Ns[-1] / np.sqrt(x), M[0], p[-2], p[-3]]
            print(f"  x={x:4.1f} N={Ns}  Lg_max={Ns[-1] / np.sqrt(x):4.1f}  M0/g={M[0]:.5f}  dM/dm={p[-2]:.4f}  curvature={p[-3]:+.3f}")
out["slope_volume"] = {f"{k[0]}|{k[1]}": v for k, v in vol.items()}
matched = []
for Lg, n9, n16 in ((24, 72, 96), (32, 96, 128)):
    if (9.0, n9) in vol and (16.0, n16) in vol:
        s9, s16 = vol[(9.0, n9)][2], vol[(16.0, n16)][2]; c9, c16 = vol[(9.0, n9)][3], vol[(16.0, n16)][3]
        s0 = s16 + (s16 - s9) * (1 / 16) / (1 / 9 - 1 / 16); c0 = c16 + (c16 - c9) * (1 / 16) / (1 / 9 - 1 / 16)
        matched.append([Lg, s9, s16, s0, c0])
        print(f"  matched volume Lg = {Lg}: dM/dm = {s9:.4f} (x=9), {s16:.4f} (x=16)  ->  a^2 -> 0: {s0:.4f}"
              f"   curvature -> {c0:+.3f}   (exact slope {np.exp(np.euler_gamma):.4f})")
out["slope_matched"] = matched

S = np.exp(np.euler_gamma) / (2 * np.pi ** 1.5)
b2 = lambda al: -(S * np.sqrt(np.pi)) ** 2 / 4 * (-8.9139) * (np.cos(4 * np.pi * al) - 1)
print("\nstring tension (DMRG, two largest sizes at each spacing)")
st = {}
for x in sorted({k[0] for k in D if k[4] == 1}):
    Ns = sorted({k[1] for k in D if k[0] == x and k[4] == 1})
    best = None
    for N1, N2 in zip(Ns[:-1], Ns[1:]):
        ms = sorted({k[2] for k in D if k[0] == x and k[1] == N2 and k[4] == 1 and k[3] == 0.0} - {0.0})
        if all((x, N, m, al, 1) in D for N in (N1, N2) for m in [0.0] + ms for al in (0.0, 0.25, 0.5)):
            best = (N1, N2, ms)
    if best is None:
        continue
    N1, N2, ms = best
    dens = lambda m, al: (D[(x, N2, m, al, 1)][0] - D[(x, N1, m, al, 1)][0]) / ((N2 - N1) / np.sqrt(x))
    for al in (0.25, 0.5):
        res0 = dens(0.0, al) - dens(0.0, 0.0)
        sig = np.array([dens(m, al) - dens(m, 0.0) for m in ms]); mm = np.array(ms)
        deg = 2 if len(mm) >= 5 else 1
        p = np.polyfit(mm, sig / mm, deg)
        st[f"{x}|{al}"] = {"N": [N1, N2], "Lg": [N1 / np.sqrt(x), N2 / np.sqrt(x)], "residual_m0": res0, "m": ms,
                           "sigma": sig.tolist(), "a": p[-1], "b": p[-2], "c3": p[-3] if deg == 2 else None}
        th = np.array([S * m * (1 - np.cos(2 * np.pi * al)) + b2(al) * m * m for m in ms])
        print(f"  x={x:4.1f} N={N1},{N2} (Lg={N1 / np.sqrt(x):.0f},{N2 / np.sqrt(x):.0f}) alpha={al}: residual(m=0) {res0:+.1e}"
              f"  a={p[-1]:.4f} b={p[-2]:+.3f}" + (f" c3={p[-3]:+.2f}" if deg == 2 else "")
              + "   lattice/theory " + " ".join(f"{v / t:.3f}" for v, t in zip(sig, th)))
out["string"] = st
for al in (0.25, 0.5):
    xs_ = sorted(float(k.split("|")[0]) for k in st if k.endswith(f"|{al}") and st[k]["c3"] is not None)
    if len(xs_) >= 3:
        u = 1 / np.array(xs_)
        res = {}
        for name in ("a", "b", "c3"):
            y = np.array([st[f"{x}|{al}"][name] for x in xs_])
            res[name] = [np.polyfit(u, y, 1)[1], np.polyfit(u[1:], y[1:], 1)[1] if len(u) >= 4 else np.nan, y.tolist()]
        out[f"string_cont_{al}"] = res
        print(f"  continuum alpha={al} (linear in a^2; all spacings / without the coarsest): "
              f"a={res['a'][0]:.4f} / {res['a'][1]:.4f} (theory {S * (1 - np.cos(2 * np.pi * al)):.4f})   "
              f"b={res['b'][0]:+.3f} / {res['b'][1]:+.3f} (theory {b2(al):+.3f})   c3={res['c3'][0]:+.2f} / {res['c3'][1]:+.2f}")
json.dump(out, open("data/theta0_summary.json", "w"), indent=1, default=float)
