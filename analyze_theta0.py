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
    try:
        M = [bulk_mass(x, m)[0] for m in (0.0, 0.02, 0.04)]
    except KeyError:
        continue
    s1, s2 = (M[1] - M[0]) / 0.02, (M[2] - M[0]) / 0.04
    rows.append([x, M[0], 2 * s1 - s2, 72 / np.sqrt(x)])
    print(f"  x={x:5.1f}  Lg_max={72 / np.sqrt(x):5.1f}  M0/g={M[0]:.5f}  dM/dm={2 * s1 - s2:.4f}")
out["gap"] = rows
r = np.array(rows)
if len(r) >= 3:
    for lab, col, exact in (("M0/g", 1, 1 / np.sqrt(np.pi)), ("dM/dm", 2, np.exp(np.euler_gamma))):
        p1 = np.polyfit(1 / r[-3:, 0], r[-3:, col], 1)[1]
        p2 = np.polyfit(1 / r[:, 0], r[:, col], 2)[2] if len(r) >= 4 else np.nan
        out[lab] = [p1, p2, exact]
        print(f"  continuum {lab}: {p1:.5f} (linear in a^2, three finest)  {p2:.5f} (quadratic, all)   exact {exact:.5f}")

S = np.exp(np.euler_gamma) / (2 * np.pi ** 1.5)
b2 = lambda al: -(S * np.sqrt(np.pi)) ** 2 / 4 * (-8.9139) * (np.cos(4 * np.pi * al) - 1)
ms = [0.025, 0.05, 0.075, 0.1]
print("\nstring tension (DMRG, N = 48, 64)")
st = {}
for x in sorted({k[0] for k in D if k[4] == 1}):
    dens = lambda m, al: (D[(x, 64, m, al, 1)][0] - D[(x, 48, m, al, 1)][0]) / (16 / np.sqrt(x))
    try:
        for al in (0.25, 0.5):
            res0 = dens(0.0, al) - dens(0.0, 0.0)
            sig = np.array([dens(m, al) - dens(m, 0.0) for m in ms])
            b, a = np.polyfit(ms, sig / np.array(ms), 1)
            st[f"{x}|{al}"] = {"residual_m0": res0, "sigma": sig.tolist(), "a": a, "b": b}
            th = np.array([S * m * (1 - np.cos(2 * np.pi * al)) + b2(al) * m * m for m in ms])
            print(f"  x={x:4.1f} alpha={al}: residual(m=0) {res0:+.1e}  a={a:.4f} b={b:+.3f}   lattice/theory "
                  + " ".join(f"{v / t:.3f}" for v, t in zip(sig, th)))
    except KeyError:
        continue
out["string"] = st
for al in (0.25, 0.5):
    xs_ = sorted(float(k.split("|")[0]) for k in st if k.endswith(f"|{al}"))
    if len(xs_) >= 3:
        A = np.polyfit(1 / np.array(xs_), [st[f"{x}|{al}"]["a"] for x in xs_], 1)[1]
        B = np.polyfit(1 / np.array(xs_), [st[f"{x}|{al}"]["b"] for x in xs_], 1)[1]
        out[f"string_cont_{al}"] = [A, B]
        print(f"  continuum alpha={al}: a={A:.4f} (theory {S * (1 - np.cos(2 * np.pi * al)):.4f})  b={B:+.3f} (theory {b2(al):+.3f})")
json.dump(out, open("data/theta0_summary.json", "w"), indent=1)
