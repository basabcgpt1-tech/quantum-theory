"""Finite-size analysis of the theta = pi critical point (data/crit.json).

Tests of the claim "Ising CFT with equal fixed boundary conditions":
  R = (E2-E0)/(E1-E0) -> 3/2          (levels h = 2, 3 of the identity tower)
  c = -48 * D2[E0] / D2[E1-E0] -> 1/2 (Casimir energy; D2 = second difference in N)
  v = 1/(2 pi d(1/gap)/dL)            (velocity of light on the lattice)
m_c(N) is where R(m) for sizes N-2 and N cross.
"""
import json
import numpy as np

res = json.load(open("data/crit.json"))
keys = [k.split("|") for k in res]
out = {}
for x in sorted({float(k[0]) for k in keys}):
    ms = np.array(sorted({float(k[2]) for k in keys if float(k[0]) == x}))
    Ns = sorted({int(k[1]) for k in keys if float(k[0]) == x
                 and all(f"{x}|{k[1]}|{m}" in res for m in ms)})
    P = {}   # cubic interpolants in m of E0, gap, E2-E0
    for N in Ns:
        E = np.array([res[f"{x}|{N}|{m}"] for m in ms])
        P[N] = [np.polyfit(ms - ms.mean(), y, 3) for y in (E[:, 0], E[:, 1] - E[:, 0], E[:, 2] - E[:, 0])]
    ev = lambda N, i, m: np.polyval(P[N][i], m - ms.mean())
    R = lambda N, m: ev(N, 2, m) / ev(N, 1, m)
    print(f"\nx = {x}  (g a = {1 / np.sqrt(x):.2f})\n   N    m_c(N)   R      c       v")
    rows = []
    for N in Ns[2:]:
        grid = np.linspace(ms[0], ms[-1], 2001)
        d = R(N, grid) - R(N - 2, grid)
        j = np.nonzero(np.sign(d[:-1]) != np.sign(d[1:]))[0]
        if len(j) == 0:
            continue
        j = j[-1]
        mc = grid[j] - d[j] * (grid[j + 1] - grid[j]) / (d[j + 1] - d[j])
        e0 = [ev(n, 0, mc) for n in (N - 4, N - 2, N)]
        gp = [ev(n, 1, mc) for n in (N - 4, N - 2, N)]
        c = -48 * (e0[2] - 2 * e0[1] + e0[0]) / (gp[2] - 2 * gp[1] + gp[0])
        v = 1 / (2 * np.pi * (1 / gp[2] - 1 / gp[1]) / (2 / np.sqrt(x)))
        rows.append([N, mc, R(N, mc), c, v])
        print(f"  {N:2d}   {mc:.4f}  {R(N, mc):.4f}  {c:.4f}  {v:.4f}")
    rows = np.array(rows)
    if len(rows) >= 4:
        lim = {}
        for name, col in (("m_c", 1), ("R", 2), ("c", 3), ("v", 4)):
            ext = [np.polyfit(rows[-4:, 0] ** -p, rows[-4:, col], 1)[1] for p in (1.0, 2.0)]
            lim[name] = ext
            print(f"   N->inf {name:3s}: {ext[0]:.4f} (linear in 1/N)   {ext[1]:.4f} (linear in 1/N^2)")
        out[str(x)] = {"rows": rows.tolist(), "limits": lim}
json.dump(out, open("data/crit_summary.json", "w"), indent=1)
