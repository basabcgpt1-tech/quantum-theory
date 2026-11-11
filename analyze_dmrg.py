"""Analysis of the DMRG runs (data/dmrg/*.json): universal quantities at ga = 1.41 and the
continuum limit of the critical mass."""
import glob, json
import numpy as np

D = {}
for f in glob.glob("data/dmrg/*.json"):
    d = json.load(open(f)); D[(d["x"], d["N"], d["m"])] = d
xs = sorted({k[0] for k in D})


def table(x):
    Ns = sorted({k[1] for k in D if k[0] == x}); out = {}
    for N in Ns:
        ms = np.array(sorted(k[2] for k in D if k[0] == x and k[1] == N))
        if len(ms) < 3:
            continue
        E = np.array([D[(x, N, m)]["E"] for m in ms])
        mid = np.array([0.5 * (D[(x, N, m)]["field"][N // 2 - 1] + D[(x, N, m)]["field"][N // 2]) for m in ms])
        deg = min(3, len(ms) - 1); m0 = ms.mean()
        out[N] = {"ms": ms, "m0": m0,
                  "E0": np.polyfit(ms - m0, E[:, 0], deg), "gap": np.polyfit(ms - m0, E[:, 1] - E[:, 0], deg),
                  "g2": np.polyfit(ms - m0, E[:, 2] - E[:, 0], deg), "mid": np.polyfit(ms - m0, mid, deg)}
    return out


def ev(T, N, q, m):
    return np.polyval(T[N][q], m - T[N]["m0"])


def R(T, N, m):
    return ev(T, N, "g2", m) / ev(T, N, "gap", m)


def cross(T, N1, N2):
    lo = max(T[N1]["ms"][0], T[N2]["ms"][0]) - 0.004; hi = min(T[N1]["ms"][-1], T[N2]["ms"][-1]) + 0.004
    g = np.linspace(lo, hi, 4001); d = R(T, N2, g) - R(T, N1, g)
    j = np.nonzero(np.sign(d[:-1]) != np.sign(d[1:]))[0]
    if len(j) == 0:
        return None
    j = j[-1]
    return g[j] - d[j] * (g[j + 1] - g[j]) / (d[j + 1] - d[j])


def three_halves(T, N):
    g = np.linspace(T[N]["ms"][0] - 0.004, T[N]["ms"][-1] + 0.004, 4001); d = R(T, N, g) - 1.5
    j = np.nonzero(np.sign(d[:-1]) != np.sign(d[1:]))[0]
    return None if len(j) == 0 else g[j[-1]] - d[j[-1]] * (g[1] - g[0]) / (d[j[-1] + 1] - d[j[-1]])


summary = {}
# ------------------------------------------------------------------ universal quantities, x = 0.5
if 0.5 in xs:
    x = 0.5; T = table(x); Ns = sorted(T)
    print("x = 0.5 (ga = 1.41): crossing of R between consecutive sizes")
    print("   N1-N2    m_c       R*       v        1/nu")
    rows = []
    for N1, N2 in zip(Ns[:-1], Ns[1:]):
        mc = cross(T, N1, N2)
        if mc is None:
            continue
        v = 1 / (2 * np.pi * (1 / ev(T, N2, "gap", mc) - 1 / ev(T, N1, "gap", mc)) / ((N2 - N1) / np.sqrt(x)))
        dR = [np.polyval(np.polyder(np.polyfit(T[N]["ms"] - mc, R(T, N, T[N]["ms"]), min(3, len(T[N]["ms"]) - 1))), 0.0)
              for N in (N1, N2)]
        inu = np.log(dR[1] / dR[0]) / np.log(N2 / N1)
        rows.append([N1, N2, mc, R(T, N2, mc), v, inu])
        print(f"  {N1:3d}-{N2:3d}  {mc:.5f}  {R(T, N2, mc):.4f}  {v:.4f}  {inu:.3f}")
    summary["cross"] = rows
    print("\n   N     m(R=3/2)")
    th = []
    for N in Ns:
        m32 = three_halves(T, N)
        if m32 is not None:
            th.append([N, m32]); print(f"  {N:3d}   {m32:.5f}")
    summary["m32"] = th
    mc = rows[-1][2]
    m32 = dict(th)
    mcr = {r[1]: r[2] for r in rows}

    def interp_m(N, vals_of_m, m):
        """Quadratic interpolation in m of a per-run quantity, using the three nearest masses."""
        ms = np.array(sorted(T[N]["ms"], key=lambda z: abs(z - m))[:3])
        return float(np.polyval(np.polyfit(ms - m, [vals_of_m(z) for z in ms], 2), 0.0))

    # central charge: E0 = e N + b - (c/48) gap, solved on consecutive triplets of sizes, evaluated
    # at the two pseudo-critical masses of the largest size of the triplet
    print("\nCasimir central charge (triplets of sizes) at m = crossing mass and at m = m(R=3/2):")
    cs = []
    for a_, b_, c_ in zip(Ns[:-2], Ns[1:-1], Ns[2:]):
        row = [c_]
        for m in (mcr[c_], m32[c_]):
            Mx = np.array([[n, 1.0, -ev(T, n, "gap", m) / 48] for n in (a_, b_, c_)])
            row.append(float(np.linalg.solve(Mx, [ev(T, n, "E0", m) for n in (a_, b_, c_)])[2]))
        cs.append(row); print(f"  N = {a_},{b_},{c_}:  {row[1]:.4f}   {row[2]:.4f}")
    summary["c_casimir"] = cs
    # order parameter exponent at the two estimates of the critical mass
    print("\nbeta/nu from the mid-chain field, local slope between consecutive sizes:")
    beta = []
    for lab, m in (("crossing", mc), ("R=3/2", th[-1][1])):
        ig = np.array([[N, 1 / ev(T, N, "gap", m)] for N in Ns])
        s_, i0 = np.polyfit(ig[-3:, 0], ig[-3:, 1], 1); n0 = i0 / s_
        mid = np.array([ev(T, N, "mid", m) for N in Ns])
        loc = [[Ns[i + 1], -np.log(mid[i + 1] / mid[i]) / np.log((Ns[i + 1] + n0) / (Ns[i] + n0)),
                -np.log(mid[i + 1] / mid[i]) / np.log(Ns[i + 1] / Ns[i])] for i in range(len(Ns) - 1)]
        beta.append({"m": m, "n0": n0, "rows": loc})
        print(f"  m = {m:.5f} ({lab}), n0 = {n0:.2f}: " + "  ".join(f"{N}:{b1:.4f}/{b0:.4f}" for N, b1, b0 in loc)
              + "   (with n0 / raw)")
    summary["beta"] = beta
    # entanglement: half-chain entropy against ln N at the critical mass; cuts of fixed parity
    print("\ncentral charge from the half-chain entanglement entropy, S = (c/6) ln(N + n0) + const:")
    ent = []
    for lab, off in (("cut at N/2 (even)", 0), ("cut at N/2 - 1 (odd)", -1)):
        for mlab, m in (("crossing", mc), ("R=3/2", th[-1][1])):
            Sv = np.array([interp_m(N, lambda z, N=N: D[(x, N, z)]["entropy"][N // 2 - 1 + off], m) for N in Ns])
            loc = [[Ns[i + 1], 6 * (Sv[i + 1] - Sv[i]) / np.log((Ns[i + 1] + n0) / (Ns[i] + n0))] for i in range(len(Ns) - 1)]
            ent.append({"cut": lab, "m": m, "S": Sv.tolist(), "rows": loc})
            print(f"  {lab}, m = {m:.5f} ({mlab}): " + "  ".join(f"{N}:{c:.3f}" for N, c in loc))
    summary["c_ent"] = ent
    # parity-averaged half-chain entropy, and the alternating part against the order parameter
    Sbar = lambda N, z: 0.5 * (D[(x, N, z)]["entropy"][N // 2 - 1]
                               + 0.5 * (D[(x, N, z)]["entropy"][N // 2 - 2] + D[(x, N, z)]["entropy"][N // 2]))
    av = []
    for mlab, m in (("crossing", mc), ("R=3/2", th[-1][1])):
        Sv = np.array([interp_m(N, lambda z, N=N: Sbar(N, z), m) for N in Ns])
        loc = [[Ns[i + 1], 6 * (Sv[i + 1] - Sv[i]) / np.log(Ns[i + 1] / Ns[i]),
                6 * (Sv[i + 1] - Sv[i]) / np.log((Ns[i + 1] + n0) / (Ns[i] + n0))] for i in range(len(Ns) - 1)]
        av.append({"m": m, "S": Sv.tolist(), "rows": loc})
        print(f"  parity-averaged, m = {m:.5f} ({mlab}): " + "  ".join(f"{N}:{c0:.3f}/{c1:.3f}" for N, c0, c1 in loc)
              + "   (raw / with n0)")
    summary["c_ent_avg"] = av
    mref = 0.358
    dS = np.array([D[(x, N, mref)]["entropy"][N // 2 - 1]
                   - 0.5 * (D[(x, N, mref)]["entropy"][N // 2 - 2] + D[(x, N, mref)]["entropy"][N // 2]) for N in Ns])
    op = np.array([0.5 * (D[(x, N, mref)]["field"][N // 2 - 1] + D[(x, N, mref)]["field"][N // 2]) for N in Ns])
    p = np.polyfit(op, dS, 1); res = np.abs(np.polyval(p, op) - dS).max()
    summary["alt_vs_op"] = {"m": mref, "dS": dS.tolist(), "field": op.tolist(), "slope": p[0], "intercept": p[1], "max_resid": res}
    print(f"  alternating part of S against mid-chain field at m = {mref}: dS = {p[1]:+.4f} {p[0]:+.4f} <L>,"
          f" max residual {res:.1e} (range of dS {np.ptp(dS):.3f})")

# ------------------------------------------------------------------ continuum limit of m_c
# At each lattice spacing: m_{3/2}(N) on every size, then N -> infinity from the two largest sizes
# assuming a correction N^-p, p = 2 and 3 (their spread is the finite-size uncertainty).
print("\ncritical mass versus lattice spacing (criterion R = 3/2)")
print("    x     ga    sizes (Lg of the two largest)         m_3/2(N)                      N->inf p=2   p=3     crossing(R*)")
cont = []
for x in xs:
    T = table(x); Ns = sorted(T)
    if len(Ns) < 2:
        continue
    vals = [(N, three_halves(T, N)) for N in Ns]
    vals = [(N, v) for N, v in vals if v is not None]
    (N1, a), (N2, b) = vals[-2], vals[-1]
    ex = [(N2 ** p * b - N1 ** p * a) / (N2 ** p - N1 ** p) for p in (2, 3)]
    mc = cross(T, Ns[-2], Ns[-1]); Rx = None if mc is None else R(T, Ns[-1], mc)
    cont.append({"x": x, "sizes": [N for N, _ in vals], "m32": [v for _, v in vals], "p2": ex[0], "p3": ex[1],
                 "cross": mc, "R_cross": Rx, "Lg": [N1 / np.sqrt(x), N2 / np.sqrt(x)]})
    print(f"  {x:5.2f}  {1 / np.sqrt(x):.3f}  {N1}-{N2} ({N1 / np.sqrt(x):.0f}, {N2 / np.sqrt(x):.0f})   "
          + " ".join(f"{v:.5f}" for _, v in vals[-3:]) + f"   {ex[0]:.5f}  {ex[1]:.5f}   "
          + ("" if mc is None else f"{mc:.5f} ({Rx:.4f})"))
summary["cont"] = cont
u = np.array([1 / c["x"] for c in cont]); mid = np.array([0.5 * (c["p2"] + c["p3"]) for c in cont])
err = np.array([max(0.5 * abs(c["p2"] - c["p3"]), 2e-5) for c in cont])
fits = {}
for lab, sel, deg in (("x >= 1, quadratic in a^2", u <= 1.0, 2), ("x >= 2, quadratic in a^2", u <= 0.5, 2),
                      ("x >= 2, linear in a^2", u <= 0.5, 1), ("x >= 4, linear in a^2", u <= 0.25, 1)):
    if sel.sum() > deg:
        vals_ = [np.polyfit(u[sel], y[sel], deg, w=1 / err[sel])[-1] for y in
                 (mid, np.array([c["p2"] for c in cont]), np.array([c["p3"] for c in cont]))]
        fits[lab] = vals_
        print(f"  continuum, {lab}: m_c/g = {vals_[0]:.5f}   (p=2: {vals_[1]:.5f}, p=3: {vals_[2]:.5f})")
summary["cont_fits"] = fits
# velocity of light at the critical mass of each lattice spacing (two largest sizes)
print("\nvelocity at criticality versus lattice spacing")
vel = []
for c in cont:
    x = c["x"]; T = table(x); Ns = sorted(T); m = 0.5 * (c["p2"] + c["p3"]); N1, N2 = Ns[-2], Ns[-1]
    v = 1 / (2 * np.pi * (1 / ev(T, N2, "gap", m) - 1 / ev(T, N1, "gap", m)) / ((N2 - N1) / np.sqrt(x)))
    vel.append([x, v]); print(f"  x={x:5.2f}  ga={1 / np.sqrt(x):.3f}  v={v:.4f}")
vv = np.array(vel); sel = (vv[:, 0] >= 2) & (vv[:, 0] <= 8)
v_a = np.polyfit(vv[sel, 0] ** -0.5, vv[sel, 1], 1)[1]; v_a2 = np.polyfit(1 / vv[sel, 0], vv[sel, 1], 1)[1]
summary["velocity"] = {"rows": vel, "linear_in_a": v_a, "linear_in_a2": v_a2}
print(f"  a -> 0 from x = 2, 4, 8: v = {v_a:.4f} (linear in a), {v_a2:.4f} (linear in a^2)")
print("  reference: 0.333561(4)  [Arguello Cruz, Tarnopolsky, Xin]")
json.dump(summary, open("data/dmrg_summary.json", "w"), indent=1, default=float)
