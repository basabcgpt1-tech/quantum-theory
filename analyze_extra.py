"""Further Ising tests at the theta = pi critical point.
  nu    : dR/dm at the crossing grows like N^(1/nu)                (Ising: nu = 1)
  beta  : mid-chain electric field ~ (N + n0)^(-beta/nu), n0 from the gap (Ising: 1/8)
  c_ent : S(l) = (c/6) ln[(2N/pi) sin(pi l/N)] + const             (Ising: 1/2)
"""
import json
import numpy as np

out = {}
# ---- nu from the slope of R(m) at the crossing point
res = json.load(open("data/crit.json"))
summ = json.load(open("data/crit_summary.json"))
print("1/nu from d R/d m at the crossing (pairs N-2, N):")
for xs, blk in summ.items():
    x = float(xs)
    ms = np.array(sorted({float(k.split("|")[2]) for k in res if float(k.split("|")[0]) == x}))
    sl = []
    for N, mc, *_ in blk["rows"]:
        N = int(N)
        E = np.array([res[f"{x}|{N}|{m}"] for m in ms])
        R = (E[:, 2] - E[:, 0]) / (E[:, 1] - E[:, 0])
        p = np.polyfit(ms - mc, R, 3)
        sl.append((N, abs(p[2])))
    inv = [(b[0], np.log(b[1] / a[1]) / np.log(b[0] / a[0])) for a, b in zip(sl[:-1], sl[1:])]
    out[f"inv_nu_x{xs}"] = inv
    print(f"  x={xs}: " + "  ".join(f"{N}:{v:.3f}" for N, v in inv))

# ---- order parameter and entanglement at x = 0.5
g = json.load(open("data/groundstate.json"))
Ns = sorted(map(int, g))
x = 0.5
ms = np.array(sorted({float(k.split("|")[2]) for k in res if k.startswith("0.5|")}))
inv_gap = []
for N in Ns:
    if all(f"{x}|{N}|{m}" in res for m in ms):
        E = np.array([res[f"{x}|{N}|{m}"] for m in ms])
        inv_gap.append((N, 1 / np.interp(0.3585, ms, E[:, 1] - E[:, 0])))
ig = np.array(inv_gap)
s, i0 = np.polyfit(ig[-4:, 0], ig[-4:, 1], 1)
n0 = i0 / s
print(f"\nextrapolation length from 1/gap: n0 = {n0:.2f} sites")
mid = np.array([0.5 * (g[str(N)]["field"][N // 2 - 1] + g[str(N)]["field"][N // 2]) for N in Ns])
print("mid-chain field (two-link average):", np.round(mid, 5))
for lab, shift in (("n0 = 0", 0.0), (f"n0 = {n0:.2f} (from gap)", n0)):
    b = -np.polyfit(np.log(np.array(Ns[-5:]) + shift), np.log(mid[-5:]), 1)[0]
    out[f"beta_{'raw' if shift == 0 else 'n0'}"] = b
    print(f"  beta/nu with {lab}: {b:.4f}   (Ising 0.125)")

print("\ncentral charge from entanglement entropy:")
for N in Ns[-3:]:
    S = np.array(g[str(N)]["entropy"]); l = np.arange(1, N)
    chord = np.log(2 * N / np.pi * np.sin(np.pi * l / N))
    row = []
    for par in (0, 1):
        sel = (l % 2 == par) & (l >= 4) & (l <= N - 4)
        row.append(6 * np.polyfit(chord[sel], S[sel], 1)[0])
    # parity-averaged entropy on bonds removes the leading oscillation
    Sa = 0.5 * (S[:-1] + S[1:]); ca = 0.5 * (chord[:-1] + chord[1:]); sel = slice(3, N - 5)
    cav = 6 * np.polyfit(ca[sel], Sa[sel], 1)[0]
    out[f"c_ent_{N}"] = [row[0], row[1], cav]
    print(f"  N={N}: even cuts {row[0]:.3f}, odd cuts {row[1]:.3f}, averaged {cav:.3f}")
json.dump(out, open("data/extra_summary.json", "w"), indent=1)
