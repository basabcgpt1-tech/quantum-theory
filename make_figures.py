"""Figures for paper/paper.md (reads data/*.json; writes paper/fig/*.svg)."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]                   # categorical, fixed order
SEQ = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0b0b0b"]   # one hue, light -> dark
MK = ["o", "s", "^", "D"]
INK, MUTED = "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e1",
                     "grid.linewidth": 0.6, "lines.linewidth": 1.6, "lines.markersize": 5,
                     "legend.frameon": False, "figure.facecolor": "white", "axes.facecolor": "white",
                     "svg.fonttype": "none", "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans"})
J = lambda f: json.load(open(f"data/{f}.json"))
save = lambda fig, name: (fig.tight_layout(), fig.savefig(f"paper/fig/{name}.svg"), plt.close(fig))
SQPI = 1 / np.sqrt(np.pi)

# ---------------------------------------------------------------- Fig 1: massless gap
raw_f = "data/gap_raw.json"
if not os.path.exists(raw_f):
    import qed2
    raw = {str(x): [float(qed2.gap(N, x, 0.0)) for N in range(8, 21, 2)] for x in (1.0, 4.0, 16.0)}
    json.dump(raw, open(raw_f, "w"))
raw = json.load(open(raw_f)); gap = J("gap")
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
Ns = np.arange(8, 21, 2)
for i, (x, g) in enumerate(raw.items()):
    x = float(x); L = Ns / np.sqrt(x); g = np.array(g)
    ax[0].plot(1 / L ** 2, g ** 2, MK[i], color=C[i], label=f"$ga={1 / np.sqrt(x):.2f}$")
ax[0].axhline(1 / np.pi, color=MUTED, lw=1, ls="--"); ax[0].text(0.002, 1 / np.pi + 0.02, r"$1/\pi$", color=MUTED)
ax[0].set_xlabel(r"$1/(Lg)^2$"); ax[0].set_ylabel(r"$(\Delta/g)^2$"); ax[0].legend(); ax[0].set_xlim(0, None)
ax[0].set_title("(a) finite-volume gap, $m=0$", loc="left", fontsize=9)
xs = np.array(sorted(map(float, gap["A"]))); M = np.array([gap["A"][str(x)] for x in xs])
ax[1].plot(1 / xs, M, "o", color=C[0], label="lattice, $L\\to\\infty$")
ax[1].axhline(SQPI, color=MUTED, lw=1, ls="--"); ax[1].text(0.6, SQPI + 0.0002, r"exact $1/\sqrt{\pi}$", color=MUTED)
ax[1].set_xlabel(r"$(ga)^2 = 1/x$"); ax[1].set_ylabel(r"$M/g$"); ax[1].set_ylim(0.560, 0.566)
ax[1].set_title("(b) infinite-volume gap vs lattice spacing", loc="left", fontsize=9)
save(fig, "fig1_gap")

# ---------------------------------------------------------------- Fig 2: string tension
d = J("string"); S = np.exp(np.euler_gamma) / (2 * np.pi ** 1.5)
b2 = lambda al: -(S * np.sqrt(np.pi)) ** 2 / 4 * (-8.9139) * (np.cos(4 * np.pi * al) - 1)
ms = np.array([0.025, 0.05, 0.075, 0.1, 0.15]); mm = np.linspace(0, 0.16, 50)
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), sharex=True)
for j, al in enumerate((0.25, 0.5)):
    a1 = S * (1 - np.cos(2 * np.pi * al))
    ax[j].plot(mm, a1 + 0 * mm, color=MUTED, ls="--", lw=1.2, label="1st order")
    ax[j].plot(mm, a1 + b2(al) * mm, color=INK, lw=1.2, label="1st + 2nd order")
    for i, x in enumerate((1.0, 2.0, 3.0)):
        s = np.array([d[f"{x}|{m}|{al}"] for m in ms])
        ax[j].plot(ms, s / ms, MK[i], color=C[i], label=f"lattice $ga={1 / np.sqrt(x):.2f}$")
    ax[j].set_xlabel("$m/g$"); ax[j].set_title(f"({'ab'[j]}) $\\alpha = {'1/4' if j == 0 else '1/2'}$", loc="left", fontsize=9)
ax[0].set_ylabel(r"$\sigma/(mg)$"); ax[0].legend(fontsize=7.5); ax[0].set_ylim(0.09, 0.2); ax[1].set_ylim(0.26, 0.37)
save(fig, "fig2_string")

# ---------------------------------------------------------------- Fig 3: effective potentials
c = S; phi = np.linspace(-1.3, 1.3, 600)
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
for i, m in enumerate((0.08, 0.158, 0.30)):
    V = phi ** 2 / (2 * np.pi) + c * m * (np.cos(2 * np.sqrt(np.pi) * phi) - 1)
    ax[0].plot(phi, V, color=C[i], label=f"$m/g={m}$")
ax[0].set_xlabel(r"$\phi$"); ax[0].set_ylabel(r"$[V(\phi)-V(0)]/g^2$"); ax[0].legend(fontsize=8)
ax[0].set_title(r"(a) classical potential, $\theta=\pi$", loc="left", fontsize=9); ax[0].set_ylim(-0.03, 0.12)
def gauss(ph, m):
    k = c * np.pi * m; a = 4 * np.pi * k * np.cos(2 * np.sqrt(np.pi) * ph)
    mu = (-a + np.sqrt(a * a + 4)) / 2
    return (mu ** 2 - 1) / (8 * np.pi) - np.log(mu ** 2) / (8 * np.pi) + 0.5 * ph ** 2 + k * mu * np.cos(2 * np.sqrt(np.pi) * ph)
for i, m in enumerate((0.20, 0.265, 0.33)):
    e = gauss(phi, m); ax[1].plot(phi, e - gauss(0.0, m), color=C[i], label=f"$m/g={m}$")
ax[1].set_xlabel(r"$\phi_0$"); ax[1].set_ylabel(r"$[\mathcal{E}(\phi_0)-\mathcal{E}(0)]/M^2$"); ax[1].legend(fontsize=8)
ax[1].set_title("(b) Gaussian effective potential", loc="left", fontsize=9); ax[1].set_ylim(-0.08, 0.12)
save(fig, "fig3_potential")

# ---------------------------------------------------------------- Fig 4: R(m) crossing, m_c(N)
res = J("crit"); summ = J("crit_summary")
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
x = 0.5; ms = np.array(sorted({float(k.split("|")[2]) for k in res if k.startswith("0.5|")}))
NN = [N for N in (12, 14, 16, 18, 20, 22, 24) if all(f"{x}|{N}|{m}" in res for m in ms)][-6:]
for i, N in enumerate(NN):
    E = np.array([res[f"{x}|{N}|{m}"] for m in ms]); R = (E[:, 2] - E[:, 0]) / (E[:, 1] - E[:, 0])
    ax[0].plot(ms, R, "-o", color=SEQ[i], ms=3.5, label=f"$N={N}$")
ax[0].axhline(1.5, color=MUTED, lw=1, ls="--"); ax[0].text(0.336, 1.503, "3/2", color=MUTED)
ax[0].set_xlabel("$m/g$"); ax[0].set_ylabel(r"$R=(E_2-E_0)/(E_1-E_0)$"); ax[0].legend(fontsize=7.5, ncol=2)
ax[0].set_title("(a) level ratio, $ga=1.41$", loc="left", fontsize=9)
for i, (xs_, blk) in enumerate(summ.items()):
    r = np.array(blk["rows"]); ax[1].plot(1 / r[:, 0] ** 2, r[:, 2], "-" + MK[i], color=C[i],
                                         label=f"$ga={1 / np.sqrt(float(xs_)):.2f}$")
ax[1].axhline(1.5, color=MUTED, lw=1, ls="--"); ax[1].set_xlim(0, None)
ax[1].set_xlabel("$1/N^2$"); ax[1].set_ylabel("$R$ at the crossing"); ax[1].legend(fontsize=8)
ax[1].set_title("(b) ratio at the crossing point", loc="left", fontsize=9)
save(fig, "fig4_crossing")

# ---------------------------------------------------------------- Fig 5: conformal towers
tw = J("tower"); sec = J("sectors")
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.2), sharex=True)
Ns = sorted(map(int, tw)); r = np.array([tw[str(N)] for N in Ns])
for j in range(1, 5):
    ax[0].plot(1 / np.array(Ns) ** 2, r[:, j], "-o", color=C[0], ms=4)
for h in (1.5, 2, 2.5):
    ax[0].axhline(h, color=MUTED, lw=1, ls="--")
ax[0].set_ylim(1.3, 2.6); ax[0].set_xlim(0, None)
ax[0].set_xlabel("$1/N^2$"); ax[0].set_ylabel(r"$(E_n-E_0)/(E_1-E_0)$")
ax[0].set_title(r"(a) $Q=0$: identity tower $\frac{3}{2},2,2,\frac{5}{2}$", loc="left", fontsize=9)
Ns = sorted(int(k.split("|")[1]) for k in sec if k.startswith("-1|"))
r = np.array([sec[f"-1|{N}"]["ratios"] for N in Ns])
for j in range(1, 5):
    ax[1].plot(1 / np.array(Ns) ** 2, r[:, j], "-s", color=C[1], ms=4)
for h in (2, 3, 4):
    ax[1].axhline(h, color=MUTED, lw=1, ls="--")
ax[1].set_ylim(1.7, 4.2); ax[1].set_xlabel("$1/N^2$")
ax[1].set_title(r"(b) $Q=-1$: $\varepsilon$ tower $2,3,4,4$", loc="left", fontsize=9)
save(fig, "fig5_towers")

# ---------------------------------------------------------------- Fig 6: c and v (separate panels)
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), sharex=True)
for i, (xs_, blk) in enumerate(summ.items()):
    r = np.array(blk["rows"]); lab = f"$ga={1 / np.sqrt(float(xs_)):.2f}$"
    ax[0].plot(1 / r[:, 0], r[:, 3], "-" + MK[i], color=C[i], label=lab)
    ax[1].plot(1 / r[:, 0], r[:, 4], "-" + MK[i], color=C[i], label=lab)
ax[0].axhline(0.5, color=MUTED, lw=1, ls="--"); ax[1].axhline(1.0, color=MUTED, lw=1, ls="--")
ax[0].set_xlim(0, None); ax[0].set_xlabel("$1/N$"); ax[1].set_xlabel("$1/N$")
ax[0].set_ylabel("central charge $c$"); ax[1].set_ylabel("velocity $v$"); ax[0].legend(fontsize=8)
ax[0].set_title("(a) Casimir central charge", loc="left", fontsize=9)
ax[1].set_title("(b) speed of light", loc="left", fontsize=9)
save(fig, "fig6_c_v")

# ---------------------------------------------------------------- Fig 7: order parameter, entanglement
g = J("groundstate"); ex = J("extra_summary")
Ns = np.array(sorted(map(int, g)))
mid = np.array([0.5 * (g[str(N)]["field"][N // 2 - 1] + g[str(N)]["field"][N // 2]) for N in Ns])
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
n0 = 3.78
ax[0].loglog(Ns + n0, mid, "o", color=C[0], label="lattice")
ref = mid[-1] * ((Ns + n0) / (Ns[-1] + n0)) ** -0.125
ax[0].loglog(Ns + n0, ref, color=MUTED, ls="--", lw=1.2, label=r"slope $-1/8$")
ax[0].set_xlabel(r"$N+n_0$"); ax[0].set_ylabel(r"mid-chain $\langle L\rangle$"); ax[0].legend(fontsize=8)
ax[0].set_xticks([12, 16, 20, 26]); ax[0].set_xticks([], minor=True); ax[0].set_xticklabels(["12", "16", "20", "26"])
ax[0].set_yticks([0.31, 0.32, 0.33, 0.34]); ax[0].set_yticks([], minor=True)
ax[0].set_yticklabels(["0.31", "0.32", "0.33", "0.34"])
ax[0].set_title("(a) order parameter", loc="left", fontsize=9)
N = int(Ns[-1]); Sx = np.array(g[str(N)]["entropy"]); l = np.arange(1, N)
chord = np.log(2 * N / np.pi * np.sin(np.pi * l / N))
ax[1].plot(chord[l % 2 == 0], Sx[l % 2 == 0], "o", color=C[0], label="even cuts")
ax[1].plot(chord[l % 2 == 1], Sx[l % 2 == 1], "s", color=C[1], label="odd cuts")
xx = np.linspace(chord.min(), chord.max(), 10)
ax[1].plot(xx, Sx[N // 2 - 1 - (N // 2) % 2] + (xx - chord.max()) / 12, color=MUTED, ls="--", lw=1.2, label="slope $c/6$, $c=1/2$")
ax[1].set_xlabel(r"$\ln[(2N/\pi)\sin(\pi\ell/N)]$"); ax[1].set_ylabel(r"$S(\ell)$"); ax[1].legend(fontsize=8)
ax[1].set_title(f"(b) entanglement entropy, $N={N}$", loc="left", fontsize=9)
save(fig, "fig7_order_ent")
print("figures written")
