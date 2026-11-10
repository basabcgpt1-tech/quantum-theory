"""dM/dm at m = 0 (prediction e^gamma = 1.781), using the same infinite-volume fit as run_gap.py."""
import json
import numpy as np
import qed2

Ns = np.arange(12, 21, 2)


def bulk(x, m):
    g2 = np.array([qed2.gap(int(N), x, m) for N in Ns]) ** 2
    best = None
    for n0 in np.linspace(-2, 8, 401):
        B = np.vstack([np.ones(len(Ns)), 1.0 / (Ns + n0) ** 2]).T
        coef, r = np.linalg.lstsq(B, g2, rcond=None)[:2]
        if best is None or r[0] < best[0]:
            best = (r[0], coef[0])
    return np.sqrt(best[1])


out = {}
for x in (2.0, 4.0, 9.0):
    M = {m: bulk(x, m) for m in (0.0, 0.02, 0.04)}
    s1, s2 = (M[0.02] - M[0.0]) / 0.02, (M[0.04] - M[0.0]) / 0.04
    out[x] = {"M": M, "slope_m->0": 2 * s1 - s2}
    print(f"x={x}: M={[round(float(v), 5) for v in M.values()]}  slopes {s1:.4f}, {s2:.4f}  -> m->0: {2 * s1 - s2:.4f}",
          flush=True)
json.dump({str(k): {"M": {str(a): float(b) for a, b in v["M"].items()}, "slope": float(v["slope_m->0"])}
           for k, v in out.items()}, open("data/slope.json", "w"), indent=1)
print("prediction", np.exp(np.euler_gamma))
