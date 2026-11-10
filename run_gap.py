"""Massless gap M/g -> 1/sqrt(pi): infinite-volume fit at each lattice spacing."""
import json, time
import numpy as np
import qed2

t0 = time.time()
out = {}

# ---- A: gap of a free massive boson in a box: gap^2 = M^2 + c/(N + n0)^2
xs = [1.0, 2.0, 4.0, 6.0, 9.0, 12.0, 16.0]
Ns = np.arange(12, 23, 2)
A = {}
for x in xs:
    g2 = np.array([qed2.gap(int(N), x, 0.0) for N in Ns]) ** 2
    best = None
    for n0 in np.linspace(-2, 8, 401):
        B = np.vstack([np.ones(len(Ns)), 1.0 / (Ns + n0) ** 2]).T
        coef, r = np.linalg.lstsq(B, g2, rcond=None)[:2]
        if best is None or r[0] < best[0]:
            best = (r[0], coef[0], n0)
    A[x] = np.sqrt(best[1])
    print(f"A  x={x:5.1f}  L_max g={Ns[-1] / np.sqrt(x):5.2f}  M/g={A[x]:.5f}  n0={best[2]:.2f}"
          f"  [{time.time() - t0:.0f}s]", flush=True)
out["A"] = {str(k): v for k, v in A.items()}
json.dump(out, open("data/gap.json", "w"), indent=1)
