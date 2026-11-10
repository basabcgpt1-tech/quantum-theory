"""Six lowest levels at the theta = pi critical point (x = 0.5), compared with the
Virasoro identity tower of the Ising CFT: (E_n - E_0)/(E_1 - E_0) = 1, 3/2, 2, 2, 5/2."""
import json, time
import numpy as np
import qed2

t0 = time.time(); out = {}
for N in (12, 14, 16, 18, 20):
    e = qed2.spectrum(N, 0.5, 0.356, 0.5, k=6, tol=1e-10)
    r = (e[1:] - e[0]) / (e[1] - e[0])
    out[N] = r.tolist()
    print(N, np.round(r, 4), f"[{time.time() - t0:.0f}s]", flush=True)
    json.dump(out, open("data/tower.json", "w"), indent=1)
