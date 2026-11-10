"""3 lowest levels at theta = pi near the critical mass, for several lattice spacings.
Resumable; writes data/crit.json with keys "x|N|m" -> [E0, E1, E2] in units of g."""
import json, os, sys, time
import numpy as np
import qed2

WINDOWS = {0.25: (0.37, 0.42), 0.5: (0.335, 0.385), 1.0: (0.31, 0.36)}
out = "data/crit.json"
res = json.load(open(out)) if os.path.exists(out) else {}
Ns = [int(a) for a in sys.argv[1].split(",")]
if len(sys.argv) > 2:      # optional: restrict to some lattice spacings
    WINDOWS = {x: WINDOWS[x] for x in map(float, sys.argv[2].split(","))}
t0 = time.time()
for N in Ns:
    for x, (lo, hi) in WINDOWS.items():
        for m in np.round(np.arange(lo, hi + 1e-9, 0.01), 4):
            key = f"{x}|{N}|{m}"
            if key not in res:
                res[key] = qed2.spectrum(N, x, m, 0.5, k=3, tol=1e-10).tolist()
                json.dump(res, open(out, "w"))
        print(f"N={N} x={x} done [{time.time() - t0:.0f}s]", flush=True)
