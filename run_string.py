"""String tension sigma(alpha, m) between external charges +-alpha g.

sigma = bulk energy density with background field alpha minus that with none,
taken from the difference of ground energies at two sizes (removes boundary terms).
Only volumes L g >= 12 are used: at m = 0 sigma must vanish, and that residual is
printed as the finite-volume diagnostic.
"""
import json, time
import numpy as np
import qed2

t0 = time.time()
ms = [0.0, 0.025, 0.05, 0.075, 0.1, 0.15]
out = {}
for x, pair in ((1.0, (16, 20)), (2.0, (18, 22)), (3.0, (18, 22))):
    for m in ms:
        dens = {}
        for al in (0.0, 0.25, 0.5):
            e = [qed2.spectrum(N, x, m, al, k=1, tol=1e-10)[0] for N in pair]
            dens[al] = (e[1] - e[0]) / ((pair[1] - pair[0]) / np.sqrt(x))
        for al in (0.25, 0.5):
            out[f"{x}|{m}|{al}"] = dens[al] - dens[0.0]
    print(f"x={x} N={pair} done [{time.time() - t0:.0f}s]", flush=True)
    json.dump(out, open("data/string.json", "w"), indent=1)
