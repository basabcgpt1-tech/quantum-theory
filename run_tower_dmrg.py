"""Conformal towers at the critical point of the x = 0.5 lattice (m/g = 0.3573) from DMRG:
six lowest levels in the charge sectors Q = 0 and Q = -1."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
import json, time
from multiprocessing import Pool
import numpy as np


def job(a):
    import dmrg
    Q, N = a; f = f"data/dmrg_tower/{Q}_{N}.json"
    if os.path.exists(f):
        return
    t = time.time(); E = dmrg.spectrum(N, 0.5, 0.3573, 0.5, k=6, Q0=Q)
    r = (E[1:] - E[0]) / (E[1] - E[0])
    json.dump({"Q": Q, "N": N, "E": E.tolist(), "ratios": r.tolist()}, open(f, "w"))
    print(f"Q={Q:+d} N={N}: {np.round(r, 4)}  [{time.time() - t:.0f}s]", flush=True)


if __name__ == "__main__":
    os.makedirs("data/dmrg_tower", exist_ok=True)
    with Pool(4) as p:
        list(p.imap_unordered(job, [(Q, N) for N in (24, 32, 48, 64) for Q in (0, -1)], chunksize=1))
    print("done")
