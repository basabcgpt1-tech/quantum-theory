"""DMRG production runs at theta = pi (three lowest neutral states, entropies, field profile).
Usage: python3 run_dmrg.py <set>   with set in {univ, cont}. Resumable; results in data/dmrg/."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["OPENBLAS_NUM_THREADS"] = "1"
import json, time
from multiprocessing import Pool
import numpy as np

SETS = {
    # universality study at ga = 1.41
    "univ": [(0.5, N, m) for N in (24, 32, 48, 64, 96, 128) for m in (0.352, 0.354, 0.356, 0.358, 0.360)],
    # continuum limit of the critical mass
    "cont": ([(1.0, N, m) for N in (48, 64) for m in (0.338, 0.342, 0.346, 0.350)]
             + [(2.0, N, m) for N in (64, 96) for m in (0.332, 0.336, 0.340, 0.344)]
             + [(4.0, N, m) for N in (96, 128) for m in (0.330, 0.334, 0.338, 0.342)]
             + [(8.0, N, m) for N in (96, 128) for m in (0.330, 0.334, 0.338, 0.342)]
             + [(16.0, N, m) for N in (96, 128) for m in (0.330, 0.334, 0.338, 0.342)]),
    # larger volumes at the finer spacings
    "cont2": ([(4.0, 192, m) for m in (0.333, 0.335, 0.337)] + [(8.0, 192, m) for m in (0.332, 0.334, 0.336)]
              + [(16.0, 192, m) for m in (0.332, 0.334, 0.336)]),
    # refined masses at ga = 1.41 for entanglement and Casimir analyses
    "fine": [(0.5, N, m) for N in (64, 96, 128) for m in (0.357, 0.3575)],
}


def job(args):
    import dmrg
    x, N, m = args
    f = f"data/dmrg/{x}_{N}_{m}.json"
    if os.path.exists(f):
        return f
    t = time.time()
    E, st = dmrg.spectrum(N, x, m, 0.5, k=3, full=True)
    json.dump({"x": x, "N": N, "m": m, "E": E.tolist(), "chi": [s.chi for s in st], "sweeps": [s.sweeps for s in st],
               "entropy": dmrg.entropies(st[0]), "field": dmrg.field(st[0], 0.5), "seconds": time.time() - t},
              open(f, "w"))
    print(f"{f}  R={(E[2] - E[0]) / (E[1] - E[0]):.5f}  [{time.time() - t:.0f}s]", flush=True)
    return f


# theta = 0 checks at fine lattice spacing: (x, N, m, eps0, k)
SETS["theta0"] = ([(x, N, m, 0.0, 2) for x in (4.0, 9.0, 16.0, 25.0) for N in (40, 56, 72) for m in (0.0, 0.02, 0.04)]
                  + [(x, N, m, al, 1) for x in (4.0, 9.0, 16.0) for N in (48, 64)
                     for m in (0.0, 0.025, 0.05, 0.075, 0.1) for al in (0.25, 0.5)]
                  + [(x, N, m, 0.0, 1) for x in (4.0, 9.0, 16.0) for N in (48, 64) for m in (0.0, 0.025, 0.05, 0.075, 0.1)])


def job0(args):
    import dmrg
    x, N, m, eps0, k = args
    f = f"data/dmrg0/{x}_{N}_{m}_{eps0}_{k}.json"
    if os.path.exists(f):
        return f
    t = time.time()
    E, st = dmrg.spectrum(N, x, m, eps0, k=k, full=True)
    json.dump({"x": x, "N": N, "m": m, "eps0": eps0, "E": E.tolist(), "chi": [s.chi for s in st],
               "seconds": time.time() - t}, open(f, "w"))
    print(f"{f}  [{time.time() - t:.0f}s]", flush=True)
    return f


if __name__ == "__main__":
    if sys.argv[1] == "theta0":
        os.makedirs("data/dmrg0", exist_ok=True)
        with Pool(4) as p:
            list(p.imap_unordered(job0, sorted(SETS["theta0"], key=lambda j: j[1] * j[4]), chunksize=1))
        print("done"); sys.exit()
    os.makedirs("data/dmrg", exist_ok=True)
    jobs = sorted(SETS[sys.argv[1]], key=lambda j: j[1] * (1 + 0.2 * np.log2(j[0] + 1)))
    with Pool(4) as p:
        list(p.imap_unordered(job, jobs, chunksize=1))
    print("done")
