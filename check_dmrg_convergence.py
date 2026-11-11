"""Convergence check of the slowest DMRG run (N = 128, m/g = 0.358): energy per sweep,
and sensitivity to the truncation threshold and bond-dimension cap."""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import json, numpy as np, dmrg
N, x, m = 128, 0.5, 0.358
ref = json.load(open(f"data/dmrg/{x}_{N}_{m}.json"))["E"][0]
Ws = dmrg.mpo(N, x, m, 0.5)
E, psi = dmrg.solve(Ws, sweeps=45, verbose=True)
print(f"default (cut 1e-13, chi<=120): E0/g = {E / (2 * np.sqrt(x)):.12f}  stored {ref:.12f}  diff {E / (2 * np.sqrt(x)) - ref:.2e}")
E2, psi2 = dmrg.solve(Ws, sweeps=45, cut=1e-11, chi=40)
print(f"cut 1e-11, chi<=40: E0/g = {E2 / (2 * np.sqrt(x)):.12f}  diff from default {(E2 - E) / (2 * np.sqrt(x)):.2e}  chi {psi2.chi}")
