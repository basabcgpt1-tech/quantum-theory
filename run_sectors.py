"""Charge sectors Q = +-1 at the theta = pi critical point (x = 0.5).
Prediction (opposite fixed boundaries, character chi_eps): levels h = 1/2, 3/2, 5/2, 7/2, 9/2, 9/2,
so (E_n - E_0)/(E_1 - E_0) = 1, 2, 3, 4, 4 within the sector that realises it."""
import json, time
import numpy as np
import qed2

t0 = time.time(); out = {}
for dq in (-1, 1):
    for N in (10, 12, 14, 16, 18, 20):
        e = qed2.spectrum(N, 0.5, 0.356, 0.5, k=6, tol=1e-10, dq=dq)
        e0 = qed2.spectrum(N, 0.5, 0.356, 0.5, k=2, tol=1e-10)
        r = (e[1:] - e[0]) / (e[1] - e[0])
        # cross-sector offset in units of the Q=0 gap (universal part would be 1/4)
        off = (e[0] - e0[0]) / (e0[1] - e0[0])
        out[f"{dq}|{N}"] = {"ratios": r.tolist(), "offset": float(off), "E": e.tolist(), "E_Q0": e0.tolist()}
        print(f"Q={dq:+d} N={N}: ratios {np.round(r, 4)}  (E0(Q)-E0(0))/gap0 = {off:.4f}  [{time.time() - t0:.0f}s]",
              flush=True)
        json.dump(out, open("data/sectors.json", "w"), indent=1)
