"""Ground-state observables at the theta = pi critical point (x = 0.5, m/g = 0.3585):
electric-field profile <L_n> (the Ising order parameter) and entanglement entropy S(l)."""
import json, time
import numpy as np
import qed2

X, M = 0.5, 0.3585
t0 = time.time(); out = {}
for N in (8, 10, 12, 14, 16, 18, 20, 22):
    E, psi = qed2.spectrum(N, X, M, 0.5, k=1, tol=1e-12, vec=True)
    psi /= np.linalg.norm(psi)
    n, tab, A, B, Mt, cs = qed2.tables(N)
    p = psi ** 2
    field = (0.5 + cs.astype(np.float64) @ p).tolist()          # <L_n>, links n = 1..N-1
    st = qed2.basis_states(N)
    full = np.zeros(1 << N); full[st] = psi
    S = []
    for l in range(1, N):                                         # cut after site l
        sv = np.linalg.svd(full.reshape(1 << (N - l), 1 << l), compute_uv=False)
        w = sv[sv > 1e-14] ** 2
        S.append(float(-(w * np.log(w)).sum()))
    out[N] = {"E0": float(E[0]), "field": field, "entropy": S}
    print(f"N={N}: <L> mid two links {field[N // 2 - 1]:.5f} {field[N // 2]:.5f}; S(N/2)={S[N // 2 - 1]:.5f}"
          f"  [{time.time() - t0:.0f}s]", flush=True)
    json.dump(out, open("data/groundstate.json", "w"))
