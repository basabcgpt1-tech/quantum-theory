"""Two-site DMRG for the lattice Schwinger Hamiltonian W of qed2.py (same units and conventions).

The Coulomb term sum_n L_n^2 is written as an exact matrix-product operator of bond
dimension 6 (including a penalty lam*(Q - Q0)^2 that selects the charge sector).
Excited states are obtained by penalising overlap with the states already found.
Dense real tensors, no symmetry blocks; validated against qed2.py in test_dmrg.py.
"""
import numpy as np

I2 = np.eye(2); SZ = np.diag([1.0, -1.0]); SP = np.array([[0.0, 1.0], [0.0, 0.0]]); SM = SP.T


def mpo(N, x, m, eps0=0.0, shift=True, lam=4.0, Q0=0):
    """List of MPO tensors W[k][wl, wr, s_out, s_in] for sites k = 1..N."""
    ml = m - 1.0 / (8 * np.sqrt(x)) if shift else m
    mu = 2.0 * ml * np.sqrt(x)
    Ws = []
    for k in range(1, N + 1):
        sg = (-1.0) ** k
        q = 0.5 * (SZ + sg * I2); s = 0.5 * SZ
        on = (0.5 * mu * (I2 + sg * SZ) + (N - k) * (q @ q + 2 * eps0 * q)
              + ((N - 1) * eps0 ** 2 / N) * I2 + lam * (0.25 * I2 - 2 * Q0 * s + (Q0 ** 2 / N) * I2))
        W = np.zeros((6, 6, 2, 2))
        W[0, 0] = I2; W[0, 1] = q; W[0, 2] = s; W[0, 3] = SP; W[0, 4] = SM; W[0, 5] = on
        W[1, 1] = I2; W[1, 5] = 2 * (N - k) * q
        W[2, 2] = I2; W[2, 5] = 2 * lam * s
        W[3, 5] = x * SM; W[4, 5] = x * SP; W[5, 5] = I2
        Ws.append(W)
    Ws[0] = Ws[0][0:1]; Ws[-1] = Ws[-1][:, 5:6]
    return Ws


def _lanczos(matvec, v0, kmax=40, tol=1e-11):
    n = v0.size; Q = np.empty((kmax, n)); al = []; be = []
    q = v0 / np.linalg.norm(v0); eprev = None
    for j in range(kmax):
        Q[j] = q; w = matvec(q); a = q @ w; al.append(a)
        w -= Q[:j + 1].T @ (Q[:j + 1] @ w); w -= Q[:j + 1].T @ (Q[:j + 1] @ w)
        b = np.linalg.norm(w)
        T = np.diag(al) + np.diag(be, 1) + np.diag(be, -1)
        ev, U = np.linalg.eigh(T)
        if b < 1e-12 or (eprev is not None and abs(ev[0] - eprev) < tol) or j == kmax - 1:
            return ev[0], U[:, 0] @ Q[:j + 1]
        eprev = ev[0]; be.append(b); q = w / b


class State:
    def __init__(self, N, rng, pattern=None):
        self.A = []
        for k in range(N):
            a = np.zeros((1, 2, 1)); up = (k % 2 == 1) if pattern is None else pattern[k]
            a[0, 0 if up else 1, 0] = 1.0
            self.A.append(a + 1e-3 * rng.standard_normal((1, 2, 1)))
        self.S = [None] * (N - 1)


def solve(Ws, prev=(), chi=120, cut=1e-13, sweeps=30, tol=1e-10, weight=20.0, seed=0, init=None, noise=0.05,
          verbose=False):
    """Lowest state of the MPO orthogonal to the states in `prev`. Returns (energy, State)."""
    N = len(Ws); rng = np.random.default_rng(seed)
    psi = State(N, rng)
    if init is not None:                                   # start from a copy of another state
        psi.A = [a.copy() for a in init.A]
        psi.A[N // 2] = psi.A[N // 2] + noise * rng.standard_normal(psi.A[N // 2].shape)
    A = psi.A
    for k in range(N - 1, 0, -1):                          # right-canonicalise
        D0, d, D1 = A[k].shape
        U, S, Vt = np.linalg.svd(A[k].reshape(D0, d * D1), full_matrices=False)
        A[k] = Vt.reshape(-1, d, D1); A[k - 1] = np.tensordot(A[k - 1], U * S, axes=(2, 0))
    A[0] /= np.linalg.norm(A[0])
    L = [None] * N; R = [None] * N
    L[0] = np.ones((1, 1, 1)); R[N - 1] = np.ones((1, 1, 1))
    for k in range(N - 1, 1, -1):
        R[k - 1] = _env_r(R[k], A[k], Ws[k])
    OL = [[None] * N for _ in prev]; OR = [[None] * N for _ in prev]
    for p, ph in enumerate(prev):
        OL[p][0] = np.ones((1, 1)); OR[p][N - 1] = np.ones((1, 1))
        for k in range(N - 1, 1, -1):
            OR[p][k - 1] = np.tensordot(A[k], np.tensordot(ph.A[k], OR[p][k], axes=(2, 1)), axes=([1, 2], [1, 2]))
    E = None; hist = []
    for sw in range(sweeps):
        chi_sw = chi if init is not None else min(chi, 16 * 2 ** sw)     # cheap first sweeps from a product state
        for direction in (+1, -1):
            rng_i = range(0, N - 1) if direction > 0 else range(N - 2, -1, -1)
            for i in rng_i:
                th = np.tensordot(A[i], A[i + 1], axes=(2, 0))
                shp = th.shape; Li, Ri, W1, W2 = L[i], R[i + 1], Ws[i], Ws[i + 1]
                P = []
                for p, ph in enumerate(prev):
                    t = np.tensordot(ph.A[i], ph.A[i + 1], axes=(2, 0))
                    t = np.tensordot(OL[p][i], t, axes=(1, 0)); t = np.tensordot(t, OR[p][i + 1], axes=(3, 1))
                    P.append(t.ravel())

                def mv(v):
                    T = np.tensordot(Li, v.reshape(shp), axes=(2, 0))            # b w s1 s2 c
                    T = np.tensordot(T, W1, axes=([1, 2], [0, 3]))               # b s2 c u t1
                    T = np.tensordot(T, W2, axes=([3, 1], [0, 3]))               # b c t1 v t2
                    T = np.tensordot(T, Ri, axes=([1, 3], [2, 1]))               # b t1 t2 d
                    out = T.ravel()
                    for pv in P:
                        out = out + weight * pv * (pv @ v)
                    return out
                E, v = _lanczos(mv, th.ravel())
                M = v.reshape(shp[0] * 2, 2 * shp[3])
                U, S, Vt = np.linalg.svd(M, full_matrices=False)
                w2 = S ** 2; keep = max(1, min(chi_sw, int(np.sum(np.cumsum(w2[::-1])[::-1] > cut * w2.sum()))))
                U, S, Vt = U[:, :keep], S[:keep], Vt[:keep]; S = S / np.linalg.norm(S)
                psi.S[i] = S
                if direction > 0:
                    A[i] = U.reshape(shp[0], 2, keep); A[i + 1] = (S[:, None] * Vt).reshape(keep, 2, shp[3])
                    L[i + 1] = _env_l(L[i], A[i], Ws[i])
                    for p, ph in enumerate(prev):
                        OL[p][i + 1] = np.tensordot(A[i], np.tensordot(OL[p][i], ph.A[i], axes=(1, 0)), axes=([0, 1], [0, 1]))
                else:
                    A[i] = (U * S).reshape(shp[0], 2, keep); A[i + 1] = Vt.reshape(keep, 2, shp[3])
                    R[i] = _env_r(R[i + 1], A[i + 1], Ws[i + 1])
                    for p, ph in enumerate(prev):
                        OR[p][i] = np.tensordot(A[i + 1], np.tensordot(ph.A[i + 1], OR[p][i + 1], axes=(2, 1)), axes=([1, 2], [1, 2]))
        hist.append(E)
        if verbose:
            print(f"  sweep {sw}: E = {E:.10f}  chi = {max(a.shape[2] for a in A)}", flush=True)
        if sw >= 3 and chi_sw == chi and abs(hist[-1] - hist[-2]) < tol and abs(hist[-2] - hist[-3]) < tol:
            break
    psi.sweeps = sw + 1; psi.chi = max(a.shape[2] for a in A)
    return E, psi


def _env_l(L, A, W):
    T = np.tensordot(L, A, axes=(2, 0))                  # b w t a'
    T = np.tensordot(T, W, axes=([1, 2], [0, 3]))        # b a' w' s
    return np.tensordot(A, T, axes=([0, 1], [0, 3])).transpose(0, 2, 1)   # b' w' a'
def _env_r(R, A, W):
    T = np.tensordot(A, R, axes=(2, 2))                  # a' t b w
    T = np.tensordot(T, W, axes=([1, 3], [3, 1]))        # a' b w' s
    return np.tensordot(A, T, axes=([1, 2], [3, 1])).transpose(0, 2, 1)   # b' w' a'


def entropies(psi):
    out = []
    for S in psi.S:
        w = S[S > 1e-15] ** 2
        out.append(float(-(w * np.log(w)).sum()))
    return out


def field(psi, eps0):
    """<L_n> on links n = 1..N-1 (state must be right-canonical from the last sweep: centre at site 0)."""
    A = psi.A; N = len(A); out = []
    # left-to-right: carry norm matrix and accumulated charge expectation
    rho = np.ones((1, 1)); acc = np.zeros((1, 1))
    for k in range(N):
        sg = (-1.0) ** (k + 1); q = 0.5 * (SZ + sg * I2)
        Tn = np.tensordot(rho, A[k], axes=(1, 0)); Ta = np.tensordot(acc, A[k], axes=(1, 0))
        new_rho = np.tensordot(A[k], Tn, axes=([0, 1], [0, 1]))
        Tq = np.tensordot(Tn, q, axes=(1, 1)).transpose(0, 2, 1)
        acc = np.tensordot(A[k], Ta, axes=([0, 1], [0, 1])) + np.tensordot(A[k], Tq, axes=([0, 1], [0, 1]))
        rho = new_rho
        if k < N - 1:
            # close with the right part: state is right-canonical beyond site k only if centre is left of k;
            # use explicit right norm environment instead
            out.append((acc, rho))
    # right norm environments
    Rn = [None] * N; Rn[N - 1] = np.ones((1, 1))
    for k in range(N - 1, 0, -1):
        Rn[k - 1] = np.tensordot(A[k], np.tensordot(A[k], Rn[k], axes=(2, 1)), axes=([1, 2], [1, 2]))
    res = []
    for k, (acc, rho) in enumerate(out):
        res.append(eps0 + float(np.sum(acc * Rn[k])) / float(np.sum(rho * Rn[k])))
    return res


def spectrum(N, x, m, eps0=0.0, k=3, Q0=0, chi=120, lam=4.0, shift=True, verbose=False, full=False,
             init_states=None, **kw):
    """k lowest energies in charge sector Q0, in units of g (same convention as qed2.spectrum)."""
    Ws = mpo(N, x, m, eps0, shift, lam, Q0)
    Es, states = [], []
    for j in range(k):
        if init_states is not None:        # warm start from the states of a neighbouring parameter point
            E, psi = solve(Ws, prev=states, chi=chi, seed=j, init=init_states[j], noise=1e-3, verbose=verbose, **kw)
        else:
            E, psi = solve(Ws, prev=states, chi=chi, seed=j, init=states[-1] if states else None, verbose=verbose, **kw)
        Es.append(E); states.append(psi)
    order = np.argsort(Es)
    E = np.array(Es)[order] / (2 * np.sqrt(x))
    return (E, [states[i] for i in order]) if full else E
