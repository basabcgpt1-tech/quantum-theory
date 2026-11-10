"""Lattice Hamiltonian of 1+1D QED (Schwinger model), exact diagonalisation.

Kogut-Susskind staggered fermions, open boundaries, Gauss law solved exactly:

  W = 2H/(a g^2) = x * sum_n (s+_n s-_{n+1} + h.c.)
                 + (mu/2) * sum_n (1 + (-1)^n sz_n)
                 + sum_{n=1}^{N-1} L_n^2
  L_n = eps0 + 1/2 * sum_{k<=n} (sz_k + (-1)^k),  x = 1/(g a)^2,  mu = 2 m_lat sqrt(x)/g

eps0 = theta/2pi is the background electric field.  Energies in units of g are
W/(2 sqrt(x)); the physical length is L g = N/sqrt(x).
m_lat = m - g^2 a/8 (Dempsey, Klebanov, Pufu, Zan 2022) removes O(a) errors.

Build the kernel once:  gcc -O3 -march=native -fopenmp -shared -fPIC kernel.c -o kernel.so
"""
import ctypes, functools, os
import numpy as np

_lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "kernel.so"))
_dp = np.ctypeslib.ndpointer(np.float64, flags="C")
_lib.matvec.argtypes = [ctypes.c_long, ctypes.c_int, np.ctypeslib.ndpointer(np.int32, flags="C"),
                        _dp, ctypes.c_double, _dp, _dp]
_lib.matvec.restype = None


@functools.lru_cache(maxsize=2)
def tables(N):
    """Zero-charge basis, hop table, and the pieces of the diagonal."""
    allst = np.arange(1 << N, dtype=np.int64)
    pc = np.zeros(1 << N, dtype=np.int8)
    for i in range(N):
        pc += ((allst >> i) & 1).astype(np.int8)
    st = allst[pc == N // 2]
    del allst, pc
    n = len(st)
    tab = np.full((n, N - 1), -1, dtype=np.int32)
    for i in range(N - 1):
        src = np.nonzero(((st >> i) & 1) != ((st >> (i + 1)) & 1))[0]
        tab[src, i] = np.searchsorted(st, st[src] ^ (3 << i))
    c = np.zeros(n)            # running charge: L_n = eps0 + c_n
    A = np.zeros(n); B = np.zeros(n); Mt = np.zeros(n)
    cs = np.zeros((N - 1, n), dtype=np.float32) if N <= 22 else None
    for i in range(N):
        sz = 2.0 * ((st >> i) & 1) - 1.0
        sg = (-1.0) ** (i + 1)
        c += 0.5 * (sz + sg)
        Mt += 0.5 * (1.0 + sg * sz)
        if i < N - 1:
            A += c * c; B += c
            if cs is not None:
                cs[i] = c
    return n, np.ascontiguousarray(tab), A, B, Mt, cs


def spectrum(N, x, m, eps0=0.0, k=2, shift=True, tol=1e-11, vec=False):
    """k lowest eigenvalues of W (Lanczos, full reorthogonalisation).
    m is the physical m/g; returns energies in units of g (and the ground vector)."""
    n, tab, A, B, Mt, _ = tables(N)
    ml = m - 1.0 / (8 * np.sqrt(x)) if shift else m
    diag = A + 2 * eps0 * B + (N - 1) * eps0 ** 2 + 2.0 * ml * np.sqrt(x) * Mt
    maxsteps = int(min(500, n, 6e9 / (8 * n)))
    Q = np.empty((maxsteps, n))
    q = np.random.default_rng(1).standard_normal(n); q /= np.linalg.norm(q)
    w = np.empty(n); al, be = [], []; prev = None
    for j in range(maxsteps):
        Q[j] = q
        _lib.matvec(n, N - 1, tab, diag, x, q, w)
        al.append(q @ w)
        for _ in range(2):
            w -= Q[:j + 1].T @ (Q[:j + 1] @ w)
        b = np.linalg.norm(w)
        m_ = j + 1
        if b < 1e-12 or (m_ >= 30 and m_ % 10 == 0) or m_ == maxsteps:
            T = np.diag(al) + np.diag(be[:m_ - 1], 1) + np.diag(be[:m_ - 1], -1)
            ev, U = np.linalg.eigh(T)
            cur = ev[:k]
            if b < 1e-12 or (prev is not None and len(prev) == len(cur)
                             and np.all(np.abs(cur - prev) < tol * max(1.0, abs(cur[0])))):
                break
            prev = cur
        be.append(b); q = w / b
    else:
        raise RuntimeError("Lanczos not converged")
    E = ev[:k] / (2 * np.sqrt(x))
    if vec:
        return E, U[:, 0] @ Q[:m_]
    return E


def gap(N, x, m, eps0=0.0, **kw):
    e = spectrum(N, x, m, eps0, **kw)
    return e[1] - e[0]
