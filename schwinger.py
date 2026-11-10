"""Mass gap of 1+1D QED (Schwinger model) from the lattice Hamiltonian.

Kogut-Susskind staggered fermions, open boundaries, Gauss law solved exactly,
so the gauge field is eliminated and only N spins remain:

  W = 2H/(a g^2) = x * sum_n (s+_n s-_{n+1} + h.c.)
                 + (mu/2) * sum_n (1 + (-1)^n sz_n)
                 + sum_{n=1}^{N-1} L_n^2
  L_n = eps0 + 1/2 * sum_{k<=n} (sz_k + (-1)^k),   x = 1/(g a)^2,  mu = 2 m sqrt(x)/g

Energies in units of g:  E/g = W / (2 sqrt(x)).
"""
import itertools, json, sys, time
import numpy as np


def basis(N):
    """All bit patterns with N/2 up spins (zero total charge), sorted."""
    s = np.array([sum(1 << i for i in c) for c in itertools.combinations(range(N), N // 2)],
                 dtype=np.int64)
    return np.sort(s)


def build(N, x, m_over_g=0.0, eps0=0.0):
    st = basis(N)
    sz = np.array([2.0 * ((st >> i) & 1) - 1.0 for i in range(N)])  # site k=i+1
    sign = np.array([(-1.0) ** (i + 1) for i in range(N)])
    L = eps0 + 0.5 * np.cumsum(sz + sign[:, None], axis=0)
    mu = 2.0 * m_over_g * np.sqrt(x)
    diag = (L[:-1] ** 2).sum(0) + 0.5 * mu * (1.0 + sign[:, None] * sz).sum(0)
    hops = []
    for i in range(N - 1):
        differ = ((st >> i) & 1) != ((st >> (i + 1)) & 1)
        src = np.nonzero(differ)[0]
        dst = np.searchsorted(st, st[src] ^ (3 << i))
        hops.append((src, dst))

    def matvec(v):
        y = diag * v
        for src, dst in hops:
            y[dst] += x * v[src]
        return y
    return len(st), matvec


def lowest(n, matvec, k=2, steps=160, seed=0):
    """Lanczos with full reorthogonalisation; returns k lowest eigenvalues."""
    steps = min(steps, n)
    rng = np.random.default_rng(seed)
    Q = np.empty((steps, n))
    q = rng.standard_normal(n); q /= np.linalg.norm(q)
    al, be = [], []
    for j in range(steps):
        Q[j] = q
        w = matvec(q)
        a = q @ w; al.append(a)
        w -= Q[:j + 1].T @ (Q[:j + 1] @ w)
        w -= Q[:j + 1].T @ (Q[:j + 1] @ w)
        b = np.linalg.norm(w)
        if b < 1e-10:
            break
        be.append(b); q = w / b
    m = len(al)
    T = np.diag(al) + np.diag(be[:m - 1], 1) + np.diag(be[:m - 1], -1)
    return np.linalg.eigvalsh(T)[:k]


def gap(N, x, m=0.0, eps0=0.0):
    n, mv = build(N, x, m, eps0)
    e = lowest(n, mv)
    return (e[1] - e[0]) / (2 * np.sqrt(x)), e[0] / (2 * np.sqrt(x))


SHIFT = True   # m_lat = m - g^2 a / 8 (Dempsey et al. 2022): removes O(a) errors


def bulk_gap(x, m, Ns):
    """N -> infinity. Relativistic particle in a box: gap^2 = M^2 + c/N^2."""
    ml = m - 1.0 / (8 * np.sqrt(x)) if SHIFT else m
    g = np.array([gap(N, x, ml)[0] for N in Ns])
    A = np.vstack([np.ones(len(Ns)), 1.0 / np.array(Ns, float) ** 2]).T
    M2 = np.linalg.lstsq(A, g ** 2, rcond=None)[0][0]
    return np.sqrt(max(M2, 0.0)), g


def continuum(m, xs, Ns):
    """x -> infinity (a -> 0): polynomial in 1/sqrt(x) = g a."""
    M = np.array([bulk_gap(x, m, Ns)[0] for x in xs])
    ga = 1.0 / np.sqrt(xs)
    lin = np.polyfit(ga ** 2, M, 1)[-1] if SHIFT else np.polyfit(ga, M, 1)[-1]
    quad = np.polyfit(ga, M, 2)[-1]
    return M, lin, quad


if __name__ == "__main__":
    t0 = time.time()
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    Ns = list(range(Nmax - 6, Nmax + 1, 2))
    xs = np.array([1.0, 1.5, 2.0, 3.0, 4.0, 6.0])
    out = {"Ns": Ns, "xs": xs.tolist(), "exact_massless": 1 / np.sqrt(np.pi), "masses": {}}
    for m in (0.0, 0.05, 0.1):
        M, lin, quad = continuum(m, xs, Ns)
        out["masses"][str(m)] = {"M_of_x": M.tolist(), "lin": lin, "quad": quad}
        print(f"m/g={m:5.2f}  M/g(x)={np.round(M, 4)}  a->0: lin={lin:.4f} quad={quad:.4f}"
              f"  [{time.time() - t0:.0f}s]", flush=True)
    # string tension for external charges +-alpha (background field eps0 = alpha)
    out["string"] = {}
    x, N1, N2 = 4.0, Nmax - 4, Nmax
    for m in (0.0, 0.05, 0.1):
        row = {}
        for al in (0.0, 0.25, 0.5):
            ml = m - 1.0 / (8 * np.sqrt(x)) if SHIFT else m
            e1, e2 = gap(N1, x, ml, al)[1], gap(N2, x, ml, al)[1]
            row[al] = (e2 - e1) / ((N2 - N1) / np.sqrt(x))   # bulk energy density / g^2
        sig = {a: row[a] - row[0.0] for a in row}
        out["string"][str(m)] = sig
        print(f"m/g={m}: string tension sigma/g^2 vs alpha: "
              + ", ".join(f"{a}: {s:.4f}" for a, s in sig.items()), flush=True)
    json.dump(out, open("results.json", "w"), indent=1)
    S = np.exp(np.euler_gamma) / (2 * np.pi ** 1.5)
    print(f"first-order prediction sigma/g^2 = (m/g) * {S:.4f} * (1 - cos 2 pi alpha)")
    print(f"exact massless M/g = 1/sqrt(pi) = {1 / np.sqrt(np.pi):.4f}; done in {time.time() - t0:.0f}s")
