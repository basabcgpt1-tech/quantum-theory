"""Checks: Lanczos vs an independent dense Kronecker-product Hamiltonian; exact gap."""
import numpy as np
import qed2


def dense(N, x, ml, eps0):
    """Full 2^N matrix built from Pauli operators, then restricted to zero charge."""
    I = np.eye(2); sz = np.diag([1.0, -1.0]); sp = np.array([[0, 1.0], [0, 0]])

    def op(o, i):
        m = np.array([[1.0]])
        for j in range(N):
            m = np.kron(m, o if j == i else I)
        return m
    Z = [op(sz, i) for i in range(N)]; P = [op(sp, i) for i in range(N)]
    Id = np.eye(2 ** N); H = np.zeros((2 ** N, 2 ** N)); L = eps0 * Id
    for i in range(N):
        sg = (-1.0) ** (i + 1)
        H += ml * np.sqrt(x) * (Id + sg * Z[i])
        L = L + 0.5 * (Z[i] + sg * Id)
        if i < N - 1:
            H += L @ L + x * (P[i] @ P[i + 1].T + P[i + 1] @ P[i].T)
    keep = np.isclose(np.diag(sum(Z)), 0)
    return np.linalg.eigvalsh(H[np.ix_(keep, keep)]) / (2 * np.sqrt(x))


def test_lanczos_matches_dense():
    for x, m, e0 in [(1.0, 0.0, 0.0), (4.0, 0.3, 0.5), (9.0, 0.1, 0.25)]:
        ref = dense(8, x, m, e0)[:4]
        got = qed2.spectrum(8, x, m, e0, k=4, shift=False)
        assert np.allclose(ref, got, atol=1e-9), (ref, got)


def test_massless_gap_near_exact():
    assert abs(qed2.gap(16, 1.5, 0.0) - 1 / np.sqrt(np.pi)) < 0.06


if __name__ == "__main__":
    test_lanczos_matches_dense(); test_massless_gap_near_exact(); print("ok")
