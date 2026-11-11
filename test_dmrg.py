"""DMRG against exact diagonalisation."""
import time
import numpy as np
import qed2, dmrg

def test_matches_ed():
    for N, x, m, e0, dq in [(12, 0.5, 0.356, 0.5, 0), (16, 1.0, 0.33, 0.5, 0), (12, 4.0, 0.1, 0.25, 0), (12, 0.5, 0.356, 0.5, -1)]:
        ref = qed2.spectrum(N, x, m, e0, k=3, dq=dq)
        t = time.time(); got = dmrg.spectrum(N, x, m, e0, k=3, Q0=dq)
        print(N, x, m, e0, dq, np.abs(got - ref).max(), f"{time.time() - t:.1f}s")
        assert np.allclose(ref, got, atol=2e-8), (ref, got)

if __name__ == "__main__":
    test_matches_ed(); print("ok")
