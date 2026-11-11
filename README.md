# quantum-theory

Mass gap, confinement and the θ = π critical point of 1+1D QED (Schwinger model):
derivations, thirteen quantitative predictions, and tests by exact diagonalisation and DMRG.

- **Paper:** [paper/paper.pdf](paper/paper.pdf) (26 pages; source `paper/paper.md`)
- **One-page summary:** [THEORY.md](THEORY.md)
- **Solvers:** `qed2.py` + `kernel.c` (exact diagonalisation), `dmrg.py` (DMRG)
- **Calculations:** `run_*.py`; **analysis:** `analyze_*.py`, `check_*.py`; **figures:** `make_figures.py`
- **Data:** `data/` (JSON)

Setup: `gcc -O3 -march=native -fopenmp -shared -fPIC kernel.c -o kernel.so`; Python needs numpy and matplotlib.
Paper build: `npm install katex` once in `paper/`, then `python3 paper/build.py` (needs pandoc and Playwright).
Tests: `python3 test_qed2.py && python3 test_dmrg.py`.
