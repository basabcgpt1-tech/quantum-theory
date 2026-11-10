# quantum-theory

Mass gap, confinement and the θ = π critical point of 1+1D QED (Schwinger model):
derivations, twelve quantitative predictions, and an exact-diagonalisation test.

- **Paper:** [paper/paper.pdf](paper/paper.pdf) (15 pages; source `paper/paper.md`, build with `python3 paper/build.py`)
- **One-page summary:** [THEORY.md](THEORY.md)
- **Code:** `qed2.py` + `kernel.c` (solver), `run_*.py` (calculations), `analyze_*.py`, `make_figures.py`, `check_formulas.py`
- **Data:** `data/*.json`

Build the kernel first: `gcc -O3 -march=native -fopenmp -shared -fPIC kernel.c -o kernel.so`.
