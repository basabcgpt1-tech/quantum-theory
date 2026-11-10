// y = diag*v + x * sum over bonds of v[flipped state]; tab[i*nb+b] = index or -1
#include <stdint.h>
void matvec(long n, int nb, const int32_t *tab, const double *diag, double x,
            const double *v, double *y) {
#pragma omp parallel for schedule(static)
    for (long i = 0; i < n; i++) {
        const int32_t *t = tab + i * nb;
        double s = 0.0;
        for (int b = 0; b < nb; b++) if (t[b] >= 0) s += v[t[b]];
        y[i] = diag[i] * v[i] + x * s;
    }
}
