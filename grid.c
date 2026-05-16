#include <math.h>
#include "grid.h"

/* Tavella-Randall mapping: concentrates near S = K */
void build_S_grid(double *S, int N_S,
                  double S_min, double S_max,
                  double K,     double alpha_S)
{
    double c1 = asinh((S_min - K) / alpha_S);
    double c2 = asinh((S_max - K) / alpha_S);

    for (int i = 0; i < N_S; i++) {
        double xi = (double)i / (N_S - 1);
        S[i] = K + alpha_S * sinh(c1 + (c2 - c1) * xi);
    }
    /* Pin exact endpoints to avoid floating-point drift */
    S[0]     = S_min;
    S[N_S-1] = S_max;
}

/* Sinh mapping: concentrates near v = 0 */
void build_v_grid(double *v, int N_v,
                  double v_max, double alpha_v)
{
    double c2 = asinh(v_max / alpha_v);

    for (int j = 0; j < N_v; j++) {
        double eta = (double)j / (N_v - 1);
        v[j] = alpha_v * sinh(c2 * eta);
    }
    /* Pin exact endpoints */
    v[0]     = 0.0;
    v[N_v-1] = v_max;
}