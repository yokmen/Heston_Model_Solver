#ifndef GRID_H
#define GRID_H

/*
 * build_S_grid  —  Tavella-Randall sinh mapping, concentrates near S = K.
 *
 *   S[i] = K + alpha_S * sinh( c1 + (c2-c1) * i/(N_S-1) )
 *   c1   = asinh( (S_min - K) / alpha_S )
 *   c2   = asinh( (S_max - K) / alpha_S )
 *
 *   alpha_S  controls density: small alpha_S → very tight around K.
 *   Rule of thumb: alpha_S = K / 5 gives ~3× more points near ATM
 *   than a uniform grid.
 */
void build_S_grid(double *S, int N_S,
                  double S_min, double S_max,
                  double K,     double alpha_S);

/*
 * build_v_grid  —  sinh mapping, concentrates near v = 0.
 *
 *   v[j] = alpha_v * sinh( c2 * j/(N_v-1) )
 *   c2   = asinh( v_max / alpha_v )
 *
 *   alpha_v  controls density near zero: alpha_v = v_max/10 gives
 *   roughly 5× more points in the lower 20% of [0, v_max].
 */
void build_v_grid(double *v, int N_v,
                  double v_max, double alpha_v);

#endif