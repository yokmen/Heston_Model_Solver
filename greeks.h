#ifndef GREEKS_H
#define GREEKS_H

/*
 * Compute all Greeks from the price surface V (at τ=T) and V_prev (at τ=T-dt)
 * on the non-uniform grids S_grid[0..N_S-1], v_grid[0..N_v-1].
 */
void greeks(double *V, double *V_prev,
            int N_S, int N_v,
            double *S_grid, double *v_grid,
            double dt,
            double *greek_delta, double *greek_gamma,
            double *greek_theta, double *greek_vega,
            double *greek_vanna, double *greek_volga);

#endif