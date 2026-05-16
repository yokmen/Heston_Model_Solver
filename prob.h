#ifndef PROB_H
#define PROB_H

/*
 * Build the CSR sparse matrix for the Heston operator L.
 * S_grid[0..N_S-1] and v_grid[0..N_v-1] are non-uniform node arrays.
 * Allocates *ia, *ja, *a internally (caller must free).
 */
int prob(int N_S, int N_v, int *n, int **ia, int **ja, double **a,
         double *S_grid, double *v_grid,
         double r, double kappa, double theta, double xi, double rho);

#endif