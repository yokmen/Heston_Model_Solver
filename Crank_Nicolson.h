#ifndef CRANK_NICOLSON_H
#define CRANK_NICOLSON_H

#include "hdf5_writer.h"

/*
 * heston_solve_all
 *
 * Solves the Heston PDE for all 4 option types in one call:
 *   European Call, European Put, American Call, American Put
 *
 * The UMFPACK factorisation of A = I - (dt/2)*L is shared across all four,
 * giving a ~4x speedup over calling Crank_Nicolson independently for each.
 *
 * Memory: only 2 time-slice arrays are kept alive at any point (ping-pong).
 *
 * Outputs: blk_eu_call, blk_eu_put, blk_am_call, blk_am_put must each point
 * to a buffer of N_S * N_v * n_cols floats, allocated by the caller.
 *
 * Returns 0 on success.
 */
int heston_solve_all(
    int N_S, int N_v, int N_t,
    double *S_grid, double *v_grid,
    double r, double kappa, double theta, double xi, double rho,
    double K, double tau, double feller, int feller_ok,
    float *blk_eu_call, float *blk_eu_put,
    float *blk_am_call, float *blk_am_put,
    hsize_t n_cols
);

#endif