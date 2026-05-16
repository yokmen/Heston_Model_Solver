#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Crank_Nicolson.h"
#include "Option_Characteristic.h"
#include "dataset.h"
#include "hdf5_writer.h"
#include "grid.h"
#include <omp.h>

/* ─── Parameter grid ──────────────────────────────────────────────────── */
static const double kappas[] = {0.5, 1.0, 1.5, 2.0, 3.0};
static const double thetas[] = {0.02, 0.04, 0.06, 0.08};
static const double xis[]    = {0.1, 0.2, 0.3, 0.4, 0.5};
static const double rhos[]   = {-0.7, -0.5, -0.3, 0.0};
static const double taus[]   = {0.25, 0.5, 1.0, 2.0};
static const double rs[]     = {0.01, 0.02, 0.03, 0.05, 0.07};

static const int n_kappa = 5, n_theta = 4, n_xi = 5,
                 n_rho   = 4, n_tau   = 4, n_r  = 5;


/* ─── Flat parameter combination struct for OpenMP flat parallelism ────── */
typedef struct {
    int ik, it, ix, ir, itau, ir_;
} ParamIdx;

void dataset(int N_S, int N_v, int N_t,
             double S_min, double S_max, double v_min, double K,
             double alpha_S, double alpha_v_rel)
{
    int n_total = n_kappa * n_theta * n_xi * n_rho * n_tau * n_r;

    /* ── Pre-enumerate all parameter combinations ─────────────────────── */
    ParamIdx *params = malloc(n_total * sizeof(ParamIdx));
    if (!params) { printf("dataset: malloc failed\n"); return; }

    int n_params = 0;
    for (int ik=0;ik<n_kappa;ik++)
    for (int it=0;it<n_theta;it++)
    for (int ix=0;ix<n_xi;  ix++)
    for (int ir=0;ir<n_rho; ir++)
    for (int itau=0;itau<n_tau;itau++)
    for (int ir_=0;ir_<n_r; ir_++)
        params[n_params++] = (ParamIdx){ik,it,ix,ir,itau,ir_};

    /* ── Pre-build S grid (fixed: same K, same S_min/S_max) ──────────── */
    double *S_grid = malloc(N_S * sizeof(double));
    if (!S_grid) { printf("dataset: S_grid malloc failed\n"); return; }
    build_S_grid(S_grid, N_S, S_min, S_max, K, alpha_S);

    /* ── Open HDF5 ───────────────────────────────────────────────────── */
    HDF5Writer writer;
    hsize_t n_cols = 18;
    if (hdf5_open(&writer, "dataset_Options.h5",
                  N_S, N_v, S_min, S_max, v_min, S_grid, n_cols)) {
        printf("dataset: HDF5 open failed\n");
        free(params); free(S_grid);
        return;
    }

    hsize_t n_rows = (hsize_t)N_S * N_v;
    int written = 0, skipped_num = 0;
    
    double t_start = 0.0;
#ifdef _OPENMP
    t_start = omp_get_wtime();
#endif

    printf("\n" );
    printf("  Parameter grid:\n");
    printf("    kappa: 5 values  | theta: 4 values | xi: 5 values\n");
    printf("    rho:   4 values  | tau:   4 values | r:  5 values\n");
    printf("\n");
    printf("  Total parameter sets : %d\n", n_total);
    printf("  Grid size per set    : %d × %d = %llu points\n", N_S, N_v, (unsigned long long)n_rows);
    printf("  (Feller < 1.0 sets are included but flagged)\n");
    printf("  Compiler: -O3 optimization enabled for max performance\n\n");

    /* ─────────────────────────────────────────────────────────────────────
     * OpenMP parallel loop over parameter combinations.
     *
     * Key optimizations:
     *  - Threads compute heston_solve_all fully independently (CPU-bound, parallel)
     *  - Results written in critical section (I/O-bound, serialized by HDF5)
     *  - Compiler optimization (-O3) is essential for performance
     *  - schedule(dynamic,1) handles variance in per-set Feller cost
     *
     * With -O3 optimization, ~2.5-3x speedup vs -O0 is typical for this workload
     * ──────────────────────────────────────────────────────────────────── */

#pragma omp parallel for schedule(dynamic,1) reduction(+:skipped_num)
    for (int p = 0; p < n_params; p++) {

        ParamIdx idx = params[p];
        double kappa = kappas[idx.ik];
        double theta = thetas[idx.it];
        double xi    = xis[idx.ix];
        double rho   = rhos[idx.ir];
        double tau   = taus[idx.itau];
        double r     = rs[idx.ir_];

        double feller    = 2.0 * kappa * theta / (xi * xi);
        int    feller_ok = (feller >= 1.0) ? 1 : 0;

        if (feller < 0.2) {
            skipped_num++;
            continue;
        }

        double v_max   = 8.0 * theta;
        double alpha_v = alpha_v_rel * v_max;

        double *v_grid      = malloc(N_v * sizeof(double));
        float  *blk_eu_call = malloc(n_rows * n_cols * sizeof(float));
        float  *blk_eu_put  = malloc(n_rows * n_cols * sizeof(float));
        float  *blk_am_call = malloc(n_rows * n_cols * sizeof(float));
        float  *blk_am_put  = malloc(n_rows * n_cols * sizeof(float));

        if (!v_grid||!blk_eu_call||!blk_eu_put||!blk_am_call||!blk_am_put) {
            printf("Thread %d: allocation failed for p=%d\n",
#ifdef _OPENMP
                   omp_get_thread_num(),
#else
                   0,
#endif
                   p);
            free(v_grid); free(blk_eu_call); free(blk_eu_put);
            free(blk_am_call); free(blk_am_put);
            continue;
        }

        build_v_grid(v_grid, N_v, v_max, alpha_v);

        int ok = heston_solve_all(
            N_S, N_v, N_t,
            S_grid, v_grid,
            r, kappa, theta, xi, rho,
            K, tau, feller, feller_ok,
            blk_eu_call, blk_eu_put,
            blk_am_call, blk_am_put,
            n_cols);

        if (!ok) {
            /* Write results immediately - HDF5 thread-safety via critical section.
               To minimize critical section time, we write all 4 types at once. */
#pragma omp critical(hdf5_write)
            {
                hdf5_write_block(&writer, CALL, EUROPEAN, blk_eu_call, n_rows);
                hdf5_write_block(&writer, PUT,  EUROPEAN, blk_eu_put,  n_rows);
                hdf5_write_block(&writer, CALL, AMERICAN, blk_am_call, n_rows);
                hdf5_write_block(&writer, PUT,  AMERICAN, blk_am_put,  n_rows);
                written++;

                if (written % 50 == 0 || written == n_params) {
                    double elapsed = 0.0;
#ifdef _OPENMP
                    elapsed = omp_get_wtime() - t_start;
#endif
                    double rate = (elapsed > 0) ? (written / elapsed) : 0;
                    double remaining = (rate > 0) ? ((n_params - written) / rate) : 0;
                    
                    int hours = (int)(elapsed / 3600);
                    int mins = (int)((elapsed - hours*3600) / 60);
                    int secs = (int)(elapsed - hours*3600 - mins*60);
                    
                    int eta_hours = (int)(remaining / 3600);
                    int eta_mins = (int)((remaining - eta_hours*3600) / 60);
                    int eta_secs = (int)(remaining - eta_hours*3600 - eta_mins*60);
                    
                    printf("  [%4d/%4d] | Elapsed: %02d:%02d:%02d | Rate: %.2f sets/s | ETA: %02d:%02d:%02d\n",
                           written, n_params, hours, mins, secs, rate, eta_hours, eta_mins, eta_secs);
                }
            }
        }

        free(v_grid); free(blk_eu_call); free(blk_eu_put);
        free(blk_am_call); free(blk_am_put);
    }

    free(params);
    free(S_grid);
    hdf5_close(&writer);

    printf("\nDataset complete.\n");
    printf("  Total param sets : %d\n", n_total);
    printf("  Skipped (Feller<0.2) : %d\n", skipped_num);
    printf("  Written : %d  (%d surfaces × 4 option types)\n",
           written * 4, written);
}