#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "prob.h"
#include "grid_index.h"

/*
 * prob  —  Build the CSR representation of the Heston spatial operator L.
 *
 * Non-uniform grids S_grid[0..N_S-1] and v_grid[0..N_v-1] are passed in.
 * All finite-difference weights are derived from the actual node spacings:
 *
 *   h_Sm = S[ix] - S[ix-1]        h_Sp = S[ix+1] - S[ix]
 *   h_vm = v[iy] - v[iy-1]        h_vp = v[iy+1] - v[iy]
 *
 * Interior 9-point stencil (second-order everywhere):
 *
 *   L = (1/2)vS² d²/dS²  +  ρξvS d²/dSdv  +  (1/2)ξ²v d²/dv²
 *       + rS d/dS  +  κ(θ-v) d/dv  -  rI
 *
 * South boundary (v=0, ix interior): reduced PDE via Feller condition.
 * North, West, East boundaries: Dirichlet identity rows.
 */
int prob(int N_S, int N_v, int *n, int **ia, int **ja, double **a,
         double *S_grid, double *v_grid,
         double r, double kappa, double theta, double xi, double rho)
{
    if (N_S <= 0 || N_v <= 0) {
        printf("Error: N_S=%d N_v=%d\n", N_S, N_v);
        return 1;
    }

    *n = N_S * N_v;

    /* nnz upper bound (same stencil topology as uniform case) */
    int nnz_max = 9*N_S*N_v - 13*N_S - 16*N_v + 28;

    *ia = malloc((*n + 1) * sizeof(int));
    *ja = malloc(nnz_max   * sizeof(int));
    *a  = malloc(nnz_max   * sizeof(double));

    if (!*ia || !*ja || !*a) {
        printf("prob: allocation failure\n");
        return 1;
    }

    int nnz = 0;

    for (int iy = 0; iy < N_v; iy++) {
        for (int ix = 0; ix < N_S; ix++) {

            int ind = grid_index(N_S, ix, iy);
            double Si = S_grid[ix];
            double vi = v_grid[iy];

            (*ia)[ind] = nnz;

            /* ── West border ix=0 : Dirichlet identity ──────────────── */
            if (ix == 0) {
                (*ja)[nnz] = ind; (*a)[nnz] = 1.0; nnz++;
            }
            /* ── East border ix=N_S-1 : Dirichlet identity ───────────── */
            else if (ix == N_S-1) {
                (*ja)[nnz] = ind; (*a)[nnz] = 1.0; nnz++;
            }
            /* ── North border iy=N_v-1, ix interior : Dirichlet ──────── */
            else if (iy == N_v-1) {
                (*ja)[nnz] = ind; (*a)[nnz] = 1.0; nnz++;
            }
            /* ── South border iy=0, ix interior : reduced PDE ───────────
               v=0  →  diffusion and vol-of-vol terms vanish.
               Remaining: rS ∂V/∂S + κθ ∂V/∂v - rV = 0
               ∂V/∂S : upwind backward (r>0, S>0)
               ∂V/∂v : forward difference                                */
            else if (iy == 0) {
                double h_Sm  = S_grid[ix] - S_grid[ix-1];
                double h_vp0 = v_grid[1]  - v_grid[0];   /* first v spacing */

                double adv_S = r * Si   / h_Sm;
                double adv_v = kappa * theta / h_vp0;

                int ind_im = grid_index(N_S, ix-1, 0);
                int ind_ip1= grid_index(N_S, ix,   1);

                (*ja)[nnz] = ind_im;  (*a)[nnz] = -adv_S;          nnz++;
                (*ja)[nnz] = ind;     (*a)[nnz] =  adv_S - adv_v - r; nnz++;
                (*ja)[nnz] = ind_ip1; (*a)[nnz] =  adv_v;          nnz++;
            }
            /* ── Interior : 9-point stencil with non-uniform weights ──── */
            else {
                double h_Sm = S_grid[ix]   - S_grid[ix-1];
                double h_Sp = S_grid[ix+1] - S_grid[ix];
                double h_vm = v_grid[iy]   - v_grid[iy-1];
                double h_vp = v_grid[iy+1] - v_grid[iy];

                double dS_c = h_Sm + h_Sp;   /* = S[ix+1] - S[ix-1] */
                double dv_c = h_vm + h_vp;   /* = v[iy+1] - v[iy-1] */

                /* PDE coefficients */
                double a_SS  = 0.5 * vi * Si * Si;
                double a_vv  = 0.5 * xi * xi * vi;
                double a_S   = r * Si;
                double a_v   = kappa * (theta - vi);
                double a_mix = rho * xi * vi * Si;

                /*
                 * Non-uniform d²f/dS² :
                 *   2*f_{-1}/(h_Sm*dS_c)  - 2*f_0/(h_Sm*h_Sp)  + 2*f_{+1}/(h_Sp*dS_c)
                 *
                 * Non-uniform df/dS (2nd-order) :
                 *   f_{-1}*(-h_Sp)/(h_Sm*dS_c)  + f_0*(h_Sp-h_Sm)/(h_Sm*h_Sp)
                 *   + f_{+1}*(h_Sm)/(h_Sp*dS_c)
                 *
                 * Same formulas in v direction.
                 *
                 * Non-uniform mixed d²f/dSdv :
                 *   [ f_{++} - f_{-+} - f_{+-} + f_{--} ] / (dS_c * dv_c)
                 */

                /* Coefficient of V_{ix-1, iy} */
                double c_im = 2.0*a_SS/(h_Sm*dS_c) + a_S*(-h_Sp)/(h_Sm*dS_c);

                /* Coefficient of V_{ix+1, iy} */
                double c_ip = 2.0*a_SS/(h_Sp*dS_c) + a_S*(h_Sm)/(h_Sp*dS_c);

                /* Coefficient of V_{ix, iy-1} */
                double c_jm = 2.0*a_vv/(h_vm*dv_c) + a_v*(-h_vp)/(h_vm*dv_c);

                /* Coefficient of V_{ix, iy+1} */
                double c_jp = 2.0*a_vv/(h_vp*dv_c) + a_v*(h_vm)/(h_vp*dv_c);

                /* Diagonal */
                double c_diag = -2.0*a_SS/(h_Sm*h_Sp)
                               + a_S*(h_Sp-h_Sm)/(h_Sm*h_Sp)
                               - 2.0*a_vv/(h_vm*h_vp)
                               + a_v*(h_vp-h_vm)/(h_vm*h_vp)
                               - r;

                /* Mixed (cross) weight */
                double c_mix = a_mix / (dS_c * dv_c);

                /* Fill 9 entries in column-index order */
                int ind_mm = grid_index(N_S, ix-1, iy-1);
                int ind_0m = grid_index(N_S, ix,   iy-1);
                int ind_pm = grid_index(N_S, ix+1, iy-1);
                int ind_m0 = grid_index(N_S, ix-1, iy  );
                int ind_p0 = grid_index(N_S, ix+1, iy  );
                int ind_mp = grid_index(N_S, ix-1, iy+1);
                int ind_0p = grid_index(N_S, ix,   iy+1);
                int ind_pp = grid_index(N_S, ix+1, iy+1);

                (*ja)[nnz]=ind_mm; (*a)[nnz]=+c_mix;  nnz++;
                (*ja)[nnz]=ind_0m; (*a)[nnz]=c_jm;    nnz++;
                (*ja)[nnz]=ind_pm; (*a)[nnz]=-c_mix;  nnz++;
                (*ja)[nnz]=ind_m0; (*a)[nnz]=c_im;    nnz++;
                (*ja)[nnz]=ind;    (*a)[nnz]=c_diag;  nnz++;
                (*ja)[nnz]=ind_p0; (*a)[nnz]=c_ip;    nnz++;
                (*ja)[nnz]=ind_mp; (*a)[nnz]=-c_mix;  nnz++;
                (*ja)[nnz]=ind_0p; (*a)[nnz]=c_jp;    nnz++;
                (*ja)[nnz]=ind_pp; (*a)[nnz]=+c_mix;  nnz++;
            }
        }
    }

    (*ia)[*n] = nnz;
    return 0;
}