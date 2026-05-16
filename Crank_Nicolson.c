#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "prob.h"
#include "umfpk.h"
#include "Option_Characteristic.h"
#include "Crank_Nicolson.h"
#include "greeks.h"

/* ──────────────────────────────────────────────────────────────────────────
 * Static helper: run the backward time-loop for ONE option type.
 *
 * Uses the already-factorised Numeric object and the prebuilt B matrix.
 * Only 2 time-slice arrays (ping-pong) are allocated here.
 *
 * After the loop:
 *   V_curr  = V[t=0]  (option price surface today)
 *   V_pingp = V[t=1]  (one dt into the future, needed for theta)
 * ─────────────────────────────────────────────────────────────────────── */
static int solve_one_type(
    int N_S, int N_v, int N_t, int n,
    double *S_grid, double *v_grid,
    double r, double K, double tau,
    double kappa, double theta_p, double xi, double rho,
    double feller, int feller_ok,
    int option_type, int option_region,
    int *ia_B, int *ja_B, double *a_B,
    int *ia_A, int *ja_A, double *a_A,
    void *Numeric,
    double dt, hsize_t n_cols,
    float *block)
{
    /* ── Ping-pong buffers ──────────────────────────────────────────── */
    double *V_curr = malloc(n * sizeof(double));
    double *V_next = malloc(n * sizeof(double));
    double *rhs    = malloc(n * sizeof(double));

    if (!V_curr || !V_next || !rhs) {
        printf("solve_one_type: allocation failure\n");
        free(V_curr); free(V_next); free(rhs);
        return 1;
    }

    /* ── Initial condition: payoff at maturity (τ=0) ─────────────────── */
    for (int iy = 0; iy < N_v; iy++) {
        for (int ix = 0; ix < N_S; ix++) {
            int ind = ix + N_S * iy;
            double S = S_grid[ix];
            V_curr[ind] = (option_type == CALL)
                         ? fmax(S - K, 0.0)
                         : fmax(K - S, 0.0);
        }
    }

    /* ── Backward time loop ──────────────────────────────────────────── */
    /* V_pingp saves V[t=1] for theta computation after the loop. */
    double *V_pingp = malloc(n * sizeof(double));
    if (!V_pingp) { free(V_curr); free(V_next); free(rhs); return 1; }
    memcpy(V_pingp, V_curr, n * sizeof(double));   /* fallback if N_t=1 */

    for (int t = N_t - 1; t >= 0; t--) {
        double tau_t = (N_t - t) * dt;   /* time-to-expiry at this step */

        /* 1. rhs = B * V_curr */
        for (int i = 0; i < n; i++) {
            rhs[i] = 0.0;
            for (int k = ia_B[i]; k < ia_B[i+1]; k++)
                rhs[i] += a_B[k] * V_curr[ja_B[k]];
        }

        /* 2. Dirichlet boundary corrections on rhs */

        /* West  S=S_min (ix=0) */
        for (int iy = 0; iy < N_v; iy++) {
            int ind = N_S * iy;
            rhs[ind] = (option_type == CALL) ? 0.0 : K * exp(-r * tau_t);
        }
        /* East  S=S_max (ix=N_S-1) */
        for (int iy = 0; iy < N_v; iy++) {
            int ind = (N_S-1) + N_S * iy;
            double S_max = S_grid[N_S-1];
            rhs[ind] = (option_type == CALL)
                      ? S_max - K * exp(-r * tau_t)
                      : 0.0;
        }
        /* North v=v_max (iy=N_v-1), interior ix only */
        for (int ix = 1; ix < N_S-1; ix++) {
            int ind = ix + N_S * (N_v-1);
            double S  = S_grid[ix];
            rhs[ind]  = (option_type == CALL)
                       ? fmax(S - K * exp(-r * tau_t), 0.0)
                       : 0.0;
        }

        /* 3. Solve A * V_next = rhs */
        if (solve_umfpack(n, ia_A, ja_A, a_A, rhs, V_next, Numeric))
            { free(V_curr); free(V_next); free(rhs); free(V_pingp); return 1; }

        /* 4. American early-exercise constraint */
        if (option_region == AMERICAN) {
            for (int i = 0; i < n; i++) {
                int ix  = i % N_S;
                double S = S_grid[ix];
                double intrinsic = (option_type == CALL)
                                  ? fmax(S - K, 0.0)
                                  : fmax(K - S, 0.0);
                V_next[i] = fmax(V_next[i], intrinsic);
            }
        }

        /* 5. Ping-pong: save V[t+1]=V_curr before overwriting */
        if (t == 0) {
            /* V_pingp will hold V[t=1] for theta */
            memcpy(V_pingp, V_curr, n * sizeof(double));
        }

        /* Swap pointers */
        double *tmp = V_curr;
        V_curr = V_next;
        V_next = tmp;
    }
    /* After loop: V_curr = V[t=0], V_pingp = V[t=1] */

    /* ── Greeks ─────────────────────────────────────────────────────── */
    double *gd = malloc(n*sizeof(double));  /* delta */
    double *gg = malloc(n*sizeof(double));  /* gamma */
    double *gt = malloc(n*sizeof(double));  /* theta */
    double *gv = malloc(n*sizeof(double));  /* vega  */
    double *gva= malloc(n*sizeof(double));  /* vanna */
    double *gvo= malloc(n*sizeof(double));  /* volga */

    if (!gd||!gg||!gt||!gv||!gva||!gvo) {
        printf("solve_one_type: greek allocation failure\n");
        free(V_curr); free(V_next); free(rhs); free(V_pingp);
        free(gd); free(gg); free(gt); free(gv); free(gva); free(gvo);
        return 1;
    }

    greeks(V_curr, V_pingp, N_S, N_v, S_grid, v_grid, dt,
           gd, gg, gt, gv, gva, gvo);

    /* ── Fill output block (N_S*N_v rows × n_cols cols) ─────────────── */
    hsize_t row = 0;
    for (int iy = 0; iy < N_v; iy++) {
        for (int ix = 0; ix < N_S; ix++) {
            int   ind   = ix + N_S * iy;
            float *rptr = block + row * n_cols;

            rptr[0]  = (float)kappa;
            rptr[1]  = (float)theta_p;
            rptr[2]  = (float)xi;
            rptr[3]  = (float)rho;
            rptr[4]  = (float)r;
            rptr[5]  = (float)tau;
            rptr[6]  = (float)K;
            rptr[7]  = (float)feller;
            rptr[8]  = (float)feller_ok;
            rptr[9]  = (float)S_grid[ix];   /* actual S, not index */
            rptr[10] = (float)v_grid[iy];   /* actual v, not index */
            rptr[11] = fmaxf((float)V_curr[ind], 0.0f);
            rptr[12] = (float)gd[ind];
            rptr[13] = (float)gg[ind];
            rptr[14] = (float)gt[ind];
            rptr[15] = (float)gv[ind];
            rptr[16] = (float)gva[ind];
            rptr[17] = (float)gvo[ind];
            row++;
        }
    }

    free(V_curr); free(V_next); free(rhs); free(V_pingp);
    free(gd); free(gg); free(gt); free(gv); free(gva); free(gvo);
    return 0;
}

/* ──────────────────────────────────────────────────────────────────────────
 * heston_solve_all
 *
 * Shared setup (L, A, B, UMFPACK factorisation) for all 4 option types.
 * ─────────────────────────────────────────────────────────────────────── */
int heston_solve_all(
    int N_S, int N_v, int N_t,
    double *S_grid, double *v_grid,
    double r, double kappa, double theta, double xi, double rho,
    double K, double tau, double feller, int feller_ok,
    float *blk_eu_call, float *blk_eu_put,
    float *blk_am_call, float *blk_am_put,
    hsize_t n_cols)
{
    double dt = tau / N_t;

    /* ── Step 1: build L (CSR) ─────────────────────────────────────── */
    int n, *ia, *ja;
    double *a;
    if (prob(N_S, N_v, &n, &ia, &ja, &a,
             S_grid, v_grid, r, kappa, theta, xi, rho))
        return 1;

    int nnz = ia[n];

    /* ── Step 2: build A = I - (dt/2)*L  and  B = I + (dt/2)*L ──── */
    int    *ia_A = malloc((n+1)*sizeof(int));
    int    *ja_A = malloc(nnz  *sizeof(int));
    double *a_A  = malloc(nnz  *sizeof(double));
    int    *ia_B = malloc((n+1)*sizeof(int));
    int    *ja_B = malloc(nnz  *sizeof(int));
    double *a_B  = malloc(nnz  *sizeof(double));

    if (!ia_A||!ja_A||!a_A||!ia_B||!ja_B||!a_B) {
        printf("heston_solve_all: A/B allocation failure\n");
        return 1;
    }

    memcpy(ia_A, ia, (n+1)*sizeof(int));
    memcpy(ja_A, ja,  nnz *sizeof(int));
    memcpy(ia_B, ia, (n+1)*sizeof(int));
    memcpy(ja_B, ja,  nnz *sizeof(int));

    for (int k = 0; k < nnz; k++) {
        a_A[k] = -dt * 0.5 * a[k];   /* A = -dt/2 * L */
        a_B[k] = +dt * 0.5 * a[k];   /* B = +dt/2 * L */
    }
    /* Add identity to both */
    for (int i = 0; i < n; i++) {
        for (int k = ia_A[i]; k < ia_A[i+1]; k++) {
            if (ja_A[k] == i) { a_A[k] += 1.0; a_B[k] += 1.0; break; }
        }
    }

    /* ── Step 3: enforce Dirichlet rows in A and B ────────────────── */
    for (int iy = 0; iy < N_v; iy++) {
        for (int ix = 0; ix < N_S; ix++) {
            if (ix == 0 || ix == N_S-1 || iy == N_v-1) {
                int ind = ix + N_S * iy;
                for (int k = ia_A[ind]; k < ia_A[ind+1]; k++) {
                    a_A[k] = (ja_A[k] == ind) ? 1.0 : 0.0;
                    a_B[k] = 0.0;
                }
            }
        }
    }

    /* ── Step 4: factorize A — done ONCE, reused 4 times ─────────── */
    void *Numeric;
    if (factor_umfpack(n, ia_A, ja_A, a_A, &Numeric)) {
        free(ia); free(ja); free(a);
        free(ia_A); free(ja_A); free(a_A);
        free(ia_B); free(ja_B); free(a_B);
        return 1;
    }

    /* ── Step 5: solve for all 4 option types ─────────────────────── */
    int ret = 0;

    ret |= solve_one_type(N_S, N_v, N_t, n, S_grid, v_grid,
                          r, K, tau, kappa, theta, xi, rho,
                          feller, feller_ok, CALL, EUROPEAN,
                          ia_B, ja_B, a_B, ia_A, ja_A, a_A,
                          Numeric, dt, n_cols, blk_eu_call);

    ret |= solve_one_type(N_S, N_v, N_t, n, S_grid, v_grid,
                          r, K, tau, kappa, theta, xi, rho,
                          feller, feller_ok, PUT, EUROPEAN,
                          ia_B, ja_B, a_B, ia_A, ja_A, a_A,
                          Numeric, dt, n_cols, blk_eu_put);

    ret |= solve_one_type(N_S, N_v, N_t, n, S_grid, v_grid,
                          r, K, tau, kappa, theta, xi, rho,
                          feller, feller_ok, CALL, AMERICAN,
                          ia_B, ja_B, a_B, ia_A, ja_A, a_A,
                          Numeric, dt, n_cols, blk_am_call);

    ret |= solve_one_type(N_S, N_v, N_t, n, S_grid, v_grid,
                          r, K, tau, kappa, theta, xi, rho,
                          feller, feller_ok, PUT, AMERICAN,
                          ia_B, ja_B, a_B, ia_A, ja_A, a_A,
                          Numeric, dt, n_cols, blk_am_put);

    /* ── Cleanup ──────────────────────────────────────────────────── */
    free_umfpack(&Numeric);
    free(ia);  free(ja);  free(a);
    free(ia_A); free(ja_A); free(a_A);
    free(ia_B); free(ja_B); free(a_B);

    return ret;
}