#include <stdio.h>
#include "grid_index.h"
#include "greeks.h"

/*
 * greeks  —  Compute all six Greeks from the option-price surface V[0..n-1]
 *            on a NON-UNIFORM grid (S_grid, v_grid).
 *
 * All spatial derivatives use 2nd-order non-uniform stencils.
 * Boundary nodes use 1st-order one-sided differences.
 *
 * Non-uniform centered 1st derivative at interior node ix:
 *   dV/dS = V_{-1}*(-h_Sp)/(h_Sm*dS_c) + V_0*(h_Sp-h_Sm)/(h_Sm*h_Sp)
 *           + V_{+1}*(h_Sm)/(h_Sp*dS_c)
 *
 * Non-uniform 2nd derivative:
 *   d²V/dS² = 2*[V_{-1}/(h_Sm*dS_c) - V_0/(h_Sm*h_Sp) + V_{+1}/(h_Sp*dS_c)]
 *
 * Non-uniform mixed derivative:
 *   d²V/dSdv = [V_{++} - V_{-+} - V_{+-} + V_{--}] / (dS_c * dv_c)
 */
void greeks(double *V, double *V_prev,
            int N_S, int N_v,
            double *S_grid, double *v_grid,
            double dt,
            double *greek_delta, double *greek_gamma,
            double *greek_theta, double *greek_vega,
            double *greek_vanna, double *greek_volga)
{
    for (int iy = 0; iy < N_v; iy++) {
        for (int ix = 0; ix < N_S; ix++) {

            int ind = grid_index(N_S, ix, iy);

            /* ── Delta  ∂V/∂S ─────────────────────────────────────── */
            if (ix == 0) {
                double h = S_grid[1] - S_grid[0];
                greek_delta[ind] = (V[grid_index(N_S,1,iy)] - V[ind]) / h;
            }
            else if (ix == N_S-1) {
                double h = S_grid[N_S-1] - S_grid[N_S-2];
                greek_delta[ind] = (V[ind] - V[grid_index(N_S,N_S-2,iy)]) / h;
            }
            else {
                double h_Sm = S_grid[ix] - S_grid[ix-1];
                double h_Sp = S_grid[ix+1] - S_grid[ix];
                double dS_c = h_Sm + h_Sp;
                greek_delta[ind] =
                    V[grid_index(N_S,ix-1,iy)] * (-h_Sp) / (h_Sm * dS_c)
                  + V[ind]                     * (h_Sp - h_Sm) / (h_Sm * h_Sp)
                  + V[grid_index(N_S,ix+1,iy)] * ( h_Sm) / (h_Sp * dS_c);
            }

            /* ── Gamma  ∂²V/∂S² ───────────────────────────────────── */
            if (ix == 0 || ix == N_S-1) {
                greek_gamma[ind] = 0.0;
            }
            else {
                double h_Sm = S_grid[ix]   - S_grid[ix-1];
                double h_Sp = S_grid[ix+1] - S_grid[ix];
                double dS_c = h_Sm + h_Sp;
                greek_gamma[ind] = 2.0 * (
                    V[grid_index(N_S,ix-1,iy)] / (h_Sm * dS_c)
                  - V[ind]                     / (h_Sm * h_Sp)
                  + V[grid_index(N_S,ix+1,iy)] / (h_Sp * dS_c)
                );
            }

            /* ── Theta  -∂V/∂τ  (calendar-time convention) ──────────
                * V_prev is V at one time-step closer to expiry (smaller τ = (N_t-1)*dt),
                * used to approximate -∂V/∂τ ≈ -(V_curr - V_prev)/dt.               
                * θ_finance = -∂V/∂τ ≈ -(V_curr - V_prev)/dt             
            */
            greek_theta[ind] = -(V[ind] - V_prev[ind]) / dt;

            /* ── Vega  ∂V/∂v ──────────────────────────────────────── */
            if (iy == 0) {
                double h = v_grid[1] - v_grid[0];
                greek_vega[ind] = (V[grid_index(N_S,ix,1)] - V[ind]) / h;
            }
            else if (iy == N_v-1) {
                double h = v_grid[N_v-1] - v_grid[N_v-2];
                greek_vega[ind] = (V[ind] - V[grid_index(N_S,ix,N_v-2)]) / h;
            }
            else {
                double h_vm = v_grid[iy]   - v_grid[iy-1];
                double h_vp = v_grid[iy+1] - v_grid[iy];
                double dv_c = h_vm + h_vp;
                greek_vega[ind] =
                    V[grid_index(N_S,ix,iy-1)] * (-h_vp) / (h_vm * dv_c)
                  + V[ind]                     * (h_vp - h_vm) / (h_vm * h_vp)
                  + V[grid_index(N_S,ix,iy+1)] * ( h_vm) / (h_vp * dv_c);
            }

            /* ── Vanna  ∂²V/∂S∂v ─────────────────────────────────── */
            if (ix == 0 || ix == N_S-1 || iy == 0 || iy == N_v-1) {
                greek_vanna[ind] = 0.0;
            }
            else {
                double dS_c = S_grid[ix+1] - S_grid[ix-1];
                double dv_c = v_grid[iy+1] - v_grid[iy-1];
                greek_vanna[ind] = (
                    V[grid_index(N_S,ix+1,iy+1)]
                  - V[grid_index(N_S,ix-1,iy+1)]
                  - V[grid_index(N_S,ix+1,iy-1)]
                  + V[grid_index(N_S,ix-1,iy-1)]
                ) / (dS_c * dv_c);
            }

            /* ── Volga  ∂²V/∂v² ───────────────────────────────────── */
            if (iy == 0 || iy == N_v-1) {
                greek_volga[ind] = 0.0;
            }
            else {
                double h_vm = v_grid[iy]   - v_grid[iy-1];
                double h_vp = v_grid[iy+1] - v_grid[iy];
                double dv_c = h_vm + h_vp;
                greek_volga[ind] = 2.0 * (
                    V[grid_index(N_S,ix,iy-1)] / (h_vm * dv_c)
                  - V[ind]                     / (h_vm * h_vp)
                  + V[grid_index(N_S,ix,iy+1)] / (h_vp * dv_c)
                );
            }
        }
    }
}