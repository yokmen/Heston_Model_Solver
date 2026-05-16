#ifndef DATASET_H
#define DATASET_H

/*
 * Generate the full Heston option dataset and write it to HDF5.
 *
 * alpha_S     : sinh concentration for S grid around K (e.g. K/5 = 0.2)
 * alpha_v_rel : sinh concentration for v grid as fraction of v_max (e.g. 0.07)
 *               → alpha_v = alpha_v_rel * v_max
 */
void dataset(int N_S, int N_v, int N_t,
             double S_min, double S_max, double v_min, double K,
             double alpha_S, double alpha_v_rel);

#endif