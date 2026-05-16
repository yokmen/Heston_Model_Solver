#ifndef HDF5_WRITER_H
#define HDF5_WRITER_H

#include <hdf5.h>

/* One handle per option type */
typedef struct {
    hid_t   file_id;
    hid_t   dset_eu_call;
    hid_t   dset_eu_put;
    hid_t   dset_am_call;
    hid_t   dset_am_put;
    hsize_t n_cols;
    hsize_t row_eu_call;
    hsize_t row_eu_put;
    hsize_t row_am_call;
    hsize_t row_am_put;
} HDF5Writer;

int  hdf5_open (HDF5Writer *w, const char *filename,
                int N_S, int N_v,
                double S_min, double S_max, double v_min,
                double *S_grid,
                hsize_t n_cols);
int  hdf5_write_block(HDF5Writer *w,
                      int option_type, int option_region,
                      float *block, hsize_t n_rows);
void hdf5_close(HDF5Writer *w);

#endif