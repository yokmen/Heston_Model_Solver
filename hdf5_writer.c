#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hdf5_writer.h"
#include "Option_Characteristic.h"

/* Column names stored as HDF5 attribute */
static const char *COL_NAMES =
    "kappa,theta,xi,rho,r,tau,K,feller,feller_ok,"
    "S,v,"
    "prix,delta,gamma,theta_gr,vega,vanna,volga";

/* Create one unlimited dataset inside a group */
static hid_t create_dataset(hid_t file_id, const char *path,
                             hsize_t n_cols)
{
    /* Chunk: 10000 rows at a time — good for sequential write */
    hsize_t chunk[2]   = {10000, n_cols};
    hsize_t maxdims[2] = {H5S_UNLIMITED, n_cols};
    hsize_t dims[2]    = {0, n_cols};

    hid_t dspace = H5Screate_simple(2, dims, maxdims);
    hid_t dcpl   = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 2, chunk);

    /* gzip compression level 4 — good balance speed/size */
    H5Pset_deflate(dcpl, 4);

    hid_t dset = H5Dcreate2(file_id, path,
                             H5T_NATIVE_FLOAT,
                             dspace, H5P_DEFAULT, dcpl,
                             H5P_DEFAULT);

    /* Store column names as string attribute */
    hid_t atype  = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, strlen(COL_NAMES) + 1);
    hid_t aspace = H5Screate(H5S_SCALAR);
    hid_t attr   = H5Acreate2(dset, "columns", atype,
                               aspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, atype, COL_NAMES);
    H5Aclose(attr);
    H5Sclose(aspace);
    H5Tclose(atype);

    H5Pclose(dcpl);
    H5Sclose(dspace);

    return dset;
}

int hdf5_open(HDF5Writer *w, const char *filename,
              int N_S, int N_v,
              double S_min, double S_max, double v_min,
              double *S_grid,
              hsize_t n_cols)
{
    memset(w, 0, sizeof(HDF5Writer));
    w->n_cols = n_cols;

    /* Create file — truncate if exists */
    w->file_id = H5Fcreate(filename, H5F_ACC_TRUNC,
                            H5P_DEFAULT, H5P_DEFAULT);
    if (w->file_id < 0) {
        printf("Error: cannot create %s\n", filename);
        return 1;
    }

    /* Create groups */
    hid_t g_eu = H5Gcreate2(w->file_id, "/european",
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t g_am = H5Gcreate2(w->file_id, "/american",
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Create datasets */
    w->dset_eu_call = create_dataset(w->file_id,
                                      "/european/call", n_cols);
    w->dset_eu_put  = create_dataset(w->file_id,
                                      "/european/put",  n_cols);
    w->dset_am_call = create_dataset(w->file_id,
                                      "/american/call", n_cols);
    w->dset_am_put  = create_dataset(w->file_id,
                                      "/american/put",  n_cols);

    H5Gclose(g_eu);
    H5Gclose(g_am);

    /* Store grid metadata */
    hid_t g_grid = H5Gcreate2(w->file_id, "/grid",
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* S grid */
    float *S_arr = malloc(N_S * sizeof(float));
    for (int i = 0; i < N_S; i++)
        S_arr[i] = (float)S_grid[i];
    hsize_t dims_S = N_S;
    hid_t sp_S  = H5Screate_simple(1, &dims_S, NULL);
    hid_t ds_S  = H5Dcreate2(w->file_id, "/grid/S",
                               H5T_NATIVE_FLOAT, sp_S,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(ds_S, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
             H5P_DEFAULT, S_arr);
    H5Dclose(ds_S);
    H5Sclose(sp_S);
    free(S_arr);

    /* Metadata attributes on root */
    hsize_t scalar = 1;
    hid_t aspace = H5Screate_simple(1, &scalar, NULL);

    int N_S_val = N_S, N_v_val = N_v;
    hid_t a;

    a = H5Acreate2(w->file_id, "N_S", H5T_NATIVE_INT,
                    aspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(a, H5T_NATIVE_INT, &N_S_val); H5Aclose(a);

    a = H5Acreate2(w->file_id, "N_v", H5T_NATIVE_INT,
                    aspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(a, H5T_NATIVE_INT, &N_v_val); H5Aclose(a);

    double S_min_d = S_min, S_max_d = S_max, v_min_d = v_min;
    a = H5Acreate2(w->file_id, "S_min", H5T_NATIVE_DOUBLE,
                    aspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(a, H5T_NATIVE_DOUBLE, &S_min_d); H5Aclose(a);

    a = H5Acreate2(w->file_id, "S_max", H5T_NATIVE_DOUBLE,
                    aspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(a, H5T_NATIVE_DOUBLE, &S_max_d); H5Aclose(a);

    a = H5Acreate2(w->file_id, "v_min", H5T_NATIVE_DOUBLE,
                    aspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(a, H5T_NATIVE_DOUBLE, &v_min_d); H5Aclose(a);

    H5Sclose(aspace);
    H5Gclose(g_grid);

    /* Init row counters */
    w->row_eu_call = 0;
    w->row_eu_put  = 0;
    w->row_am_call = 0;
    w->row_am_put  = 0;

    return 0;
}

int hdf5_write_block(HDF5Writer *w,
                     int option_type, int option_region,
                     float *block, hsize_t n_rows)
{
    hid_t    dset;
    hsize_t *row_counter;

    /* Select target dataset */
    if (option_type == CALL && option_region == EUROPEAN) {
        dset        = w->dset_eu_call;
        row_counter = &w->row_eu_call;
    } else if (option_type == PUT  && option_region == EUROPEAN) {
        dset        = w->dset_eu_put;
        row_counter = &w->row_eu_put;
    } else if (option_type == CALL && option_region == AMERICAN) {
        dset        = w->dset_am_call;
        row_counter = &w->row_am_call;
    } else {
        dset        = w->dset_am_put;
        row_counter = &w->row_am_put;
    }

    /* Extend dataset */
    hsize_t new_dims[2] = {*row_counter + n_rows, w->n_cols};
    H5Dset_extent(dset, new_dims);

    /* Select hyperslab for new rows */
    hid_t   fspace = H5Dget_space(dset);
    hsize_t offset[2] = {*row_counter, 0};
    hsize_t count[2]  = {n_rows,       w->n_cols};
    H5Sselect_hyperslab(fspace, H5S_SELECT_SET,
                        offset, NULL, count, NULL);

    /* Memory space */
    hid_t mspace = H5Screate_simple(2, count, NULL);

    /* Write */
    H5Dwrite(dset, H5T_NATIVE_FLOAT, mspace,
             fspace, H5P_DEFAULT, block);

    H5Sclose(mspace);
    H5Sclose(fspace);

    *row_counter += n_rows;
    return 0;
}

void hdf5_close(HDF5Writer *w)
{
    H5Dclose(w->dset_eu_call);
    H5Dclose(w->dset_eu_put);
    H5Dclose(w->dset_am_call);
    H5Dclose(w->dset_am_put);
    H5Fclose(w->file_id);
}