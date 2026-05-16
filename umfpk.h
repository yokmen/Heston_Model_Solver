/* Prototypes */
int  factor_umfpack(int n, int *ia, int *ja, double *a, void **Numeric);
int  solve_umfpack (int n, int *ia, int *ja, double *a, double *b, double *x, void *Numeric);
void free_umfpack  (void **Numeric);
