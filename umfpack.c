#include <stdio.h>
#include "umfpack.h"
#include "umfpk.h"

// Factorisation : appelée UNE FOIS
int factor_umfpack(int n, int *ia, int *ja, double *a, void **Numeric)
{
    void *Symbolic;
    double Info[UMFPACK_INFO], Control[UMFPACK_CONTROL];

    int status = umfpack_di_symbolic(n, n, ia, ja, a, &Symbolic, Control, Info);
    if (status < 0) {
        printf("ERREUR : umfpack_di_symbolic a échoué\n");
        return 1;
    }

    status = umfpack_di_numeric(ia, ja, a, Symbolic, Numeric, Control, Info);
    if (status < 0) {
        printf("ERREUR : umfpack_di_numeric a échoué\n");
        return 1;
    }

    umfpack_di_free_symbolic(&Symbolic);
    return 0;
}

// Solve : appelé à CHAQUE PAS DE TEMPS
int solve_umfpack(int n, int *ia, int *ja, double *a,
                  double *b, double *x, void *Numeric)
{
    double Info[UMFPACK_INFO], Control[UMFPACK_CONTROL];

    int status = umfpack_di_solve(UMFPACK_At, ia, ja, a, x, b,
                                  Numeric, Control, Info);
    if (status < 0) {
        printf("ERREUR : umfpack_di_solve a échoué\n");
        return 1;
    }
    return 0;
}

// Libération : appelée UNE FOIS après la boucle
void free_umfpack(void **Numeric)
{
    umfpack_di_free_numeric(Numeric);
}