#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "norme.h"

double norme(int n, int *ia, int *ja, double *a, double *b, double *x)
/*
  Goal
  ===
  Calculation of the residual norm of the system Ax = b 
 
                  ||b - Ax|| 
   norme_residu = ==========
                    ||b||
            
  Arguments
  =========                
   (input)  n  - pointer of the number of unknowns
   (input)  ia - pointer of the vector 'ia' of the matrix A (fine grid)
   (input)  ja - pointer of the vector 'ja' of the matrix A (fine grid)
   (input)  a  - pointer of the vector 'a' of the matrix A (fine grid)
   (input)  b  - pointer of the vector 'b' of the matrix b (fine grid)
   (input)  x  - pointer of the solution vector  
   
  (ouput) norme_residu - residual norm of the system Ax = b               
*/ 
{
	int i,j;
	double *y = malloc(n * sizeof(double)); /* y = A * x */
	double *r = malloc(n * sizeof(double)); /* r = y - b */
	
 	double norme_carre_r = 0.0; /*Initialisation du carré de la norme de r */  
	double norme_carre_b = 0.0; /*Initialisation du carré de la norme de b */
	double norme_residu; /*Initialisation de la norme du residu */ 
	
	
	for (i=0; i<n; i++){
		y[i] = 0.0; /*Initialisation de chaque élément de y à 0 */
	}
	
	/* Produit de la matrice A au format CSR par le vecteur solution x */
	
	for (i=0; i<n; i++){
		for (j=ia[i]; j< (ia[i+1]); j++){
			y[i] += a[j] * x[ja[j]];   
		}
		
		r[i] = y[i] - b[i]; 
		norme_carre_r +=  r[i] * r[i]; 
		norme_carre_b +=  b[i] * b[i];
	}
	
	/* Calcul de la norme du résidu */ 
	norme_residu = sqrt(norme_carre_r) / sqrt(norme_carre_b);
	
	
	/* Libération de la mémoire */
	free(y); free(r);
	
	return norme_residu;
}