#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "time.h"
#include "dataset.h"

int main(void)
{
   /* ------------------------------------------------------------------
      Parametres de la grille
   ------------------------------------------------------------------ */

   // Conditions pour les paramètres de grille : 
   // 1) (N_S-1) divisible par 3
   // 2) (N_v-1) divisible par 8
   // 3) N_t >= N_S /10
   int    N_S = 100;
   int    N_v = 81;
   int    N_t = 100;                   
   double S_min = 0.0;
   double v_min = 0.0;
   double K     = 1.0;
   double S_max = 3.0 * K;
   
   /* ------------------------------------------------------------------
      Mesure du temps
   ------------------------------------------------------------------ */
   double tc1 = mytimer_cpu();
   double tw1 = mytimer_wall();
   
   time_t start_time = time(NULL);
   char time_str[100];
   strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S", localtime(&start_time));
   
   printf("\n========================================\n");
   printf("  HESTON MODEL - DATASET GENERATION\n");
   printf("========================================\n");
   printf("Start time: %s\n", time_str);
   printf("\nGrid parameters:\n");
   printf("  N_S = %d, N_v = %d, N_t = %d\n", N_S, N_v, N_t);
   printf("  S range: [%.2f, %.2f], K = %.2f\n", S_min, S_max, K);
   printf("  alpha_S = %.2f, alpha_v_rel = %.2f\n", 0.2, 0.1);
   printf("\nGenerating dataset...\n\n");
   
   //int n, *ia, *ja;
   //double *a;

   /* Crank-Nicolson */

   /* ------------------------------------------------------------------
      Type d'option
      CALL 0 / PUT 1
      EUROPEAN 0 / AMERICAN 1
   ------------------------------------------------------------------ */

   double alpha_S = 0.2; // Tavella-Randall grid density parameter for S
   double alpha_v_rel = 0.1; // Relative grid density parameter for v (
   
   dataset( N_S,  N_v,  N_t, S_min,  S_max, v_min, K, alpha_S, alpha_v_rel);

   /* ------------------------------------------------------------------
      Mesure du temps
   ------------------------------------------------------------------ */
   double tc2 = mytimer_cpu();
   double tw2 = mytimer_wall();
   
   time_t end_time = time(NULL);
   strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S", localtime(&end_time));
   
   printf("\n========================================\n");
   printf("  GENERATION COMPLETE\n");
   printf("========================================\n");
   printf("End time  : %s\n", time_str);
   printf("\nTiming Summary:\n");
   printf("  Wall-clock time  : %.2f seconds\n", tw2 - tw1);
   printf("  CPU time         : %.2f seconds\n", tc2 - tc1);
   printf("  Parallelism ratio: %.2f x\n", (tc2 - tc1) / (tw2 - tw1));
   printf("\nOutput file: dataset_Options.h5\n");
   printf("========================================\n\n");
   
   return 0;
}