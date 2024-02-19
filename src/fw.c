/* fw.c
   Fisher-Wright simulation program for forward simualtion of reproduction, 
   recombination, and introgression.
   ziheng yang, 29 March 2022

   gcc -o fw -Wno-unused-result -O3 fw.c genetics.c threads.c tools-sub.c -lm -lpthread
   icc -o fw -O3 fw.c genetics.c threads.c tools-sub.c -lm -lpthread
   cl -Ox fw.c genetics.c threads.c tools-sub.c -link libpthreadVC3.lib
   ./fw fw-control0.txt
*/
#include "fw.h"

int debug = 1, noisy = 3;
int init_pop = 1, print_opt = 1;
int N = 1000, nchromo = 16, nloci = 200, ngen = 1000, nrepl=1000;
double rec_rate0 = 0.057, rec_rate1 = 0.057, rec_rate2 = 0.057 / 100, rec_rate[4 * 4];
double prob_misseg = 0.45;
int g_completeness;
double prob_selfing = 0.15, prob_asexual = 0.05, migration_rate = 1e-4;

/* pop[0] & pop[1] are curr and next generations */
char curr_gen = 0, ** pop[2], * gamete[2];

/* two-locus genotype for viability selection. */
int ngt_sloci = 4 * 4;
double gt2_freqs[4 * 4];
int sloci[2][2] = { {0,0}, {1,0} };            /* selected loci: chr0, locus0 & chr1, locus0 */
double gt_fitness[4 * 4] = { 1,  1,  1,  1,    /* AABB AABb AAbB AAbb */
                             1,  1,  1,  1,    /* AaBB AaBb AabB Aabb */
                             1,  1,  1,  1,    /* aABB aABb aAbB aAbb */
                             1,  1,  1,  1};   /* aaBB aaBb aabB aabb */

double* f_allele;
FILE* fout;
double* fitness[2], * fitness_Falias;
int* fitness_Lalias;

int nthreads = 1, thread_start = 0, thread_step = 1;
extern thread_data_t thread_data[];

double fixation_probability(FILE* fout, int nr)
{
   int ir;
   double h = 0, f_alien = 0;  /* average heterozygosity at all loci*/

   if (print_opt == 0 && debug == 0) return(0);
   update_pop_features(curr_gen);

   for (ir = 0; ir < nr; ir++) {
      ;
   }
   return(0);
}

int initial_population(void)
{
   int i, j, ind;

   curr_gen = 0;

   if (migration_rate > 0) {  /* make 1 ind an immigrant (or genotype 1/1) */
      ind = (int)(N * rndu(0));
      memset(pop[curr_gen][ind], 3, nchromo * nloci); /* 1 immigrant of pure type 1/1 */
   }
   else if (init_pop == 1)     /* F1 heterozygote 0/1 */
      for (ind = 0; ind < N; ind++)
         memset(pop[curr_gen][ind], 1, nchromo * nloci);
   else                        /* HW with p=0.5 */
      for (ind = 0; ind < N; ind++)
         for (j = 0; j < nchromo * nloci; j++)
            pop[curr_gen][ind][j] = (char)(rndu(0) * 4);
   return(0);
}

int termination(int option)
{
/* returns 1 if population is fixed and simulation should be terminated.
*/
   int i;
   int loc0 = sloci[0][0] * nloci + sloci[0][1];  /* 1st selected locus */
   int loc1 = sloci[1][0] * nloci + sloci[1][1];  /* 2nd selected locus */
   double small = 0.1 / N;

   switch (option) {
   case 0:   /* stop if every locus is fixed. */
      for (i = 0; i < nchromo * nloci; i++) {
         if (f_allele[i] > small && 1 - f_allele[i] < small)
            return(0);
      }
      return (-1);
   case 1:  /* stop if both selected loci are fixed. */
      if ((f_allele[loc0] < small || 1 - f_allele[loc0] < small) &&
          (f_allele[loc1] < small || 1 - f_allele[loc1] < small))
         return 1;
   }
   return(0);
}

int simulate_1replicate(int model, int replicate)
{
   int i, igen;
   char timestr[64];

   for (igen = 0; igen < ngen; igen++) {
      printf("\r*** Replicate %4d gen %4d ..", replicate+1, igen);
      fflush(stdout);

      //if (debug || (noisy > 3 && N >= 100 && ngen > 100 && (igen + 1) % 5 == 0))
      //   printf("\n*** Gen %4d (%s): ***", igen + 1, printtime(timestr));
      if (print_opt)  fprintf(fout, "%d", igen);
      individual_fitness(curr_gen);
      update_pop_features(curr_gen);
      if (print_opt)
         print_pop_features(fout, curr_gen);
      if (termination(model)) 
         break;
      if (nthreads > 1) {
         /* work: -1: end; 0: idle waiting for work; 1: reproduction */
         threads_wakeup((int)1, NULL);
      }
      else
         reproduction(0, N, 0);
      curr_gen = !curr_gen;
   }

   if (print_opt) {
      fprintf(fout, "\ngenotypes in generation #%d (N = %d)\n", ngen, N);
      for (i = 0; i < N; i++) {
         fprintf(fout, "indiv %2d: ", i + 1);
         print_str(fout, pop[!curr_gen][i]);
      }
   }
   return(igen);
}

int main(int argc, char* argv[])
{
   char ctlf[4096] = "fw-control0.txt", outf[4096] = "out.txt", resf[4096] = "result.txt";
   char timestr[64];
   int model = 1, ir, i, AAbb;
   double meantime = 0, t, pAAbb = 0, small = 0.1/N;
   FILE* fresult;

   starttimer();
   if (argc > 1) strcpy(ctlf, argv[1]);
   if (argc > 2) strcpy(outf, argv[2]);
   if (argc > 3) strcpy(resf, argv[3]);
   get_options(ctlf);
   initialize(outf);

   fresult = zopen(resf, "w");
   fprintf(fresult, "replicate\ttime\n");

   for (ir = 0; ir < nrepl; ir++) {
      //printf("\n*** Replicate %4d", ir + 1);
      initial_population();
      t = simulate_1replicate(model, ir);
      AAbb = (1 - gt2_freqs[0 * 4 + 3] < small);
      if (1-gt2_freqs[0 * 4 + 3] < small)
         meantime += (t - meantime) / ++pAAbb;  /* running mean */
      printf("\r*** Replicate %4d (%s): %6.0f generations, %s\n", ir + 1, printtime(timestr), t, (AAbb ? "AAbb" : "  "));
      fprintf(fresult, "%d\t%.0f\n", ir+1, t);
   }

   printf("\nP_AAbb = %7.4f", pAAbb/nrepl);
   printf("\nmean waitign time = %7.1f\n", meantime);

   freemem();
   fclose(fresult);
   if (print_opt) fclose(fout);
   if (nthreads > 1) threads_exit();
}
