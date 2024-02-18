/* genetics.c
*/
#include "fw.h"

extern int debug, noisy;
extern int init_pop, print_opt;
extern int N, nchromo, nloci, ngen, nrepl;
extern double rec_rate0, rec_rate1, rec_rate2, rec_rate[4 * 4];
extern double prob_misseg;
extern int g_completeness;
extern double prob_selfing, prob_asexual, migration_rate;

/* pop[0] & pop[1] are curr and next generations */
extern char curr_gen, ** pop[2], * gamete[2];

/* two-locus genotype for viability selection. */
extern int ngt_sloci;
extern double gt2_freqs[4 * 4];
extern int sloci[2][2];             /* selected loci: chr0, locus0 & chr1, locus0 */
extern double gt_fitness[4 * 4];

extern double* f_allele;
extern FILE* fout;
extern double* fitness[2], * fitness_Falias;
extern int* fitness_Lalias;

extern int nthreads, thread_start, thread_step;
extern thread_data_t thread_data[];

int get_options(char* ctlf)
{
   int iopt, i, nopt = 19, lline = 4096;
   int print_GT = 0, print_h = 0, print_allelef = 0, seed = -1;
   char line[4096], * pline, opt[32], * comment = "*#/";
   char* optstr[] = { "debug", "noisy", "seed", "print", "rec_rate", "prob_misseg", "N",
      "nchromo", "nloci", "ngen", "nrepl", "selected_loci", "gt_fitness", "init_pop", "g_completeness",
      "prob_selfing", "prob_asexual", "migration_rate", "threads" };
   double t, * fit = gt_fitness;
   FILE* fctl;

   fctl = zopen(ctlf, "r");
   if (noisy) printf("Reading options from %s..\n", ctlf);
   for (; ; ) {
      if (fgets(line, lline, fctl) == NULL) break;
      for (i = 0, t = 0, pline = line; i < lline && line[i]; i++) {
         if (isalnum(line[i])) {
            t = 1; break;
         }
         else if (strchr(comment, line[i]))
            break;
      }
      if (t == 0) continue;
      sscanf(line, "%s%*s%lf", opt, &t);
      if ((pline = strstr(line, "=")) == NULL)
         zerror("option file.");

      for (iopt = 0; iopt < nopt; iopt++) {
         if (strncmp(opt, optstr[iopt], 8) == 0) {
            if (noisy >= 9)
               printf("\n%3d %15s | %-20s %6.2f", iopt + 1, optstr[iopt], opt, t);
            switch (iopt) {
            case (0): debug = (int)t;    break;
            case (1): noisy = (int)t;    break;
            case (2): seed = (int)t;     break;
            case (3):
               sscanf(pline + 1, "%d%d%d", &print_GT, &print_h, &print_allelef);
               print_opt = 1 & print_GT;
               if (print_h) print_opt += 2;
               if (print_allelef) print_opt += 4;
               break;
            case (4):
               sscanf(pline + 1, "%lf%lf%lf", &rec_rate0, &rec_rate1, &rec_rate2);
               break;
            case (5): prob_misseg = t;  break;
            case (6): N = (int)t;  break;
            case (7): nchromo = (int)t;  break;
            case (8): nloci = (int)t;  break;
            case (9): ngen = (int)t;  break;
            case (10): nrepl = (int)t;  break;
            case (11):
               sscanf(pline + 1, "%d%d%d%d", &sloci[0][0], &sloci[0][1], &sloci[1][0], &sloci[1][1]);
               if (--sloci[0][0] < 0 || sloci[0][0] >= nchromo) zerror("selected_loci outside range");
               if (--sloci[1][0] < 0 || sloci[1][0] >= nchromo) zerror("selected_loci outside range");
               if (--sloci[0][1] < 0 || sloci[0][1] >= nloci) zerror("selected_loci outside range");
               if (--sloci[1][1] < 0 || sloci[1][1] >= nloci) zerror("selected_loci outside range");
               break;
            case (12): /* gt_fitness */
               sscanf(pline + 1, "%lf%lf%lf%lf", fit + 0, fit + 1, fit + 2, fit + 3);
               for (i = 1; i < 4; i++) {
                  fscanf(fctl, "%lf%lf%lf%lf", fit + i * 4 + 0, fit + i * 4 + 1, fit + i * 4 + 2, fit + i * 4 + 3);
                  fgets(line, 1024, fctl);
               }
               break;
            case (13): init_pop = (int)t;  break;
            case (14):
               g_completeness = (int)t;  break;
            case (15):
               prob_selfing = t;  break;
            case (16):
               prob_asexual = t;  break;
            case (17):
               migration_rate = t;  break;
            case (18):
               sscanf(pline + 1, "%d%d%d", &nthreads, &thread_start, &thread_step);
               /* thread_step is not used right now. */
               break;
            }
            break;
         }
      }
      if (iopt == nopt) {
         printf("\noption %s in %s not recognised\n", opt, ctlf);
         exit(-1);
      }
   }
   fclose(fctl);

   if (nthreads < 0 || nthreads > NTHREADS || nthreads > N)
      zerror("nthreads = %d, max %d, for %d individuals..", nthreads, NTHREADS, N);
   if (--thread_start < 0) zerror("thread_start wrong..");
#if(defined(__linux__) && defined(PIN_THREADS_CORE))
   pin_to_core(thread_start);  /* this pins the master/initial thread */
#endif

   for (i = 0; i < nthreads; ++i)
      SetSeed(seed < 0 ? -1 : seed+2*i, i);

   return (0);
}

int initialize(char outf[])
{
   int i, j, ind;
   char* GT_names[4 * 4] = { "AABB", "AABb", "AAbB", "AAbb",
                             "AaBB", "AaBb", "AabB", "Aabb",
                             "aABB", "aABb", "aAbB", "aAbb",
                             "aaBB", "aaBb", "aabB", "aabb" };

   printf("2-loci selection scheme: ");
   matout2(stdout, gt_fitness, 4, 4, 10, 5);

   if (print_opt) {
      fout = zopen(outf, "w");
      fprintf(fout, "gen");

      if (print_opt & 1) /* 4x4 GT freqs at selected loci */
         for (i = 0; i < 4 * 4; i++)
            fprintf(fout, "\tf_%s", GT_names[i]);
      if (print_opt & 2) /* average heterozygosity */
         fprintf(fout, "\th");
      if (print_opt & 4) {  /* allele freq at nchromo*nloci loci */
         for (i = 0; i < nchromo; i++)
            for (j = 0; j < nloci; j++)
               fprintf(fout, "\tf_L%d_a%d", i + 1, j + 1);
      }
      fprintf(fout, "\n");
   }

   /* recombination rates, using r0 r2 */
   for (i = 0; i < 4; i++) for (j = 0; j < 4; j++)
      rec_rate[i * 4 + j] = rec_rate0;
   rec_rate[1 * 4 + 1] = rec_rate[1 * 4 + 2] = rec_rate[2 * 4 + 1] = rec_rate[2 * 4 + 2] = rec_rate2;

   /* allocate memory for current generation [0] and next generation [1], initialize as F1 */
   for (i = 0; i < 2; i++) {
      if ((pop[i] = malloc(N * sizeof(char*))) == NULL) zerror("oom");
      for (ind = 0; ind < N; ind++) {
         if ((pop[i][ind] = malloc(nchromo * nloci * sizeof(char))) == NULL)
            zerror("oom");
      }
   }

   /* gamete[2] are used to form the new zygote */
   if (nthreads > 1)
      for (i = 0; i < nthreads; i++) {
         if ((thread_data[i].gamete[0] = malloc(2 * nchromo * nloci * sizeof(char))) == NULL)
            zerror("oom");
         thread_data[i].gamete[1] = thread_data[i].gamete[0] + nchromo * nloci;
      }
   else {
      if ((gamete[0] = malloc(2 * nchromo * nloci * sizeof(char))) == NULL) zerror("oom");
      gamete[1] = gamete[0] + nchromo * nloci;
   }

   if ((f_allele = malloc(nchromo * nloci * sizeof(double))) == NULL) zerror("oom f_allele ");
   if ((fitness[0] = malloc(N * 4 * sizeof(double))) == NULL) zerror("oom fitness ");
   fitness[1] = fitness[0] + N;
   fitness_Falias = fitness[0] + N * 2;
   fitness_Lalias = (int*)(fitness[0] + N * 3);

   thread_data[0].ind_start = 0;  thread_data[0].nind = N;
   if (nthreads > 1)  threads_init();
   return(0);
}

void freemem(void)
{
   int i, ind;
   for (i = 0; i < 2; i++) {
      for (ind = 0; ind < N; ind++)
         free(pop[i][ind]);
      free(pop[i]);
   }
   if (nthreads > 1)
      for (i = 0; i < nthreads; i++)  free(thread_data[i].gamete[0]);
   else
      free(gamete[0]);
   free(f_allele);
   free(fitness[0]);
}

int meiosis(char* gamete, char* zygote, int thread_id)
{
   /* zygote -> gamete.  This returns 1 for a viable gamete and 0 for mis-segregation
   */
   int tid = thread_id, nrec, chromo, locus, i, side;  /* side = 0 or 1 */
   char* g = gamete, * z = zygote, gt0, gt;  /* gt0, gt are genotypes at adjacent loci */
   char breakpoints[1000] = { 0 };  /* crash if more than 1000 crossovers */

   /* if(missegregation), some chromosomes are not generated. */
   if (debug && thread_id == 0)
      memset(g, '.' - '0', nchromo * nloci * sizeof(char));
   for (chromo = 0; chromo < nchromo; chromo++) {
      nrec = 0;
      gt0 = *z++;
      side = (rndu(tid) < 0.5);
      *g++ = side ? (gt0 >> 1) : (gt0 & 1);
      for (locus = 1; locus < nloci; locus++) {
         gt = *z++;
         if (rndu(tid) < rec_rate[gt0 * 4 + gt]) { /* crossover with prob rec_rate */
            side = !side;
            if (debug && thread_id == 0) {
               if (nrec > 1000) zerror("many crossovers");
               breakpoints[nrec] = locus;
            }
            nrec++;
         }
         *g++ = side ? (gt >> 1) : (gt & 1);
         gt0 = gt;
      }
      if (debug && thread_id == 0) {
         printf("  chromo %d, nrec = %2d:", chromo+1, nrec);
         for (i = 0; i < nrec; i++) printf("%2d ", breakpoints[i]);
         printf("\n");
      }
      if (nrec == 0 && rndu(tid) < prob_misseg)  /* mis-segregation or gamete death */
         break;
   }
   if (debug && thread_id == 0) {
      print_zygote(stdout, zygote);
      printf("gamete: ");
      print_str(stdout, gamete);
      if (chromo != nchromo)
         printf("aha, a dead gamete..\n");
   }
   return(chromo == nchromo);  /* return 1 if a gamete is generated */
}

int individual_fitness(int curr_gen)
{
/* this calculates the fitness values for all individuals.
*/
   int i, j;
   int loc0 = sloci[0][0] * nloci + sloci[0][1];  /* 1st selected locus */
   int loc1 = sloci[1][0] * nloci + sloci[1][1];  /* 2nd selected locus */
   int gt2loci, gt0, gt1;
   char* z;
   double t, h = 0;  /* average heterozygosity at all loci */
   double* space, c[2];

   if ((space = malloc(N * sizeof(double))) == NULL) zerror("oom space");

   for (i = 0; i < N; i++) fitness[curr_gen][i] = 1;
   if (g_completeness) {
      for (i = 0; i < N; i++) {
         z = pop[curr_gen][i];
         c[0] = c[1] = 0;  /* fitness based on genome completeness */
         for (j = 0; j < nchromo * nloci; j++) {
            if (z[j] != 3) c[0] ++;
            if (z[j] != 0) c[1] ++;
         }
         fitness[curr_gen][i] = max2(c[0], c[1]) / ((double)nchromo * nloci);
      }
   }

   /* use 2-loci genotype viability to modify individual fitness */
   for (i = 0, t = 0; i < N; i++) {
      gt2loci = pop[curr_gen][i][loc0] * 4 + pop[curr_gen][i][loc1];
      t += (fitness[curr_gen][i] *= gt_fitness[gt2loci]);
   }
   for (i = 0; i < N; i++) fitness[curr_gen][i] /= t;

   MultiNomialAliasSetTable(N, fitness[curr_gen], fitness_Falias, fitness_Lalias, space);
   free(space);
   return(0);
}

void reproduction(int ind_start, int nind, int thread_id)
{
   /* pop[curr_gen] -> pop[next_gen] */
   int tid = thread_id, i, ind, igamete, parent, selfing, next_gen = !curr_gen;
   char* z, * g[2];

   if (debug && thread_id==0) {
      for (i = 0; i < N; i++) {
         printf("indiv %2d: ", i + 1);
         print_str(stdout, pop[curr_gen][i]);
      }
      printf("fitness (probs) of N = %d individuals in parent pop:", N);
      matout2(stdout, fitness[curr_gen], 1, N, 8, 4);
   }
   if (nthreads <= 1) {
      g[0] = gamete[0]; g[1] = gamete[1];
   }
   else {
      g[0] = thread_data[tid].gamete[0];
      g[1] = thread_data[tid].gamete[1];
   }

   for (ind = ind_start; ind < ind_start + nind; ind++) {
      if (debug && thread_id == 0) printf("\n**Reproducing, individual %2d**\n", ind + 1);
      z = pop[next_gen][ind];
      /* replace ind by an immigrant */
      if (migration_rate > 0 && rndu(0) < migration_rate) /* ind is immigrant */
         memset(z, 3, nchromo * nloci);
      else {

         if (prob_asexual > 0 && rndu(tid) < prob_asexual) {  /* asexual */
            parent = rndDiscreteAlias(N, fitness_Falias, fitness_Lalias);
            memmove(z, pop[curr_gen][parent], nchromo * nloci * sizeof(char));
            fitness[next_gen][ind] = fitness[curr_gen][parent];
         }
         else {
            /* if selfing, 1 gamete is produced from meiosis, the second copied. */
            selfing = (prob_selfing > 0 && rndu(tid) < prob_selfing);
            for (igamete = 0; igamete < 2 - selfing; ) {
               parent = rndDiscreteAlias(N, fitness_Falias, fitness_Lalias);
               if (meiosis(g[igamete], pop[curr_gen][parent], tid))
                  igamete++;
               if (igamete && selfing)  /* the first gamete is already produced */
                  memmove(g[1], g[0], nchromo * nloci * sizeof(char));
            }
            for (i = 0; i < nchromo * nloci; i++)
               z[i] = g[0][i] * 2 + g[1][i];
         }
      }
   }
}

void print_zygote(FILE* fout, char* zygote)
{
   int i, j, pos;
   char* g[2], gt;  /* g[2] are two gametes */

   if ((g[0] = malloc(2 * nchromo * (nloci + 1) * sizeof(char))) == NULL) zerror("oom");
   g[1] = g[0] + nchromo * (nloci + 1);

   fprintf(fout, "zygote: ");
   print_str(fout, zygote);

   for (i = 0, pos = 0; i < nchromo; i++, pos++) {
      for (j = 0; j < nloci; j++) {
         gt = zygote[i * nloci + j];
         g[0][pos] = '0' + (gt >> 1);
         g[1][pos] = '0' + (gt & 1);
         pos++;
      }
      g[0][pos] = g[1][pos] = ' ';
   }
   g[0][pos - 1] = g[1][pos - 1] = '\0';
   fprintf(fout, "        %s\n        %s\n", g[0], g[1]);
   free(g[0]);
}

void print_str(FILE* fout, char* str)
{
   /* gametes are coded 0 1 and zygotes are coded 0 1 2 3. */
   int i, j, pos = 0;

   for (i = 0, pos = 0; i < nchromo; i++) {
      for (j = 0; j < nloci; j++)
         fprintf(fout, "%c", '0' + str[i * nloci + j]);
      fprintf(fout, " ");
   }
   fprintf(fout, "\n");
}

int update_pop_features(int curr_gen)
{
   /* returns 1 if all loci are fixed */
   int nfixed = 0, i, j;
   int loc0 = sloci[0][0] * nloci + sloci[0][1];  /* 1st selected locus */
   int loc1 = sloci[1][0] * nloci + sloci[1][1];  /* 2nd selected locus */
   char gt;
   double small = 0.1/(nchromo * nloci);

   /* count 2-loci genotypes, ind i has genotype gt */
   memset(gt2_freqs, 0, ngt_sloci * sizeof(double));
   for (i = 0; i < N; i++) {
      gt = pop[curr_gen][i][loc0] * 4 + pop[curr_gen][i][loc1]; /* genotype at 2 loci */
      gt2_freqs[gt] ++;
   }
   for (i = 0; i < ngt_sloci; i++) gt2_freqs[i] /= (double)N;

   /* count alleles at each of the nchromo*nloci loci */
   memset(f_allele, 0, nchromo * nloci * sizeof(double));
   for (i = 0; i < N; i++)
      for (j = 0; j < nchromo * nloci; j++) {
         gt = pop[curr_gen][i][j];
         if (gt >> 1) f_allele[j]++;
         if (gt & 1)  f_allele[j]++;
      }
   for (j = 0; j < nchromo * nloci; j++) {
      f_allele[j] /= 2.0 * N;
      nfixed += (f_allele[j] < small || f_allele[j] > 1 - small);
   }
   return(nfixed == nchromo * nloci);
}

int print_pop_features(FILE* fout, int curr_gen)
{
   int i, j;
   double h = 0, f_alien = 0;  /* average heterozygosity at all loci*/

   if (print_opt == 0 && debug == 0) return(0);

   if (print_opt & 1) /* 4x4 GT freqs at selected loci */
      for (i = 0; i < 4 * 4; i++)
         fprintf(fout, "\t%.6f", gt2_freqs[i]);
   if (print_opt & 2) { /* average heterozygosity */
      for (i = 0; i < N; i++) for (j = 0; j < nchromo * nloci; j++)
         if (pop[curr_gen][i][j] == 1 || pop[curr_gen][i][j] == 2) h++;
      for (j = 0; j < nchromo * nloci; j++)
         f_alien += f_allele[j] / (nchromo * nloci);
      fprintf(fout, "\t%.6f\t%.6f", h / (N * nchromo * nloci), f_alien);
   }
   if (print_opt & 4) {  /* allele freq at nchromo*nloci loci */
      for (i = 0; i < nchromo * nloci; i++)
         fprintf(fout, "\t%.6f", f_allele[i]);
   }
   if (print_opt) fprintf(fout, "\n");

   if (debug) {
      printf("\nfreqs for 4x4 genotypes at 2 selected loci\n");
      for (i = 0; i < 4; i++) {  /* */
         for (j = 0; j < 4; j++)
            printf("%10.6f", gt2_freqs[i * 4 + j]);
         printf("\n");
      }
   }
   return(0);
}
