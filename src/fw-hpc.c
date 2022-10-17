﻿/* Fisher-Wright simulation program for forward simualtion 
   of reproduction, recombination, and introgression.
   ziheng yang, 29 March 2022

   gcc -o fw-hpc -Wno-unused-result -O3 fw-hpc.c fw-sub.c -lm -lpthread
   icc -o fw-hpc -O3 fw-hpc.c fw-sub.c -lm -lpthread
   cl -Ox fw-hpc.c fw-sub.c -link libpthreadVC3.lib
   ./fw-hpc fw-control0.txt
*/

#if(defined(__linux__))
#define _GNU_SOURCE
#include<sched.h>
#include<unistd.h>
//#include <sys/resource.h>
//#include <sys/sysinfo.h>
#endif

#include<pthread.h>
#include "fw.h"
int meiosis(char* gamete, char* zygote, int thread_id);
void reproduction(int ind_start, int nind, int thread_id);

int debug = 1, noisy = 3;
int N = 1000, nchromo = 16, nloci = 200, ngen = 1000;
double rec_rate0 = 0.057, rec_rate1 = 0.057, rec_rate2 = 0.057 / 100, rec_rate[4 * 4];
double prob_misseg = 0.45;
int init_pop=1, print_opt = 1;

/* pop[0] & pop[1] are curr and next generations */
char curr_pop = 0, ** pop[2], * gamete[2];

/* two-locus genotype for viability selection.
   gt_list[16] lists the individuals for each of the 16 GTs */
int ngt = 4 * 4, gt_counts[4 * 4], * gt_list[4 * 4], gt_Lalias[4 * 4];
double gt_probs[4 * 4], gt_freqs[4 * 4], gt_Falias[4 * 4];
int sloci[2][2] = { {0,0}, {1,0} };             /* selected loci: chr0, locus0 & chr1, locus0 */
double gt_fitness[4 * 4] = { 1,  1,  1,  1,    /* AABB AABb AAbB AAbb */
                            50, 50, 50, 50,    /* AaBB AaBb AabB Aabb */
                            50, 50, 50, 50,    /* aABB aABb aAbB aAbb */
                            100,50, 50, 1 };   /* aaBB aaBb aabB aabb */
double* f_allele;
FILE* fout;
double gt_Falias[4 * 4];
int gt_Lalias[4 * 4];

int nthreads = 4, thread_start = 0;

struct thread_data_s {
   pthread_t thread;
   pthread_mutex_t mutex;
   pthread_cond_t cond;

   /* work: -1: end; 0: idle waiting for work; 1: meiosis */
   volatile int work;
   int id, ind_start, nind;
   char *gamete[2];
}  thread_data[NTHREADS];


#if(defined(__linux__) && defined(PIN_THREADS_CORE))
void pin_to_core(int t)
{
   cpu_set_t cpuset;
   CPU_ZERO(&cpuset);
   CPU_SET(t, &cpuset);
   if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset)) {
      printf("Error while pinning thread to core. ");
      exit(-1);
   }
}
#endif

void* thread_worker(void* arg)
{
   int tid = (int)arg;
   int i0 = thread_data[tid].ind_start, i1 = i0 + thread_data[tid].nind;

#if(defined(__linux__) && defined(PIN_THREADS_CORE))
   pin_to_core(thread_start + tid);
#endif

   pthread_mutex_lock(&thread_data[tid].mutex);
   /* loop until signalled to quit */
   while (thread_data[tid].work >= 0) {
      /* wait for work available */
      if (thread_data[tid].work == 0)
         pthread_cond_wait(&thread_data[tid].cond, &thread_data[tid].mutex);
      if (thread_data[tid].work > 0) {
         if (debug) printf("Thread %d working, ind %5d -%5d\n", tid, i0 + 1, i1);

         /* work: -1: end; 0: idle waiting for work; 1: baby boom */
         if (thread_data[tid].work == 1)
            reproduction(thread_data[tid].ind_start, thread_data[tid].nind, tid);
         else
            zerror("thread work not recognised...");

         if (debug) printf("\t\t\t\tThread %d finished, loci %5d -%5d\n", tid, i0 + 1, i1);
         thread_data[tid].work = 0;
         pthread_cond_signal(&thread_data[tid].cond);
      }
   }
   pthread_mutex_unlock(&thread_data[tid].mutex);
   pthread_exit(NULL);

   return(NULL);
}

void threads_init(void)
{
   int tid;
   for (tid = 0; tid < nthreads; ++tid) {
      thread_data[tid].work = 0;
      pthread_mutex_init(&thread_data[tid].mutex, NULL);
      pthread_cond_init(&thread_data[tid].cond, NULL);
      if (pthread_create(&thread_data[tid].thread, NULL, thread_worker, (void*)tid))
         puts("Cannot create thread");
   }
}

void threads_wakeup(int work, void* data)
{
   int tid;
   /* put threads to work */
   for (tid = 0; tid < nthreads; ++tid) {
      pthread_mutex_lock(&thread_data[tid].mutex);
      thread_data[tid].work = work;
      pthread_cond_signal(&thread_data[tid].cond);
      pthread_mutex_unlock(&thread_data[tid].mutex);
   }
   /* wait for threads to finish work */
   for (tid = 0; tid < nthreads; ++tid) {
      pthread_mutex_lock(&thread_data[tid].mutex);
      while (thread_data[tid].work > 0)
         pthread_cond_wait(&thread_data[tid].cond, &thread_data[tid].mutex);
      pthread_mutex_unlock(&thread_data[tid].mutex);
   }
   /* TODO: threads have now finished, master thread collects results */
}

void threads_exit()
{
   for (int tid = 0; tid < nthreads; ++tid) {
      pthread_mutex_lock(&thread_data[tid].mutex);
      thread_data[tid].work = -1;
      pthread_cond_signal(&thread_data[tid].cond);
      pthread_mutex_unlock(&thread_data[tid].mutex);
      /* wait for worker to quit */
      if (pthread_join(thread_data[tid].thread, 0))
         puts("Cannot join thread");
      pthread_cond_destroy(&thread_data[tid].cond);
      pthread_mutex_destroy(&thread_data[tid].mutex);
   }
}


int get_options(char* ctlf)
{
   int iopt, i, nopt = 13, lline = 4096;
   int print_GT = 0, print_h = 0, print_allelef = 0;
   char line[4096], * pline, opt[32], * comment = "*#/";
   char* optstr[] = { "debug", "noisy", "rec_rate", "prob_misseg", "N", "nchromo", "nloci", "ngen",
      "selected_loci", "gt_fitness", "print", "init_pop", "threads" };
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
            case (2):
               sscanf(pline + 1, "%lf%lf%lf", &rec_rate0, &rec_rate1, &rec_rate2);
               break;
            case (3): prob_misseg = t;  break;
            case (4): N = (int)t;  break;
            case (5): nchromo = (int)t;  break;
            case (6): nloci = (int)t;  break;
            case (7): ngen = (int)t;  break;
            case (8):
               sscanf(pline + 1, "%d%d%d%d", &sloci[0][0], &sloci[0][1], &sloci[1][0], &sloci[1][1]);
               if (--sloci[0][0] < 0 || sloci[0][0] >= nchromo) zerror("selected_loci outside range");
               if (--sloci[1][0] < 0 || sloci[1][0] >= nchromo) zerror("selected_loci outside range");
               if (--sloci[0][1] < 0 || sloci[0][1] >= nloci) zerror("selected_loci outside range");
               if (--sloci[1][1] < 0 || sloci[1][1] >= nloci) zerror("selected_loci outside range");
               break;
            case (9): /* gt_fitness */
               sscanf(pline + 1, "%lf%lf%lf%lf", fit + 0, fit + 1, fit + 2, fit + 3);
               for (i = 1; i < 4; i++) {
                  fscanf(fctl, "%lf%lf%lf%lf", fit + i * 4 + 0, fit + i * 4 + 1, fit + i * 4 + 2, fit + i * 4 + 3);
                  fgets(line, 1024, fctl);
               }
               break;
            case (10):
               sscanf(pline + 1, "%d%d%d", &print_GT, &print_h, &print_allelef);
               print_opt = 1 & print_GT;
               if (print_h) print_opt += 2;
               if (print_allelef) print_opt += 4;
               break;
            case (11): init_pop = (int)t;  break;
            case (12):
               sscanf(pline + 1, "%d%d", &nthreads, &thread_start);
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
   if (nthreads <0 || nthreads > NTHREADS || nthreads > N)
      zerror("nthreads = %d, max %d, for %d individuals..", nthreads, NTHREADS, N);
   if (--thread_start < 0) zerror("thread_start wrong..");
   return (0);
}

int initialize(char outf[])
{
   int i, j, ind;
   char* GT_names[4 * 4] = { "AABB", "AABb", "AAbB", "AAbb",
                             "AaBB", "AaBb", "AabB", "Aabb",
                             "aABB", "aABb", "aAbB", "aAbb",
                             "aaBB", "aaBb", "aabB", "aabb" };

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
               fprintf(fout, "\tf_L%d_a%d", i+1, j+1);
      }
      fprintf(fout, "\n");
   }

   /* recombination rates, using r0 r2 */
   for (i = 0; i < 4; i++) for (j = 0; j < 4; j++)
      rec_rate[i * 4 + j] = rec_rate0;
   rec_rate[1 * 4 + 1] = rec_rate[1 * 4 + 2] = rec_rate[2 * 4 + 1] = rec_rate[2 * 4 + 2] = rec_rate2;

   if ((gt_list[0] = malloc(ngt * N * sizeof(int))) == NULL)
      zerror("oom");
   memset(gt_list[0], -1, ngt * N * sizeof(int));
   for (i = 1; i < ngt; i++)
      gt_list[i] = gt_list[i - 1] + N;

   /* allocate memory for current generation [0] and next generation [1], initialize as F1 */
   for (i = 0; i < 2; i++) {
      if ((pop[i] = malloc(N * sizeof(char*))) == NULL) zerror("oom");
      for (ind = 0; ind < N; ind++)
         if ((pop[i][ind] = malloc(nchromo * nloci * sizeof(char))) == NULL) zerror("oom");
   }
   for (ind = 0; ind < N; ind++) {
      for (j = 0; j < nchromo * nloci; j++)
         if (init_pop == 1)
            pop[curr_pop][ind][j] = (char)1;              /* F1 heterozygote 0/1 */
         else
            pop[curr_pop][ind][j] = (char)(rndu(0) * 4);  /* HW with p=0.5 */
   }

   /* gamete[2] are used to form the new zygote */
   if(nthreads>1) 
      for (i = 0; i < nthreads; i++) {
         if ((thread_data[i].gamete[0] = malloc(2 * nchromo * nloci * sizeof(char))) == NULL)
            zerror("oom");
         thread_data[i].gamete[1] = thread_data[i].gamete[0] + nchromo * nloci;
      }
   else {   /* gamete[2] are used to form the new zygote */
      if ((gamete[0] = malloc(2 * nchromo * nloci * sizeof(char))) == NULL) zerror("oom");
      gamete[1] = gamete[0] + nchromo * nloci;
   }

   if ((f_allele = malloc(nchromo * nloci * sizeof(double))) == NULL) zerror("oom");

   return(0);
}

void freemem()
{
   int i, ind;
   free(gt_list[0]);
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
}

int meiosis(char* gamete, char* zygote, int thread_id)
{
   /* zygote -> gamete.  This returns 1 for a viable gamete and 0 for mis-segregation
   */
   int tid = thread_id, nrec, chromo, locus, i, side;  /* side = 0 or 1 */
   char* g = gamete, * z = zygote, gt0, gt;  /* gt0 gt are genotypes at adjacent loci */
   char breakpoints[1000] = { 0 };  /* crash if more than 1000 crossovers */

   /* if(missegregation), some chromosomes are not generated. */
   if (debug) memset(g, 'x' - '0', nchromo * nloci * sizeof(char));
   for (chromo = 0; chromo < nchromo; chromo++) {
      nrec = 0;
      gt0 = *z++;
      side = (rndu(tid) < 0.5);
      *g++ = side ? (gt0 >> 1) : (gt0 & 1);
      for (locus = 1; locus < nloci; locus++) {
         gt = *z++;
         if (rndu(tid) < rec_rate[gt0 * 4 + gt]) { /* crossover with prob rec_rate */
            side = !side;
            if (debug) {
               if (nrec > 1000) zerror("many crossovers");
               breakpoints[nrec] = locus;
            }
            nrec++;
         }
         *g++ = side ? (gt >> 1) : (gt & 1);
         gt0 = gt;
      }
      if (debug) {
         printf("chromo %d, nrec = %2d:", chromo, nrec);
         for (i = 0; i < nrec; i++) printf("%2d ", breakpoints[i]);
         printf("\n");
      }
      if (nrec == 0 && rndu(tid) < prob_misseg)  /* mis-segregation or gamete loss */
         break;
   }
   if (debug) {
      print_zygote(stdout, zygote);
      print_gamete(stdout, gamete);
   }
   return(chromo == nchromo);  /* return 1 if a gamete is generated */
}

void reproduction(int ind_start, int nind, int thread_id)
{
   int tid=thread_id, i, ind, ngamete, parent;
   char gt, *z, *g[2];

   g[0] = gamete[0]; g[1] = gamete[1];
   if (nthreads > 1) {
      g[0] = thread_data[tid].gamete[0]; g[1] = thread_data[tid].gamete[1];
   }

   for (ind = ind_start; ind < ind_start + nind; ind++) {
      for (ngamete = 0; ngamete < 2; ) {
         gt = rndDiscreteAlias(ngt, gt_Falias, gt_Lalias);
         parent = (int)(gt_counts[gt] * rndu(tid));
         if (meiosis(g[ngamete], pop[curr_pop][parent], tid))  ngamete++;
      }
      z = pop[1 - curr_pop][ind];
      for (i = 0; i < nchromo * nloci; i++)
         z[i] = (g[0][i] << 1) + g[1][i];
   }
}


void print_zygote(FILE* fout, char* zygote)
{
   int i, j, pos;
   char* g[2], * z, gt;  /* g[2] are two gametes */

   if ((g[0] = malloc(2 * nchromo * (nloci + 1) * sizeof(char))) == NULL) zerror("oom");
   g[1] = g[0] + nchromo * (nloci + 1);
   for (i = 0, z = g[0], pos = 0; i < nchromo; i++) {
      for (j = 0; j < nloci; j++)
         z[pos++] = '0' + zygote[i * nloci + j];
      z[pos++] = ' ';
   }
   z[pos - 1] = '\0';
   fprintf(fout, "zygote: %s\n", z);

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

void print_gamete(FILE* fout, char* gamete)
{
   int i, j, pos = 0;
   char* g;

   if ((g = malloc(nchromo * (nloci + 1) * sizeof(char))) == NULL) zerror("oom");
   for (i = 0, pos = 0; i < nchromo; i++) {
      for (j = 0; j < nloci; j++)
         g[pos++] = '0' + gamete[i * nloci + j];
      g[pos++] = ' ';
   }
   g[pos - 1] = '\0';
   fprintf(fout, "gamete: %s\n", g);
   free(g);
}

int update_pop_features(char** popu)
{
   /* returns 1 if all loci are fixed */
   int nfixed, i, j, loc0, loc1;
   char gt;
   double space[4*4], t, h = 0;  /* average heterozygosity at all loci*/

   /* count 2-loci genotypes, ind i has genotype gt */
   memset(gt_counts, 0, ngt * sizeof(int));
   memset(gt_freqs, 0, ngt * sizeof(double));
   for (i = 0; i < N; i++) {
      loc0 = sloci[0][0] * 4 + sloci[0][1];
      loc1 = sloci[1][0] * 4 + sloci[1][1];
      gt = popu[i][loc0] * 4 + popu[i][loc1];
      gt_list[gt][gt_counts[gt]++] = i;
   }
   for (i = 0; i < ngt; i++) gt_freqs[i] = gt_counts[i] / (double)N;
   t = 0;
   for (i = 0; i < ngt; i++) t += gt_probs[i] = gt_freqs[i] * gt_fitness[i];
   for (i = 0; i < ngt; i++) gt_probs[i] /= t;

   MultiNomialAliasSetTable(ngt, gt_probs, gt_Falias, gt_Lalias, space);
   //MultiNomialAlias(N, ngt, gt_Falias, gt_Lalias, gt_counts);

   /* count alleles at each of the nchromo*nloci loci */
   memset(f_allele, 0, nchromo * nloci * sizeof(double));
   for (i = 0; i < N; i++) for (j = 0; j < nchromo * nloci; j++) {
      gt = popu[i][j];
      if (gt >> 1) f_allele[j]++;
      if (gt & 1)  f_allele[j]++;
   }
   for (j = 0, nfixed = 0; j < nchromo * nloci; j++) {
      f_allele[j] /= 2.0 * N;
      nfixed += (f_allele[j] < 0.5 / (nchromo * nloci) || f_allele[j]>1 - 0.5 / (nchromo * nloci));
   }
   return(nfixed == nchromo * nloci);
}

int print_pop_features(FILE* fout, char** popu, int update)
{
   int i, j;
   double h = 0;  /* average heterozygosity at all loci*/

   if (print_opt == 0 && debug == 0) return(0);
   if (update)  update_pop_features(popu);

   if (print_opt & 1) /* 4x4 GT freqs at selected loci */
      for (i = 0; i < 4 * 4; i++)
         fprintf(fout, "\t%.6f", gt_freqs[i]);
   if (print_opt & 2) { /* average heterozygosity */
      for (i = 0; i < N; i++) for (j = 0; j < nchromo * nloci; j++)
         if (popu[i][j] == 1 || popu[i][j] == 2) h++;
      fprintf(fout, "\t%.6f", h / (N * nchromo * nloci));
   }
   if (print_opt & 4) {  /* allele freq at nchromo*nloci loci */
      for (i = 0; i < nchromo * nloci; i++)
         fprintf(fout, "\t%.6f", f_allele[i]);
   }
   if (print_opt) fprintf(fout, "\n");

   if (debug) {
      printf("\n4x4 freqs & sampling-probs for 4x4 genotypes at selected loci\n");
      for (i = 0; i < 4; i++) {  /* */
         for (j = 0; j < 4; j++)
            printf("%10.6f", gt_freqs[i*4+j]);
         printf("   ");
         for (j = 0; j < 4; j++)
            printf("%10.6f", gt_probs[i*4+j]);
         printf("\n");
      }
      matout2(stdout, gt_probs, 4, 4, 10, 6);
   }

   return(0);
}

int main(int argc, char* argv[])
{
   char ctlf[4096] = "fw-control0.txt", outf[4096] = "out.txt", timestr[64];
   int i, igen, status;

   if (argc > 1) strcpy(ctlf, argv[1]);
   if (argc > 2) strcpy(outf, argv[2]);
   get_options(ctlf);
   starttimer();

   for (i = 0; i < nthreads; ++i)
      SetSeed(-123, i);
   printf("2-loci selection scheme: ");
   matout2(stdout, gt_fitness, 4, 4, 10, 5);

   initialize(outf);
   // thread_data[0].ind_start = 0; thread_data[0].nind = N;
   if (nthreads > 1) {
      int start = 0, remaining = N% nthreads, per_thread = N / nthreads;
      for (i = 0; i < nthreads; ++i) {
         thread_data[i].ind_start = start;
         thread_data[i].nind = per_thread + (remaining > 0 ? 1 : 0);
         if (remaining) remaining--;
         start += thread_data[i].nind;
         printf("Thread %2d : ind %5d -- %5d\n", i + 1, thread_data[i].ind_start + 1, thread_data[i].ind_start + thread_data[i].nind);
      }
      threads_init();
   }

   for (igen = 0; igen < ngen; igen++) {
      if (noisy >= 3 || (N >= 100 && ngen > 100 && (igen+1) % 10 == 0))
         printf("\n* Gen %4d (%s): *", igen+1, printtime(timestr));
      if (print_opt)
         fprintf(fout, "%d", igen);
      status = update_pop_features(pop[curr_pop]);
      if(print_opt) print_pop_features(fout, pop[curr_pop], 0);
      if (status == 1) break;  /* all loci are fixed, no point of continuing */

      if (nthreads > 1) {
         /* work: -1: end; 0: idle waiting for work; 1: update_times; 2: loglike calculation */
         threads_wakeup((int)1, NULL);
      }
      else
         reproduction(0, N, 0);

      curr_pop = 1 - curr_pop;
   }
   freemem();
   if (print_opt) fclose(fout);
   if (nthreads > 1) threads_exit();
}