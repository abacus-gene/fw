/* fw-sub.c
*/

#include "fw.h"

unsigned int z_rndu[NTHREADS] = { 1237,1237,1237,1237,1237,1237,1237 };

double rndu(int thread_id)
{
   /* From Ripley (1987) p. 46 or table 2.4 line 2.  32-bit integer assumed.  */
   z_rndu[thread_id] = z_rndu[thread_id] * 69069 + 1;
   if (z_rndu[thread_id] == 0)  z_rndu[thread_id] = 12345671;
   return ldexp((double)(z_rndu[thread_id]), -32);
}

void SetSeed(int seed, int thread_id)
{
   /* Note seed is of type int with -1 meaning "please find a seed". */
   if (seed <= 0)
      seed = abs(2 * (int)time(NULL) + 1);
   z_rndu[thread_id] = (unsigned int)seed;
}

int zero(double x[], int n)
{
   int i;
   for (i = 0; i < n; i++)
      x[i] = 0;
   return (0);
}

double sum(double x[], int n)
{
   int i;
   double t = 0;
   for (i = 0; i < n; i++)
      t += x[i];
   return (t);
}

double rndNormal(void)
{
   /* Standard normal variate, using the Box-Muller method (1958), improved by
      Marsaglia and Bray (1964).  The method generates a pair of N(0,1) variates,
      but only one is used.
      Johnson et al. (1994), Continuous univariate distributions, vol 1. p.153.
   */
   double u, v, s;

   for ( ; ; ) {
      u = 2 * rndu(0) - 1;
      v = 2 * rndu(0) - 1;
      s = u * u + v * v;
      if (s > 0 && s < 1)
         break;
   }
   s = sqrt(-2 * log(s) / s);
   return (u * s); /* (v*s) is the other N(0,1) variate, wasted. */
}

void zerror(const char* format, ...)
{
   va_list argptr;
   va_start(argptr, format);
   fprintf(stderr, "\nerror: ");
   vfprintf(stderr, format, argptr);
   va_end(argptr);
   exit(1);
}

static time_t time_start;
void starttimer(void)
{
   time_start = time(NULL);
}

char* printtime(char timestr[])
{
   /* print time elapsed since last call to starttimer()
    */
   time_t t;
   int h, m, s;

   t = time(NULL) - time_start;
   h = (int)t / 3600;
   m = (int)(t % 3600) / 60;
   s = (int)(t - (t / 60) * 60);
   if (h)
      sprintf(timestr, "%d:%02d:%02d", h, m, s);
   else
      sprintf(timestr, "%2d:%02d", m, s);
   return (timestr);
}


FILE* zopen(char* filename, char* mode)
{
   FILE* fp;

   if (filename == NULL || filename[0] == 0)
      zerror("file name empty.");

   fp = (FILE*)fopen(filename, mode);
   if (fp == NULL) {
      printf("\nerror when opening file %s\n", filename);
      if (!strchr(mode, 'r')) exit(-1);
      printf("tell me the full path-name of the file? ");
      scanf("%s", filename);
      if ((fp = (FILE*)fopen(filename, mode)) != NULL)  return(fp);
      puts("Can't find the file.  I give up.");
      exit(-1);
   }
   return(fp);
}


int matout(FILE* fout, double x[], int n, int m)
{
   int i, j;
   fprintf(fout, "\n");
   for (i = 0; i < n; i++) {
      for (j = 0; j < m; j++)
         fprintf(fout, " %11.6f", x[i * m + j]);
      fprintf(fout, "\n");
   }
   return (0);
}

int matout2(FILE* fout, double x[], int n, int m, int wid, int deci)
{
   int i, j;
   fprintf(fout, "\n");
   for (i = 0; i < n; i++) {
      for (j = 0; j < m; j++)
         fprintf(fout, " %*.*g", wid - 1, deci, x[i * m + j]);
      fprintf(fout, "\n");
   }
   return (0);
}


int matIout(FILE* fout, int x[], int n, int m)
{
   int i, j;
   fprintf(fout, "\n");
   for (i = 0; i < n; i++) {
      for (j = 0; j < m; j++)
         fprintf(fout, "  %4d", x[i * m + j]);
      fprintf(fout, "\n");
   }
   return (0);
}

int MultiNomialAliasSetTable(int ncat, double prob[], double F[], int L[], double space[])
{
/* This sets up the tables F and L for the alias algorithm for generating samples from the
   multinomial distribution MN(ncat, p) (Walker 1974; Kronmal & Peterson 1979).

   F[i] has cutoff probabilities, L[i] has aliases.
   I[i] is an indicator: -1 for F[i]<1; +1 for F[i]>=1; 0 if the cell is now empty.

   Should perhaps check whether prob[] sums to 1.
*/
   signed char* I = (signed char*)space;
   int i, j, k, nsmall;

   for (i = 0; i < ncat; i++)  L[i] = -9;
   for (i = 0; i < ncat; i++)  F[i] = ncat * prob[i];
   for (i = 0, nsmall = 0; i < ncat; i++) {
      if (F[i] >= 1)  I[i] = 1;
      else { I[i] = -1; nsmall++; }
   }
   for (i = 0; nsmall > 0; i++) {
      for (j = 0; j < ncat; j++)  if (I[j] == -1) break;
      for (k = 0; k < ncat; k++)  if (I[k] == 1)  break;
      if (k == ncat)  break;

      L[j] = k;
      F[k] -= 1 - F[j];
      if (F[k] < 1) { I[k] = -1; nsmall++; }
      I[j] = 0;  nsmall--;
   }
   return(0);
}

int MultiNomialAlias(int n, int ncat, double F[], int L[], int nobs[])
{
/* This generates multinomial samples using the F and L tables set up before,
   using the alias algorithm (Walker 1974; Kronmal & Peterson 1979).

   F[i] has cutoff probabilities, L[i] has aliases.
   I[i] is an indicator: -1 for F[i]<1; +1 for F[i]>=1; 0 if the cell is now empty.
*/
   int i, k;
   double r;

   for (i = 0; i < ncat; i++)  nobs[i] = 0;
   for (i = 0; i < n; i++) {
      r = rndu(0) * ncat;
      k = (int)r;
      r -= k;
      if (r > F[k]) k = L[k];
      nobs[k]++;
   }
   return (0);
}

int rndDiscreteAlias(int ncat, double F[], int L[])
{
/* This is modified from MultiNomialAlias()
*/
   int k;
   double r;

   r = rndu(0) * ncat;
   k = (int)r;
   r -= k;
   if (r > F[k]) k = L[k];
   return (k);
}
