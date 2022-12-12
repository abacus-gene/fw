#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<stdarg.h>


#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))
#define Pi  3.1415926535897932384626433832795

#define NTHREADS 18
#define PIN_THREADS_CORE

double rndu(int thread_id);
void SetSeed(int seed, int thread_id);
int zero(double x[], int n);
double sum(double x[], int n);
double rndNormal(void);
void zerror(const char* format, ...);
void starttimer(void);
char* printtime(char timestr[]);
FILE* zopen(char* filename, char* mode);
int matout(FILE* fout, double x[], int n, int m);
int matout2(FILE* fout, double x[], int n, int m, int wid, int deci);
int matIout(FILE* fout, int x[], int n, int m);
int MultiNomialAliasSetTable(int ncat, double prob[], double F[], int L[], double space[]);
int MultiNomialAlias(int n, int ncat, double F[], int L[], int nobs[]);
int rndDiscreteAlias(int ncat, double F[], int L[]);

int get_options(char* ctlf); 
int initialize();
void freemem();
void print_zygote(FILE* fout, char* zygote);
void print_gamete(FILE* fout, char* gamete);
