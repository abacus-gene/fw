#if(defined(__linux__))
#define _GNU_SOURCE
#include<sched.h>
#include<unistd.h>
#endif
#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<stdarg.h>

#define NTHREADS 18
#define PIN_THREADS_CORE
#include<pthread.h>

typedef struct thread_data_s {
   pthread_t thread;
   pthread_mutex_t mutex;
   pthread_cond_t cond;

   /* work: -1: end; 0: idle waiting for work; 1: meiosis */
   volatile int work;
   int id, ind_start, nind;
   char* gamete[2];
}  thread_data_t;

#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))
#define Pi  3.1415926535897932384626433832795

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

void pin_to_core(int t);
void* thread_worker(void* arg);
void threads_init(void);
void threads_wakeup(int work, void* data);
void threads_exit();

int get_options(char* ctlf);
int initialize();
void freemem();
int meiosis(char* gamete, char* zygote, int thread_id);
void reproduction(int ind_start, int nind, int thread_id);
int individual_fitness(int curr_gen);
int update_pop_features(int curr_gen);
int print_pop_features(FILE* fout, int curr_gen);
void print_zygote(FILE* fout, char* zygote);
void print_str(FILE* fout, char* str);
