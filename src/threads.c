#include "fw.h"

extern int debug;
extern int nthreads, thread_start;
thread_data_t thread_data[NTHREADS];

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
