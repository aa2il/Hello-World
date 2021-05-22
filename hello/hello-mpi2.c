#include "stdio.h"
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
  int tid,nthreads;
  char *cpu_name;
  double time_initial,time_current,time;
  
  /* add in MPI startup routines */
  /* 1st: launch the MPI processes on each node */
  MPI_Init(&argc,&argv);
  time_initial  = MPI_Wtime();

  /* 2nd: request a thread id, sometimes called a "rank" from
     the MPI master process, which has rank or tid == 0   */
  MPI_Comm_rank(MPI_COMM_WORLD, &tid);
  
  /* 3rd: this is often useful, get the number of threads
     or processes launched by MPI, this should be NCPUs-1   */
  MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
  
  /* on EVERY process, allocate space for the machine name */
  cpu_name    = (char *)calloc(80,sizeof(char));

  /* get the machine name of this particular host ... well
     at least the first 80 characters of it ... */
  gethostname(cpu_name,80);
  time_current  = MPI_Wtime();
  time  = time_current - time_initial;
  printf("%.3f tid=%i : hello MPI user: machine=%s [NCPU=%i]\n",
         time, tid, cpu_name, nthreads);

  MPI_Finalize();
  return(0);
}
