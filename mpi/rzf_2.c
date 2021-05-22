#include <stdio.h>
#include <math.h>
#include <sys/timeb.h>
#include <time.h>
#include <stdlib.h>
/* include MPI definitions */
#include <mpi.h>
/*   */
#include <unistd.h>
#include "string.h"

#define NUMBER_OF_CALIPER_POINTS 10

struct timeb t_initial,t_final,caliper[NUMBER_OF_CALIPER_POINTS],node_timing;
int in_parallel[NUMBER_OF_CALIPER_POINTS];

int main(int argc, char **argv)
  {
    int i,j,milestone,n,rc, nthreads;
    int NMAX=5000000, MMAX=10;
    double *array,total,p_sum,sum,delta_t,pi, pi_exact;
    int true = (1==1), false = (1==0),inf;
    int name_length=0;
    char *cpu_name,c; 

    /* beginning of parallelization bits */
    int loop_min, loop_max, NCPUs, tid;

    /* add in MPI startup routines */
    /* 1st: launch the MPI processes on each node */
    MPI_Init(&argc,&argv);	

    /* 2nd: request a thread id, sometimes called a "rank" from 
            the MPI master process, which has rank or tid == 0 */
    MPI_Comm_rank(MPI_COMM_WORLD, &tid);  

    /* 3rd: this is often useful, get the number of threads
            or processes launched by MPI, this should be NCPUs-1 */
    MPI_Comm_size(MPI_COMM_WORLD, &nthreads);

    
    /* this code is run on EVERY thread */
    milestone	= 0;
    sum		= 0.0;
    p_sum	= 0.0;   /* partial sums, will be used in a bit */
    pi_exact	= 4.0*atan(1.0);
    n		= 2;
    inf         = 1000000000 ;
    rc=ftime(&caliper[milestone]);
    in_parallel[milestone]	= false;

    /*
       Riemann zeta function of integer argument
       (c.f. http://mathworld.wolfram.com/RiemannZetaFunction.html )
       
       zeta(n) = sum[k=1;k<=inf;k++] (1.0/pow(k,n))
    
                  inf 
                  ----
                  \        1
       zeta(n) =   >       -
                  /         n
		  ----     k
                  i=1
    

	this code will compute zeta(2) in order to calculate
	pi.  pi*pi/6 = zeta(2), so pi = sqrt(6*zeta(2))

	to run this code, type

		./rzf.exe -l INFINITY -n n -N NCPUs

	where INFINITY is an integer value of how many terms you would 
	like to take for your sum, and n is the argument to the Reimann
	zeta function.  If you use 2 for n, then it will calculate pi 
	for you as well.

     */
   
     /* remember, we want to enclose as much work as possible
        in the parallel region ... as long as you are careful,
        this can include argument processing.  You need to be
	careful about file IO and related, but argument processing
	should be fine */
 
       printf("D[tid=%i]: checking arguments: N_args=%i \n",tid,argc);
       if (argc!=5 & tid==1) {
         printf("\n     USAGE:     %s -n <N> -l <LARGE_INTEGER>\n\n",argv[0]);
       }
       for(i=0;i<argc;i++)
          {
	    printf("D[tid=%i]: arg[%i] = %s\n",tid,i,argv[i]);
            if (strncmp(argv[i],"-n",2)==0)
               {
		n = atoi(argv[i+1]);
 		printf("D[tid=%i]: n found to be = %i\n",tid,n);
		printf("D[tid=%i]: should be %s\n",tid,argv[i+1]);	
               }
            if (strncmp(argv[i],"-l",2)==0)
               {
		inf = atoi(argv[i+1]);
 		printf("D[tid=%i]: infinity found to be = %i\n",tid,inf);
		printf("D[tid=%i]: should be %s\n",tid,argv[i+1]);
               }
	/* now here we change this, as we are getting NCPUs from 
	   nthreads as above */
	/*
            if (strncmp(argv[i],"-N",2)==0)
               {
                NCPUs = atoi(argv[i+1]);
                printf("D: NCPUs found to be = %i\n",NCPUs);
                printf("D: should be %s\n",argv[i+1]);  
               }
	 */
		NCPUs = nthreads ;
		printf("D[tid=%i]: NCPUs found to be = %i\n",tid,NCPUs);
          }
       printf("D[tid=%i]: N = %d\nINFINITY = %d\n",tid,n,inf);

    /* so far, so good.  Now we start getting tricky ... */
    sum=0.0;
    
    /* it is entirely possible/likely that cpu_name will be
       completely different on each thread 

	First, cpu_name character array is allocated on 
	EVERY thread ... remember, NOTHING is shared, 
	every process is independent.

     */
    cpu_name	= (char *)calloc(80,sizeof(char));
    gethostname(cpu_name,80);
    printf("D[tid=%i]: running on machine = %s\n",tid,cpu_name);
   
    /* this loop is no longer needed, as it is IMPLICIT */
/*    for(tid=0;tid<NCPUs;tid++)   */
    {
    
    milestone++;
    rc=ftime(&caliper[milestone]);
    in_parallel[milestone]	= true;
    loop_min 	= 1 +  (int)((tid + 0)  *  (inf-1)/NCPUs);
    loop_max	=      (int)((tid + 1) 	*  (inf-1)/NCPUs); 
    printf("D[tid=%i]: loop_max=%i, loop_min=%i\n",tid,loop_max,loop_min);
    /* now here is a problem, each thread updates its own copy of 
       p_sum, and its own copy of caliper ...
       how do we get these back to the main process ? */
    for(i=loop_max-1;i>=loop_min;i--)
       {
          p_sum += 1.0/pow(i,(double)n);
       }
    }
    /* we will use something called a REDUCTION */
    MPI_Reduce(&p_sum,&sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);   
    /* this tells the nodes to apply the MPI_SUM operator to the
       values pointed at by &p_sum, where the & operator tells
       the program to return the "address of" the item to its 
       right ... it is a pointer, which we hand to the MPI_Reduce
       call.  From there, MPI_Reduce applies this operator to all
       threads or ranks, and returns the value in the variable 
       pointed to by &sum

	We tell the system that these variables are of MPI_DOUBLE
	type, that the thread requesting this reduction is 0,
	or the master thread, and that we are using a particular
	MPI communicator handle ... think of the last one like a 
	file handle, you need it in all MPI operations

     */
    
    /* 
	but what about the caliper? ... we do the same thing, though
        we will need a little slight of hand as caliper is a structure
	Unfortunately, caliper is a set of time stamps, so we are less
	interested in reducing them as compared to printing them ...
     */

    milestone++;
    rc=ftime(&caliper[milestone]);
    in_parallel[milestone]	= false;
    /* the rest of this should be done only on the master process, tid == 0 */
    
    if (tid == 0) {
    printf("zeta(%i)  = %-.15f \n",n,sum);
    if (n == 2)
     {
       pi = sqrt(sum*6.0);
       printf("pi = %-.15f \n",pi);
    
       printf("error in pi = %-.15f \n",pi_exact-pi);
       printf("relative error in pi = %-.15f \n",(pi_exact-pi)/pi_exact);
     }
    }

    /* now report the milestone time differences */
    for (i=0;i<=(milestone-1);i++)
       {
         delta_t = (double)(caliper[i+1].time-caliper[i].time);
         delta_t += (double)(caliper[i+1].millitm-caliper[i].millitm)/1000.0;
         if ((tid == 0) || (in_parallel[i] == true)) 
 	  {
	   printf("[tid=%i] Milestone %i to %i: time = %-.3fs\n",tid,i,i+1,delta_t);
          }
       }
   
    /* ok, time to tear down the parallel region */
    MPI_Finalize();
  }
