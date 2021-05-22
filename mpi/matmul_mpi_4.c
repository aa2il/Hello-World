// Lets see if we can do the unrolling correctly & generally
// The final answer is probablt not correct bx the final gather is hosed.

#include <stdio.h>
#include <math.h>
#include <sys/timeb.h>
#include <time.h>
#include <stdlib.h>

#include <unistd.h>
#include "string.h"
/* include MPI definitions */
#include <mpi.h>
/*   */
#define NUMBER_OF_CALIPER_POINTS 10
#define UNROLL 1

struct timeb t_initial,t_final,caliper[NUMBER_OF_CALIPER_POINTS];
double genrand(void);
double round(double x);

int main(int argc, char **argv)
  {
    int i,j,milestone,rc,iteration,k,istep;
    int DIM=1000,DIM2;
    double **a,**b,**c,total,sum_a,sum_b,delta_t,sanity,dot[4];
    double dot0;
    int true = (1==1), false = (1==0);
        /* beginning of parallelization bits */
    int loop_min, loop_max, NCPUs, tid, nthreads, number_transfered;
    char *cpu_name; 
    
    MPI_Status stat;
    
    /* add in MPI startup routines */
    /* 1st: launch the MPI processes on each node */
    MPI_Init(&argc,&argv);	

    /* 2nd: request a thread id, sometimes called a "rank" from 
            the MPI master process, which has rank or tid == 0 */
    MPI_Comm_rank(MPI_COMM_WORLD, &tid);  

    /* 3rd: this is often useful, get the number of threads
            or processes launched by MPI, this should be NCPUs-1 */
    MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
    NCPUs = nthreads ;
    
    cpu_name	= (char *)calloc(80,sizeof(char));
    gethostname(cpu_name,80);
    printf("D[tid=%i]: running on machine = %s\n",tid,cpu_name);

    printf("D: checking arguments: N_args=%i \n",argc);
    if (argc!=3 & tid==1) {
      printf("\n     USAGE:     %s -n <DIMENSION>\n\n",argv[0]);
    }
    for(i=0;i<argc;i++)
     {
       printf("D[%d]: arg[%i] = %s\n",tid,i,argv[i]);
      if (strncmp(argv[i],"-n",2)==0)
       {
        DIM = atoi(argv[i+1]);
        printf("D[%d]: DIM found to be = %i\n",tid,DIM);
        printf("D[%d]: should be %s\n",tid,argv[i+1]);  
       }
     }
    printf("D[%d]: DIM=%d\n",tid,DIM);    
    
    /* allocate memory for matrices */
    milestone       = 0;
    printf("D[%d]: Allocating memory (milestone %d) ... \n",tid,milestone);
    rc=ftime(&caliper[milestone]);
    
    /* first allocate a pointer to an array of pointers */
    DIM2=DIM + 3*UNROLL;             // Add some space to overrun when we unroll
    a   = (double **)calloc(DIM2,sizeof(double *));   
    b   = (double **)calloc(DIM2,sizeof(double *));
    
    /* make sure these allocations returned a pointer to
       an array, otherwise quit with an error message  */
    if ((a == NULL) || (b == NULL) )
       {
         fprintf(stderr,"memory allocation for array or new_array failed\n");
         exit (-3);
       }

    /* second, for each pointer to pointer, 
       allocate an array */
    for(i=0;i<DIM2;i++)
      {
        a[i]   = (double *)calloc(DIM,sizeof(double));
        b[i]   = (double *)calloc(DIM,sizeof(double));
         
         /* again, if there was an error in any of these
            allocations, then exit and print a message */
         if ((a[i] == NULL) || (b[i] == NULL) )
            {
              fprintf(stderr,"array column allocation did not succeed\n");
              exit(-6);
            }
       }    

    // No - this does not guarentee that the data is contiguous -should instead alloc a 
    // big block & compute index nto each row
    c   = (double **)calloc(DIM2,sizeof(double *));
    if (c == NULL)
       {
         fprintf(stderr,"memory allocation for c array or new_array failed\n");
         exit (-3);
       }
#if 0
    for(i=0;i<DIM2;i++)
       {      
         c[i]   = (double *)calloc(DIM,sizeof(double));
         
         /* again, if there was an error in any of these
            allocations, then exit and print a message */
         if (c[i] == NULL)
            {
              fprintf(stderr,"array column allocation did not succeed\n");
              exit(-6);
	    }
       }    
#else
    c[0]  = (double *)calloc(DIM2*DIM,sizeof(double));
    if (c[0] == NULL)
      {
         fprintf(stderr,"2D c array allocation failed\n");
         exit(-6);
       }
    for(i=1;i<DIM2;i++)
      c[i]   = c[i-1] + DIM;
#endif

    milestone++;
    rc=ftime(&caliper[milestone]);
    total=(double)(DIM*DIM)*sizeof(double)/(1024.0*1024.0);
    printf("array size in MB = %-.3f MB\n (remember, you have 2 of these - milestone %d)",total,milestone);

    /* put random numbers into a and b, clear c 
       This is the SPIF we will allow, we will explain later on
     */
    for(i=0;i<DIM;i++)
      for(j=0;j<DIM;j++)
       {        
         a[i][j]=genrand();
         b[i][j]=genrand();
	 c[i][j]=0.0;
       }
    milestone++;
    rc=ftime(&caliper[milestone]);
    printf("Filling with random data (milestone %d) \n",milestone);

    /* normalize the rows of a and the columns of b */
     for(i=0;i<DIM;i++)
      {
       sum_a=0.0;
       for(j=0;j<DIM;j++)
        {        
         sum_a += a[i][j]*a[i][j];
         sum_b += b[j][i]*b[j][i];
        }
       for(j=0;j<DIM;j++)
        {        
         a[i][j]/=sqrt(sum_a);
	 b[j][i]/=sqrt(sum_b);
        }       
      }

    milestone++;
    rc=ftime(&caliper[milestone]);
    printf("normalization a: %-.5f,  b: %-.5f (milestone %d)\n",1.0/sqrt(sum_a),1.0/sqrt(sum_b),milestone);

   /* matrix multiply 
    *
    *  c[i][j]= a_row[i] dot b_col[j]  for all i,j
    *           a_row[i] -> a[i][0 .. DIM-1]
    *		b_col[j] -> b[0 .. DIM-1][j]
    *
    */
    
    
    /*
        Each process will operate on a "panel" of i = loop_min .. loop_max
	in the a matrix.  Every process will iterate over the full b
	matrix, and calculate their panel of the c matrix.
	
	Notice that no messages have been sent yet
     */
    loop_min	=(int)((long)(tid + 0) *  (long)(DIM)/(long)NCPUs);
    loop_max    =(int)((long)(tid + 1) *  (long)(DIM)/(long)NCPUs);
    printf("%i : loop_min = %i, loop_max = %i \n",tid,loop_min,loop_max);
   
    if(UNROLL)
      istep=4;
    else
      istep=1;
    for(i=loop_min;i<loop_max;i+=istep)
     {     
      for(j=0;j<DIM;j++)
       {  
#if UNROLL
         // Unroll but be very careful!!  Author of original code obviously wasn't
         // This does seem to make a big difference in execution time (2x) so its worth getting it right
         dot[0]=dot[1]=dot[2]=dot[3]=0.0;     
         for(k=0;k<DIM;k++) 
          {
            dot[0] += a[i+0][k]*b[k][j];
            dot[1] += a[i+1][k]*b[k][j];
            dot[2] += a[i+2][k]*b[k][j];
            dot[3] += a[i+3][k]*b[k][j];
	  }
         c[i+0][j]=dot[0];
         c[i+1][j]=dot[1];
         c[i+2][j]=dot[2];
         c[i+3][j]=dot[3];
#else
         // Dont unroll the loop - let the compiler do this!
         dot0=0.0;      
         for(k=0;k<DIM;k++) 
           dot0 += a[i][k]*b[k][j];
         c[i][j]=dot0;
#endif
       }

     }
 
    milestone++;
    rc=ftime(&caliper[milestone]);
    printf("D[%d] : Done with matrix multiply (milestone %d)\n",tid,milestone);
    
    /* 
       now the master process has to get (gather) the c matrix 
       back to itself.  Rather than using the most efficient method
       with the MPI_Gather call, we will use Send and Recv pairs to send
       rows back to the master process.       
       
       First, lets wait until all processes hit the same spot
     */
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* then loop over each process, and send the indices to the master */
    for(i=1;i<NCPUs;i++)
     {
       /* have the master process receive the row */
       // This probably doesn't work correctly bx, in C, we have no assurance that 2D data is
       // stored contiguously
       loop_min =(int)((long)(i + 0) *  (long)(DIM)/(long)NCPUs);
       loop_max =(int)((long)(i + 1) *  (long)(DIM)/(long)NCPUs);
       number_transfered = ((loop_max ) - loop_min ) * DIM;   // +1;       // Not sure why +1 ???????
       if(tid==1)
         printf("D[%d]: Gathering from %d - min,max,n=\t%d\t%d\t%d\n",tid,i,loop_min,loop_max,number_transfered);
       
       /* receive the rows from the i_th_ process to the master process */
       if (tid == 0)
        {
          MPI_Recv(&c[loop_min][0],number_transfered,
                MPI_DOUBLE,i,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);           
        }
        
       /* send the rows from the i_th_ process to the master process */ 
       if (tid == i)
        {
         MPI_Send(&c[loop_min][0],number_transfered,MPI_DOUBLE,0,i,MPI_COMM_WORLD);
        }
     }
    
 
    milestone++;
    rc=ftime(&caliper[milestone]);
    printf("D[%d]: Done gathering (milestone %d)\n",tid,milestone);
    
    /* print a vector to prevent the optimizer from making the calcuation go away... */
    dot0=0.0;
    
    /* only run this on the master process */
    if (tid == 0)
     {
       for(i=0;i<DIM-1;i++)
         {
           printf("%-.4f, ",c[i][DIM/2]);
           dot0+=c[i][DIM/2];
         }   
       printf("%-.4f\n ",c[DIM-1][DIM/2]);
       printf("%-.9f\n",dot0);
       
       milestone++;
       rc=ftime(&caliper[milestone]);

        /* now report the milestone time differences */
        for (i=0;i<milestone;i++)
           {
             delta_t = (double)(caliper[i+1].time-caliper[i].time);
             delta_t += (double)(caliper[i+1].millitm-caliper[i].millitm)/1000.0;
             printf("milestone %i to %i time=%-.3f seconds\n",i,i+1,delta_t);
           }
      }
   
    MPI_Finalize();
  }

#define N 25
#define M 7

double
genrand()
{
    unsigned long y;
    static int k = 0;
    static unsigned long x[N]={ /* initial 25 seeds, change as you wish */
        0x95f24dab, 0x0b685215, 0xe76ccae7, 0xaf3ec239, 0x715fad23,
        0x24a590ad, 0x69e4b5ef, 0xbf456141, 0x96bc1b7b, 0xa7bdf825,
        0xc1de75b7, 0x8858a9c9, 0x2da87693, 0xb657f9dd, 0xffdc8a9f,
        0x8121da71, 0x8b823ecb, 0x885d05f5, 0x4e20cd47, 0x5a9ad5d9,
        0x512c0c03, 0xea857ccd, 0x4cc1d30f, 0x8891a8a1, 0xa6b7aadb
    };
    static unsigned long mag01[2]={
        0x0, 0x8ebfd028 /* this is magic vector `a', don't change */
    };
    if (k==N) { /* generate N words at one time */
      int kk;
      for (kk=0;kk<N-M;kk++) {
        x[kk] = x[kk+M] ^ (x[kk] >> 1) ^ mag01[x[kk] % 2];
      }
      for (; kk<N;kk++) {
        x[kk] = x[kk+(M-N)] ^ (x[kk] >> 1) ^ mag01[x[kk] % 2];
      }
      k=0;
    }
    y = x[k];
    y ^= (y << 7) & 0x2b5b2500; /* s and b, magic vectors */
    y ^= (y << 15) & 0xdb8b0000; /* t and c, magic vectors */
    y &= 0xffffffff; /* you may delete this line if word size = 32 */
/*
   the following line was added by Makoto Matsumoto in the 1996 version
   to improve lower bit's corellation.
   Delete this line to o use the code published in 1994.
*/
    y ^= (y >> 16); /* added to the 1994 version */
    k++;
    return( (double) y / (unsigned long) 0xffffffff);
}

  

