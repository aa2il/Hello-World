/*
  Copyright (c) 2003-2007 Scalable Informatics
 */

#include <stdio.h>
#include <math.h>
#include <sys/timeb.h>
#include <time.h>
#include <stdlib.h>

#include <unistd.h>
#include "string.h"

#define NUMBER_OF_CALIPER_POINTS 10

struct timeb t_initial,t_final,caliper[NUMBER_OF_CALIPER_POINTS];
double genrand(void);

int main(int argc, char **argv)
  {
    int i,j,milestone,rc,iteration,k;
    int DIM=4000;
    double **a,**b,**c,total,sum_a,sum_b,delta_t,dot;
    int true = (1==1), false = (1==0);
    
    printf("D: checking arguments: N_args=%i \n",argc);
    if (argc!=3) {
      printf("\n     USAGE:     %s -n <DIMENSION>\n\n",argv[0]);
    }
    for(i=0;i<argc;i++)
     {
      printf("D: arg[%i] = %s\n",i,argv[i]);
      if (strncmp(argv[i],"-n",2)==0)
       {
	DIM = atoi(argv[i+1]);
 	printf("D: DIM found to be = %i\n",DIM);
	printf("D: should be %s\n",argv[i+1]);	
       }
     }
    printf("DIM=%d\n",DIM);    
    
    /* allocate memory for matrices */
    milestone       = 0;
    printf("Allocating memory - DIM = %d ... \n",DIM);
    rc=ftime(&caliper[milestone]);
    
    /* first allocate a pointer to an array of pointers */
    a	= (double **)calloc(DIM,sizeof(double *));
    b	= (double **)calloc(DIM,sizeof(double *));
    c	= (double **)calloc(DIM,sizeof(double *));

    /* make sure these allocations returned a pointer to
       an array, otherwise quit with an error message  */
    if ((a == NULL) || (b == NULL) || (c == NULL))
       {
         fprintf(stderr,"memory allocation for array or new_array failed\n");
         exit (-3);
       }

    /* second, for each pointer to pointer, 
       allocate an array */
    for(i=0;i<DIM;i++)
       {
         a[i]	= (double *)calloc(DIM,sizeof(double));
         b[i]	= (double *)calloc(DIM,sizeof(double));
         c[i]	= (double *)calloc(DIM,sizeof(double));
	 
	 /* again, if there was an error in any of these
	    allocations, then exit and print a message */
	 if ((a[i] == NULL) || (b[i] == NULL) || (c[i] == NULL))
            {
              fprintf(stderr,"array column allocation did not succeed\n");
              exit(-6);
	    }
       }    
    
    milestone++;
    rc=ftime(&caliper[milestone]);
    total=(double)(DIM*DIM)*sizeof(double)/(1024.0*1024.0);
    printf("array size in MB = %-.3f MB\n (remember, you have 2 of these)",total);

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

    printf("normalization a: %-.5f,  b: %-.5f\n",1.0/sqrt(sum_a),1.0/sqrt(sum_b));
    milestone++;
    rc=ftime(&caliper[milestone]);

   /* matrix multiply 
    *
    *  c[i][j]= a_row[i] dot b_col[j]  for all i,j
    *           a_row[i] -> a[i][0 .. DIM-1]
    *		b_col[j] -> b[0 .. DIM-1][j]
    *
    */
    for(i=0;i<DIM;i++)
     {     
      for(j=0;j<DIM;j++)
       {  
         dot=0.0;      
         for(k=0;k<DIM;k++) dot += a[i][k]*b[k][j];
         c[i][j]=dot;
       }

     }
 
    milestone++;
    rc=ftime(&caliper[milestone]);
    
    /* print a vector to prevent the optimizer from making the calcuation go away... */
    dot=0.0;
    for(i=0;i<DIM-1;i++)
      {
        printf("%-.4f, ",c[i][DIM/2]);
        dot+=c[i][DIM/2];
      }   
     printf("%-.4f\n ",c[DIM-1][DIM/2]);
     printf("%-.9f\n",dot);

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

  

