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

int main(int argc, char **argv)
  {
    int i,j,k,milestone,n,rc,inf;
    int NMAX=5000000, MMAX=10;
    double *array,total,p_sum,sum,delta_t,pi, pi_exact;
    int true = (1==1), false = (1==0);
    int name_length=0;
    char *cpu_name,c; 

    milestone	= 0;
    sum		= 0.0;
    pi_exact	= 4.0*atan(1.0);
    n		= 2;
    inf         = 1000000000 ;
    
    rc=ftime(&caliper[milestone]);

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
                  k=1
    

	this code will compute zeta(2) in order to calculate
	pi.  pi*pi/6 = zeta(2), so pi = sqrt(6*zeta(2))

	to run this code, type

		./rzf.exe -l INFINITY -n n

	where INFINITY is an integer value of how many terms you would 
	like to take for your sum, and n is the argument to the Reimann
	zeta function.  If you use 2 for n, then it will calculate pi 
	for you as well.

     */
    
       printf("D: checking arguments: N_args=%i \n",argc);
       if (argc!=5) {
         printf("\n     USAGE:     %s -n <N> -l <LARGE_INTEGER>\n\n",argv[0]);
       }
       for(i=0;i<argc;i++)
          {
	    printf("D: arg[%i] = %s\n",i,argv[i]);
            if (strncmp(argv[i],"-n",2)==0)
               {
		n = atoi(argv[i+1]);
 		printf("D: N found to be = %i\n",n);
		printf("D: should be %s\n",argv[i+1]);	
               }
            if (strncmp(argv[i],"-l",2)==0)
               {
		inf = atoi(argv[i+1]);
 		printf("D: infinity found to be = %i\n",inf);
		printf("D: should be %s\n",argv[i+1]);
               }
          }
       printf("N = %d\nINFINITY = %d\n",n,inf);

    sum=0.0;
    cpu_name	= (char *)calloc(80,sizeof(char));
    gethostname(cpu_name,80);
    printf("D: running on machine = %s\n",cpu_name);
   
    milestone++;
    rc=ftime(&caliper[milestone]);

 
    for(k=inf-1;k>=1;k--)
       {
          sum += 1.0/pow(k,(double)n);
       }

    milestone++;
    rc=ftime(&caliper[milestone]);

    printf("zeta(%i)  = %-.15f \n",n,sum);
    if (n == 2)
     {
       pi = sqrt(sum*6.0);
       printf("pi = %-.15f \n",pi);
    
       printf("error in pi = %-.15f \n",pi_exact-pi);
       printf("relative error in pi = %-.15f \n",(pi_exact-pi)/pi_exact);
     }
   
    
    /* now report the milestone time differences */
    for (i=0;i<=(milestone-1);i++)
       {
         delta_t = (double)(caliper[i+1].time-caliper[i].time);
         delta_t += (double)(caliper[i+1].millitm-caliper[i].millitm)/1000.0;
	 printf("Milestone %i to %i: time = %-.3fs\n",i,i+1,delta_t);
       }
   

  }
