#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define ISPRIME 1
#define NOTPRIME 0

int isPrime(int n)
{
	int squareroot, i;
	if ( n>10 ){
		squareroot = (int)sqrt(n);
		for ( i=3; i<=squareroot; i+=2 ){
			if ( (n%i) == 0 )
				return NOTPRIME;
		}
		return ISPRIME;
	}
	else
		return NOTPRIME;
}

int main(int argc, char *argv[])
{
	int rank, tasks;
	int myStart, gap, primeCount,  countSum, prime, largestPrime;
	int i;
	int maxNumber;

	if ( argc != 2 ){
		printf("prime_mpi MAXNUMBER\n");
		exit(EXIT_FAILURE);
	}

	maxNumber = atoi(argv[1]);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &tasks);
	
	myStart = rank + rank + 1;
	gap = tasks + tasks;
	primeCount = 0;
	prime = 0;

	if ( rank == 0 ){
		primeCount = 4;	/* 2, 3, 5, 7 primes */
		for ( i=myStart; i<=maxNumber; i+=gap ){
			if ( isPrime(i) == ISPRIME ){
				primeCount++;
				prime = i;
			}
		}
                MPI_Reduce(&prime, &largestPrime, 1, MPI_LONG, MPI_MAX,
                                0, MPI_COMM_WORLD);
                MPI_Reduce(&primeCount, &countSum, 1, MPI_LONG, MPI_SUM,
                                0, MPI_COMM_WORLD);
		printf("The largest prime is %d and total primes is %d\n",
				largestPrime, countSum);

	}
	else{
		for ( i=myStart; i<=maxNumber; i+=gap ){
			if ( isPrime(i) == ISPRIME ){
				primeCount++;
				prime = i;
			}
		}
		MPI_Reduce(&prime, &largestPrime, 1, MPI_LONG, MPI_MAX,
				0, MPI_COMM_WORLD);
		MPI_Reduce(&primeCount, &countSum, 1, MPI_LONG, MPI_SUM,
				0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return (0);
}
