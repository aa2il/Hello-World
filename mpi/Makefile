# Language 	Script Name 	Underlying Compiler
#
#	C 	mpicc 			gcc
#		mpigcc 			gcc
#		mpiicc 			icc
#		mpipgcc 		pgcc
#
#	C++ 	mpiCC 			g++
#		mpig++ 			g++
#		mpiicpc 		icpc
#		mpipgCC 		pgCC
#
#	Fortran mpif77 			g77
#               mpif90                  also seems to point to gfortran
#		mpigfortran 		gfortran
#		mpiifort 		ifort
#		mpipgf77 		pgf77
#		mpipgf90 		pgf90
#
##########################################################################
#
# Many of these codes are from  http://www.linux-mag.com/id/5759/
#
# Single threaded programs:
#
#       hello
#       hello-f
#	time rzf -n 2 -l 1000000000
#	time matmul -n 1000
#
# How to run multi-threaded programs:
#
#	mpirun -np 4 hello-mpi
#	mpirun -np 4 hello-mpi-2
#	mpirun -np 4 rzf_2 -n 2 -l 1000000000
#	mpirun -np 4 rzf_3 -n 2 -l 1000000000
#	mpirun -np 4 matmul_mpi_2 -n 1000
#	mpirun -np 4 matmul_mpi_3 -n 1000
#	mpirun -np 4 matmul_mpi_4 -n 1000
#	mpirun -n 4 mpi_hello_world
#
# Some other hints:
#
# 1)	mpirun --gdb.
# 2)   	mpirun -np 4 ddd my_program    ---> Gives a different ddd window for each process
#
# 3)   	a. Get the program running
#	b. ps ax | grep my_prog         --- Get process IDs
#       c. In source, insert:     while (DebugWait) ;   ---> Make prog wait for you
#	d. gdb my_prog process_number   ---> Attach gdb to process
#	e.  (gdb) set DebugWait = 0     ---> Get prog on its way
#
##########################################################################

FC	= mpif90
FFLAGS	= -g 

CC	= mpicc 
CFLAGS	= -g -O0 

LD	= mpicc 
LDFLAGS	= -g 

PROGRAMS = hello hello-f rzf matmul \
		hello-mpi hello-mpi-2 rzf_2 rzf_3 \
		matmul_mpi matmul_mpi_2 matmul_mpi_3 matmul_mpi_4 \
		prime_mpi mpi_hello_world

##########################################################################

all:	${PROGRAMS}

# A few single threaded programs to check out gcc & gfortran
hello:	hello.c
	gcc hello.c -o hello

hello-f:	hello.f
	gfortran hello.f -o hello-f

rzf:		rzf.c
	gcc rzf.c -o rzf -lm

matmul:		matmul.c
	gcc matmul.c -o matmul -lm

##########################################################################

# Multi-threaded programs
mpi_hello_world: mpi_hello_world.o
	${LD} ${LDFLAGS} $< -o $@ 

matmul_mpi_4:		matmul_mpi_4.o
	${LD} ${LDFLAGS} $< -o $@ -lm

matmul_mpi_3:		matmul_mpi_3.o
	${LD} ${LDFLAGS} $< -o $@ -lm

matmul_mpi_2:		matmul_mpi_2.o
	${LD} ${LDFLAGS} $< -o $@ -lm

matmul_mpi:		matmul_mpi.o
	${LD} ${LDFLAGS} $< -o $@ -lm

rzf_3:		rzf_3.o
	${LD} ${LDFLAGS} $< -o $@ -lm

rzf_2:		rzf_2.o
	${LD} ${LDFLAGS} $< -o $@ -lm

prime_mpi:		prime_mpi.o
	${LD} ${LDFLAGS} $< -o $@ -lm

hello-mpi:		hello-mpi.o
	${LD} ${LDFLAGS} $< -o $@

hello-mpi-2:		hello-mpi-2.o
	${LD} ${LDFLAGS} $< -o $@

##########################################################################

#.c.o:
#	${CC} ${CFLAGS} -c $*.c -o $*.o

.c.o:
	$(CC) -c $(CFLAGS) $<

.f.o:
	$(F77) -c $(FFLAGS) $<

clean:
	rm -f *.o ${PROGRAMS}
