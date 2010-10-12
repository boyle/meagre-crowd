MACHINE=$(shell gcc -dumpmachine)
CC=../$(MACHINE)/bin/mpicc
FC=../$(MACHINE)/bin/mpif77
EXEC=../$(MACHINE)/bin/mpirun -n 10
INCL=../$(MACHINE)/include
RPATH=../$(MACHINE)/lib

all: helloworld-mpi helloworld-mumps

clean:
	rm -f helloworld-mpi helloworld-mumps

helloworld-mpi: helloworld-mpi.c
	@echo "** TEST: $@"
	${CC} $^ -o $@
	${EXEC} $@

helloworld-mumps: helloworld-mumps.c
	@echo "** TEST: $@"
	# need fortran compiler to compile C that is linking into fortran
	${FC} -I${INCL} -Wl,-rpath,${RPATH} $^ -ldmumps -o $@
	${EXEC} $@	

