MACHINE=$(shell gcc -dumpmachine)
#CC=../$(MACHINE)/bin/mpicc
#FC=../$(MACHINE)/bin/mpif77
#EXEC=../$(MACHINE)/bin/mpirun -n 10
CC=mpicc
FC=mpif77
EXEC=mpirun -n 10
INCL=../$(MACHINE)/include
RPATH=../$(MACHINE)/lib

all: helloworld-mpi helloworld-mumps file

clean:
	rm -f helloworld-mpi helloworld-mumps file tags out-test.*

helloworld-mpi: helloworld-mpi.c
	@echo "** TEST: $@"
	${CC} $^ -o $@
	${EXEC} $@

file: file.c
	@echo "** TEST: $@"
	${CC} -I${INCL} -L${RPATH} -Wl,-rpath,${RPATH} -Wall -lbebop_util -lsparse_matrix_converter $^ -o $@

helloworld-mumps: helloworld-mumps.c
	@echo "** TEST: $@"
	# need fortran compiler to compile C that is linking into fortran
	${FC} -L${RPATH} -I${INCL} -Wl,-rpath,${RPATH} $^ -ldmumps -o $@
	${EXEC} $@	

ctags:
	ctags --c-kinds=pe --recurse -R file.c ../x86_64-pc-linux-gnu/include/bebop/
