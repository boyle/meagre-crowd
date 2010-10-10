MACHINE=$(shell gcc -dumpmachine)
CC=../$(MACHINE)/bin/mpicc
EXEC=../$(MACHINE)/bin/mpirun -np 10
all: helloworld-mpi

clean:
	rm -f helloworld-mpi

%: %.c
	@echo "** TEST: $@"
	${CC} $^ -o $@
	${EXEC} $@
