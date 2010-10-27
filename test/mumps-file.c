/* Example program using the C interface to the
* double precision version of MUMPS, dmumps_c.
* We solve the system A x = RHS with
* A = diag(1 2) and RHS = [1 4]ˆT
* Solution is [1 2]ˆT */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "mpi.h"
#include "dmumps_c.h"

#include <bebop/util/init.h>
#include <bebop/util/enumerations.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h>
#include <bebop/smc/coo_matrix.h>

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
int main(int argc, char ** argv) {
  DMUMPS_STRUC_C id;
  int myid, ierr;

  ierr = MPI_Init(&argc, &argv);               assert(ierr == 0);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid); assert(ierr == 0);
  
  /* Initialize a MUMPS instance. Use MPI_COMM_WORLD. */
  id.job=JOB_INIT; id.par=1; id.sym=0; id.comm_fortran=USE_COMM_WORLD;
  dmumps_c(&id);
  /* Define the problem on the host */
  if (myid == 0) {
    bebop_default_initialize (argc, argv, &ierr); assert (ierr == 0);

    struct sparse_matrix_t* A;
    A = load_sparse_matrix (MATRIX_MARKET, "test.mm"); assert (A != NULL);
    ierr = sparse_matrix_convert(A, COO); assert(ierr == 0);

    // check somethings are as expected
    struct coo_matrix_t* Acoo = A->repr;
    assert(Acoo->index_base == ZERO); // index zero is the first entry
    assert(Acoo->symmetry_type == UNSYMMETRIC);
    assert(Acoo->value_type == REAL); // don't handle complex... yet TODO

    // TODO do soemthing with A.ownership, so we can tell bebop to clean itself up, but not have to copy the elements

    // mumps: irn=row indices, jcn=column idices, a=values, rhs=right-hand side, n = matrix order (on-a-side?) nz=non-zeros?
    id.n   = Acoo->m; // A.m: rows, A.n: columns
    id.nz  = Acoo->nnz; // non-zeros
    id.irn = Acoo->II; // row    indices
    id.jcn = Acoo->JJ; // column indices
    id.a   = Acoo->val;
    // show what we got
    printf("Input A is\n");
    int i;
    for(i=0;i<id.nz;i++) {
      printf("  (%i,%i)=%.2f\n",id.irn[i],id.jcn[i],id.a[i]);
    }
    // allocate an all-zeros right-hand side of A.m rows
    id.rhs = malloc(id.n * sizeof(double)); assert(id.rhs != NULL);
    printf("Right-hand side is\n");
    for(i=0;i<id.n;i++) {
      id.rhs[i] = i;
      printf("  %.2f\n", id.rhs[i]);
    }
  
    //destroy_sparse_matrix (A); // TODO can't release it unless we're copying it...
  }
  #define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
  /* No outputs */
  id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
  /* Call the MUMPS package. */
  id.job=6;       dmumps_c(&id);
  id.job=JOB_END; dmumps_c(&id); /* Terminate instance */
  if (myid == 0) {
    printf("Solution is\n");
    int i;
    for(i=0;i<id.n;i++) {
      printf("  %.2f\n", id.rhs[i]);
      /* Octave says the answer should be
       * ans =
       *	   1.000000
       *	   0.434211
       *	   0.250000
       *	  -0.092105
       *	   0.250000
       */
    }
    printf("FAIL\n");
  }
  ierr = MPI_Finalize();
  return 10;
}
