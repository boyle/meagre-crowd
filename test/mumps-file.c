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

int results_match(const double * const expected, const double const * result, const int n, const double precision);
int results_match(const double * const expected, const double const * result, const int n, const double precision) {
  int i;
  for(i=0;i<n;i++) {
    if( (result[i] < (expected[i] - precision)) ||
        (result[i] > (expected[i] + precision)) ) {
      return 0;
    }
  }
  return 1;
}


#define JOB_INIT -1
#define JOB_END -2
#define JOB_ANALYSE 1
#define JOB_FACTORIZE 2
#define JOB_SOLVE 3
#define JOB_ALL 6
#define USE_COMM_WORLD -987654
int main(int argc, char ** argv) {
  DMUMPS_STRUC_C id;
  int myid, ierr;
  int retval = 0;

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
    coo_c_to_fortran(Acoo); assert(Acoo != NULL);
    assert(Acoo->index_base == ONE); // index zero is the first entry
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
  if(1) { // no debug
    id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
  }
  else { // debug
    id.ICNTL(1)=6; // err output stream
    id.ICNTL(2)=6; // warn/info output stream
    id.ICNTL(3)=6; // global output stream
    id.ICNTL(4)=4; // debug level 0:none, 1: err, 2: warn/stats 3:diagnostics, 4:parameters
  }
  /* Call the MUMPS package. */
  id.job=JOB_ALL; dmumps_c(&id);
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
    // test result
    const double const e[5] = {1.0, 0.434211, 0.25, -0.092105, 0.25};
    if( results_match(e,id.rhs,id.n,0.001) ) {
      retval = 0; printf("PASS\n");
    }
    else {
      retval = 1; printf("FAIL\n");
    }
  }
  ierr = MPI_Finalize(); assert(ierr == 0);
  return retval;
}
