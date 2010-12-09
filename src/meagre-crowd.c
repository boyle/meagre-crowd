/* Example program using the C interface to the
* double precision version of MUMPS, dmumps_c.
* We solve the system A x = RHS with
* A = diag(1 2) and RHS = [1 4]ˆT
* Solution is [1 2]ˆT */
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <libgen.h>
#include <string.h>

#include <mpi.h>

#include <bebop/util/init.h>
#include <bebop/util/enumerations.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h>
#include <bebop/smc/coo_matrix.h>

#include "args.h"
#include "perftimer.h"
#include "solver_mumps.h"
#include "file.h"




// test result of matrix computations
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


int main(int argc, char ** argv) {
  int ierr, retval;

  // handle command-line arguments
  struct parse_args* args = calloc(1,sizeof(struct parse_args));
  assert(args != NULL); // calloc failure
  if((retval = parse_args(argc, argv, args)) != EXIT_SUCCESS) {
    free(args);
    return retval;
  }
  

  // initialize MPI
  perftimer_t* timer = perftimer_malloc();
  if(args->rep == 0) {
    perftimer_inc(timer,"initialization",-1);
    perftimer_adjust_depth(timer,+1);
    perftimer_inc(timer,"MPI init",-1);
  }
  ierr = MPI_Init(&argc, &argv);               assert(ierr == 0);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &(args->mpi_rank)); assert(ierr == 0);

  
  // Define the problem on the host
  double* b = NULL;
  struct sparse_matrix_t* A;
  if (args->mpi_rank == 0) {

    if(args->rep == 0)
      perftimer_inc(timer,"sparse file handling",-1);
    bebop_default_initialize (argc, argv, &ierr); assert (ierr == 0);

    if(args->rep == 0) {
      perftimer_adjust_depth(timer,-1);
      perftimer_inc(timer,"input",-1);
      perftimer_adjust_depth(timer,+1);
    }

    if(args->rep == 0)
      perftimer_inc(timer,"load",-1);

    if((retval = load_matrix(args->input, &A)) != 0) {
      return retval;
    }

    if(args->rep == 0) {
      perftimer_inc(timer,"solver",-1);
      perftimer_inc(timer,"rhs",-1);
    }

    // verbose output
    if(args->verbosity >= 2) { // an additional line showing MPI/OpenMP info
      int c_mpi, c_omp;
      ierr = MPI_Comm_size(MPI_COMM_WORLD,&c_mpi); assert(ierr == 0);
      c_omp = 0; // TODO omp_get_num_threads();
      printf("MPI cores=%d, openMP threads=%d\n",c_mpi,c_omp);
    }

    if(args->verbosity >= 1) {
      assert(A != NULL);
      int ierr = sparse_matrix_convert(A, COO); assert(ierr == 0);
      struct coo_matrix_t* Acoo = A->repr;
      printf("Ax=b: A is %dx%d, nz=%d, b is %dx%d, nz=%d\n",
             Acoo->m, Acoo->n, Acoo->nnz,
             Acoo->m, 1, Acoo->m);

      if(args->verbosity >= 2) {
        int i;
        int n = Acoo->nnz;
        double* v = Acoo->val;
        for(i=0;i<n;i++) {
          printf("  A(%i,%i)=%.2f\n", Acoo->II[i], Acoo->JJ[i], v[i]);
        }
      }
    }

    if(args->rep == 0)
      perftimer_adjust_depth(timer,-1);
    //destroy_sparse_matrix (A); // TODO can't release it unless we're copying it...

    // allocate an all-zeros right-hand side of A.m rows
    // TODO or load from a file
    
    unsigned int n = matrix_rows(A);
    b = malloc(n*sizeof(double));
    assert(b != NULL); // malloc failure
    { // initialize right-hand-side (b)
      int i;
      for(i=0;i<n;i++)
        b[i] = i;
      if(args->verbosity >= 2) { // show the rhs matrix
        for(i=0;i<n;i++)
          printf("  b(%i,1)=%.2f\n", i, b[i]);
      }
    }
  }

  DMUMPS_STRUC_C* id = solver_init_dmumps(args, timer, NULL); // TODO figure out passing arguments for initialization to all clients (re: A)

  if(args->mpi_rank == 0) {
    solver_data_prep_dmumps(id, A, NULL); // Note: b is set w/in the execution loop
  }

  int i, r = 0;
  do {
    if(r != 0)
      perftimer_restart(&timer);

    if(args->mpi_rank == 0)
      solver_data_prep_dmumps(id, NULL, b); // reset rhs = b

    solver_solve_dmumps(id, args, timer);

    r++;
  } while (r < args->rep);

  if(args->rep == 0)
    perftimer_inc(timer,"clean up",-1);
  solver_finalize_dmumps(id);

  if(args->rep == 0) {
    perftimer_adjust_depth(timer,-1);
    perftimer_inc(timer,"output",-1);
  }
/*  if (args->mpi_rank == 0) {
// TODO this code needs reworking to make it flexible: compare answers with-in some percentage?
    printf("Solution is\n");
    int i;
    for(i=0;i<id.n;i++) {
      printf("  %.2f\n", id.rhs[i]);
      // Octave says the answer should be
      // ans =
      //     1.000000
      //     0.434211
      //     0.250000
      //    -0.092105
      //     0.250000
      //
    }
    // test result
    // TODO make this a command-line option with file to compare
    perftimer_inc(timer,"test",-1);
    const double const e[5] = {1.0, 0.434211, 0.25, -0.092105, 0.25};
    if( results_match(e,id.rhs,id.n,0.001) ) {
      retval = 0; printf("PASS\n");
    }
    else {
      retval = 1; printf("FAIL\n");
    }
  }
*/
  // clean up
  free(b);
  if(args->rep == 0) {
    perftimer_inc(timer,"clean up",-1);
    perftimer_adjust_depth(timer,+1);

    perftimer_inc(timer,"MPI",-1);
  }
  ierr = MPI_Finalize(); assert(ierr == 0);

  // show timing info, if requested, to depth N
  if(args->rep == 0) {
    perftimer_adjust_depth(timer,-1);
    perftimer_inc(timer,"finished",-1);
  }
  if(args->mpi_rank == 0) {
    if(args->timing_enabled == 1) {
      perftimer_printf_csv_header(timer,2);
      perftimer_printf_csv_body(timer,2);
    }
    else if(args->timing_enabled != 0) {
      perftimer_printf(timer,args->timing_enabled-2);
    }
  }
  perftimer_free(timer);
  free(args);
  return retval;
}




