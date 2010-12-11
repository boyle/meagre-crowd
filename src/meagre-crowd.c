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
#include "solver_umfpack.h"
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
  const int false = 0;
  const int extra_timing = 0; // allow extra timing: initialization, matrix load, etc.

  // handle command-line arguments
  struct parse_args* args = calloc(1,sizeof(struct parse_args));
  assert(args != NULL); // calloc failure
  if((retval = parse_args(argc, argv, args)) != EXIT_SUCCESS) {
    free(args);
    return retval;
  }


  // initialize MPI
  perftimer_t* timer = perftimer_malloc();
  if(extra_timing && args->rep == 0) {
    perftimer_inc(timer,"initialization",-1);
    perftimer_adjust_depth(timer,+1);
    perftimer_inc(timer,"MPI init",-1);
  }
  ierr = MPI_Init(&argc, &argv);               assert(ierr == 0);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &(args->mpi_rank)); assert(ierr == 0);


  // Define the problem on the host
  double* b = NULL;
  struct sparse_matrix_t* A;
  unsigned int m = 0; // rows
  if (args->mpi_rank == 0) {

    if(extra_timing && args->rep == 0)
      perftimer_inc(timer,"sparse file handling",-1);
    bebop_default_initialize (argc, argv, &ierr); assert (ierr == 0);

    if(extra_timing && args->rep == 0) {
      perftimer_adjust_depth(timer,-1);
      perftimer_inc(timer,"input",-1);
      perftimer_adjust_depth(timer,+1);
    }

    if(extra_timing && args->rep == 0)
      perftimer_inc(timer,"load",-1);

    if((retval = load_matrix(args->input, &A)) != 0) {
      return retval;
    }

    if(extra_timing && args->rep == 0) {
      perftimer_inc(timer,"solver",-1);
      perftimer_inc(timer,"rhs",-1);
    }

    // verbose output
    if(args->verbosity >= 1) {
      assert(A != NULL);
      int c_mpi, c_omp;
      ierr = MPI_Comm_size(MPI_COMM_WORLD,&c_mpi); assert(ierr == 0);
      c_omp = 0; // TODO omp_get_num_threads();
      int ierr = sparse_matrix_convert(A, COO); assert(ierr == 0);
      struct coo_matrix_t* Acoo = A->repr;
      const char* sym, *location, *type;
      switch(Acoo->symmetry_type) {
	case UNSYMMETRIC:    sym = "unsymmetric"; break;
	case SYMMETRIC:      sym = "symmetric"; break;
	case SKEW_SYMMETRIC: sym = "skew symmetric"; break;
	case HERMITIAN:      sym = "hermitian"; break;
	default: assert(false); // fell through
      }
      if((args->verbosity < 2) || (Acoo->symmetry_type == UNSYMMETRIC)) {
	location = "";
      }
      else {
	switch(Acoo->symmetric_storage_location) {
	  case UPPER_TRIANGLE: location = " (upper)"; break;
	  case LOWER_TRIANGLE: location = " (lower)"; break;
	  default: assert(false);
	}
      }
      switch(Acoo->value_type) {
	case REAL:    type = "real"; break;
	case COMPLEX: type = "complex"; break;
	case PATTERN: type = "pattern"; break;
	default: assert(false); // fell through
      }
      const char* solver;
      switch(args->solver) {
        case MUMPS:   solver = "mumps"; break;
        case UMFPACK: solver = "umfpack"; break;
	default: assert(false); // fell through
      }
      printf("Ax=b: A is %dx%d, nz=%d, %s%s, %s, b is %dx%d, nz=%d\nsolved with %s on %d core%s, %d thread%s\n",
             Acoo->m, Acoo->n, Acoo->nnz,
	     sym, location, type,
             Acoo->m, 1, Acoo->m,
	     solver,
	     c_mpi, c_mpi==1?"":"s",c_omp,c_omp==1?"":"s");

      if(args->verbosity >= 3) {
        int i;
        int n = Acoo->nnz;
        double* v = Acoo->val;
        for(i=0;i<n;i++) {
          printf("  A(%i,%i)=%.2f\n", Acoo->II[i], Acoo->JJ[i], v[i]);
        }
      }
    }

    if(extra_timing && args->rep == 0)
      perftimer_adjust_depth(timer,-1);
    //destroy_sparse_matrix (A); // TODO can't release it unless we're copying it...

    // allocate an all-zeros right-hand side of A.m rows
    // TODO or load from a file

    m = matrix_rows(A);
    b = malloc(m*sizeof(double));
    assert(b != NULL); // malloc failure
    { // initialize right-hand-side (b)
      int i;
      for(i=0;i<m;i++)
        b[i] = i;
      if(args->verbosity >= 3) { // show the rhs matrix
        for(i=0;i<m;i++)
          printf("  b(%i)=%.2f\n", i, b[i]);
      }
    }
  }

  DMUMPS_STRUC_C* mumps_p = NULL;
  solve_system_umfpack_t* umfpack_p = NULL;
  switch(args->solver) {
    case MUMPS:
      mumps_p = solver_init_dmumps(args, timer, NULL); // TODO figure out passing arguments for initialization to all clients (re: A)
      break;
    case UMFPACK:
      umfpack_p = solver_init_umfpack(args, timer, NULL); // TODO figure out passing arguments for initialization to all clients (re: A)
      break;
    default: assert(false); // should have caught unknown solver before now!
  }

  if(args->mpi_rank == 0) {
    switch(args->solver) { // Note: b is set w/in the execution loop
      case MUMPS:
	solver_data_prep_dmumps(mumps_p, A, NULL);
	break;
      case UMFPACK:
	solver_data_prep_umfpack(umfpack_p, A, NULL);
	break;
      default: assert(false); // should have caught unknown solver before now!
    }
  }

  int r = 0;
  double* rhs = NULL;
  do {
    if(r != 0)
      perftimer_restart(&timer);

    if(args->mpi_rank == 0) {
      switch(args->solver) {
        case MUMPS:
          solver_data_prep_dmumps(mumps_p, NULL, b); // reset rhs = b
          break;
        case UMFPACK:
          solver_data_prep_umfpack(umfpack_p, NULL, b); // reset rhs = b
          break;
        default: assert(false); // should have caught unknown solver before now!
      }
    }

    switch(args->solver) {
      case MUMPS:
        rhs = solver_solve_dmumps(mumps_p, args, timer);
        break;
      case UMFPACK:
        rhs = solver_solve_umfpack(umfpack_p, args, timer);
        break;
      default: assert(false); // should have caught unknown solver before now!
    }

    r++;
  } while (r < args->rep);

  if(extra_timing && args->rep == 0) {
    perftimer_inc(timer,"clean up",-1);
    perftimer_adjust_depth(timer,-1);
    perftimer_inc(timer,"output",-1);
  }

  if((args->mpi_rank == 0) && (args->verbosity >= 2)) {
    int i;
    for(i=0;i<m;i++)
      printf("  x(%d)=%.2f\n",i,rhs[i]);
  }

  switch(args->solver) {
    case MUMPS:
      solver_finalize_dmumps(mumps_p);
      break;
    case UMFPACK:
      solver_finalize_umfpack(umfpack_p);
      break;
    default: assert(false); // should have caught unknown solver before now!
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
  if(extra_timing && args->rep == 0) {
    perftimer_inc(timer,"clean up",-1);
    perftimer_adjust_depth(timer,+1);

    perftimer_inc(timer,"MPI",-1);
  }
  ierr = MPI_Finalize(); assert(ierr == 0);

  // show timing info, if requested, to depth N
  if(extra_timing && args->rep == 0) {
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
