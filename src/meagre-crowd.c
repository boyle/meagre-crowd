/* Meagre-Crowd: A sparse distributed matrix solver testbench for performance benchmarking.
 * Copyright (C) 2010 Alistair Boyle <alistair.js.boyle@gmail.com>
 *
 *     This file is part of Meagre-Crowd.
 *
 *     Meagre-Crowd program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 *
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
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
#include <float.h> // machine epsilon: LDBL_EPSILON, DBL_EPSILON, FLT_EPSILON

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
// returns: 1=match w/in precision, 0=non-matching
// TODO refactor mv to matrix.h, operate on matrix_t objects
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
  // TODO should default to appropriate epsilon for solver, may need to be *2 or some larger value given numerical instability? should print out epsilon of solution in verbose mode
  args->expected_precision = 5e-15; // TODO was DBL_EPSILON=1.11e-16 but not stored with enough digits? // default to machine epsilon for 'double'
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

  // if this solver is not thread/mpi safe, report an error if its been launched that way
  int single_threaded = 0;
  int c_mpi;
  switch(args->solver) {
    case UMFPACK: single_threaded = 1; break;
    default: ; // otherwise do nothing
  }
  ierr = MPI_Comm_size(MPI_COMM_WORLD,&c_mpi); assert(ierr == 0);
  if(single_threaded && c_mpi > 1) {
    if(args->mpi_rank == 0)
      fprintf(stderr, "error: selected solver is single threaded, but launched with %d threads\n",c_mpi);
    ierr = MPI_Finalize(); assert(ierr == 0);
    return 10;
  }

  // Define the problem on the host
  double* b = NULL;
  double* expected = NULL;
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
    m = matrix_rows(A);
    // TODO load a rhs from the --input matrix file

    // allocate an sequentially numbered right-hand side of A.m rows
    // TODO warn if there is already a rhs loaded
    if(args->rhs != NULL) {
      struct sparse_matrix_t* b_loaded;
      if((retval = load_matrix(args->rhs, &b_loaded)) != 0) {
        return retval;
      }

      int ierr = sparse_matrix_convert(b_loaded, COO); assert(ierr == 0);
      struct coo_matrix_t* bcoo = b_loaded->repr;

      assert(bcoo->m == m); // rows must match
      assert(bcoo->n == 1); // TODO can only handle single column vector currently

      b=malloc((bcoo->m)*sizeof(double));
      assert(b!=NULL);

      // convert from COO vector to dense format vector
      size_t i;
      size_t j = 0;
      double* data = bcoo->val;
      for(i=0; i < m; i++) {
        assert( bcoo->II[j] >= i ); // or somehow the row indicies weren't ordered properly
        if( bcoo->II[j] == i ) { // are we at the next row yet?
	  b[i] = data[j]; // good, store the value
	  j++;
	}
	else {
	  b[i] = 0; // or we've got another zero
	}
      }

      // clean up
      //destroy_coo_matrix(bcoo);
    }
    else {
      b = malloc(m*sizeof(double));
      assert(b != NULL); // malloc failure
      { // initialize right-hand-side (b)
	int i;
	for(i=0;i<m;i++)
	  b[i] = i;
      }
    }
    if(args->expected != NULL) {
      // TODO refactor: this is a cut and paste of the loader for 'b'
      struct sparse_matrix_t* expected_loaded;
      if((retval = load_matrix(args->expected, &expected_loaded)) != 0) {
        return retval;
      }

      int ierr = sparse_matrix_convert(expected_loaded, COO); assert(ierr == 0);
      struct coo_matrix_t* ecoo = expected_loaded->repr;

      assert(ecoo->m == m); // rows must match
      assert(ecoo->n == 1); // TODO can only handle single column vector currently

      expected=malloc((ecoo->m)*sizeof(double));
      assert(expected!=NULL);

      // convert from COO vector to dense format vector
      size_t i;
      size_t j = 0;
      double* data = ecoo->val;
      for(i=0; i < m; i++) {
        assert( ecoo->II[j] >= i ); // or somehow the row indicies weren't ordered properly
        if( ecoo->II[j] == i ) { // are we at the next row yet?
          expected[i] = data[j]; // good, store the value
          j++;
        }
        else {
          expected[i] = 0; // or we've got another zero
        }
      }

      // clean up
      //destroy_coo_matrix(ecoo);
    }


    if(extra_timing && args->rep == 0) {
      perftimer_inc(timer,"solver",-1);
      perftimer_inc(timer,"rhs",-1);
    }

    // verbose output
    if(args->verbosity >= 1) {
      assert(A != NULL);
      int c_omp;
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

      if(args->verbosity >= 2) {
        int i;
        int n = Acoo->nnz;
        double* v = Acoo->val;
        for(i=0;i<n;i++) {
          printf("  A(%i,%i)=%.2f\n", Acoo->II[i], Acoo->JJ[i], v[i]);
        }
      }

      if(args->verbosity >= 2) { // show the rhs matrix
	int i;
	for(i=0;i<Acoo->m;i++)
	  printf("  b(%i)=%.2f\n", i, b[i]);
      }
    }

    if(extra_timing && args->rep == 0)
      perftimer_adjust_depth(timer,-1);
    //destroy_sparse_matrix (A); // TODO can't release it unless we're copying it...
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

  // test result?
  if((args->mpi_rank == 0) && (expected != NULL)){
    perftimer_inc(timer,"test",-1);
    if( results_match(expected, rhs, m, args->expected_precision) ) {
      retval = 0;
      printf("PASS\n");
    }
    else {
      retval = 100;
      printf("FAIL\n");
    }
  }

  if((args->mpi_rank == 0) && (args->output != NULL)) {
    if(strncmp(args->output,"-",2) == 0) {
      int i;
      for(i=0;i<m;i++)
        printf("  x(%d)=%.2f\n",i,rhs[i]);
    }
    else {
      // convert from vector to matrix format
      struct coo_matrix_t Bcoo = {0};
      struct sparse_matrix_t B = {0};
      B.format = COO;
      B.repr = &Bcoo;
      Bcoo.m = m;
      Bcoo.n = 1;
      Bcoo.nnz = m;
      Bcoo.val = rhs;
      Bcoo.II = malloc(m*sizeof(int));
      Bcoo.JJ = malloc(m*sizeof(int));
      int i;
      for(i=0;i<m;i++) {
	Bcoo.II[i] = i;
	Bcoo.JJ[i] = 0;
      }
      Bcoo.index_base = ZERO;
      Bcoo.symmetry_type = UNSYMMETRIC;
      Bcoo.value_type = REAL;
      Bcoo.ownership = USER_DEALLOCATES;

      // TODO test args->output to decide if its an acceptable filename before now
      // TODO really, we shouldn't care what the file name is, just the file format...
      if(save_matrix(&B, args->output) != 0)
	retval = 12; // TODO normalize error codes so they are defined in one place (defines?)

      free(Bcoo.II);
      free(Bcoo.JJ);
    }
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
