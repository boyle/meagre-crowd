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

#include "config.h"
#include "solvers.h"
#include "perftimer.h"
#include "matrix.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include <mpi.h>
#include <stdlib.h>

#include "solver_lookup.h"

static inline int _convert_matrix_A( const int solver, matrix_t* A );
static inline int _convert_matrix_A( const int solver, matrix_t* A ) {
  const unsigned int c = solver_lookup[solver].capabilities;
  // symmetry: don't convert unless we have to
  switch(A->sym) {
    case SM_UNSYMMETRIC:
      assert(c & SOLVES_UNSYMMETRIC); // nothing we can do...
      break;
    case SM_SYMMETRIC:
      if(c & SOLVES_SYMMETRIC) {
        switch(A->location) {
	  case MC_STORE_BOTH:
	    if(!(c & SOLVES_SYMMETRIC_BOTH)) {
	      if(c & SOLVES_SYMMETRIC_UPPER_TRIANGULAR) {
	        int ierr = convert_matrix_symmetry( A, UPPER_TRIANGULAR );
		assert(ierr == 0);
	      }
	      else if(c & SOLVES_SYMMETRIC_LOWER_TRIANGULAR) {
	        int ierr = convert_matrix_symmetry( A, LOWER_TRIANGULAR );
		assert(ierr == 0);
	      }
	      else {
	        assert(0); // unhittable? TODO
	      }
	    }
	    break;
	  case UPPER_TRIANGULAR:
	    if(!(c & SOLVES_SYMMETRIC_UPPER_TRIANGULAR)) {
	      // choose lower triangular first: a cheaper conversion and faster calculations later
	      if(c & SOLVES_SYMMETRIC_LOWER_TRIANGULAR) {
	        int ierr = convert_matrix_symmetry( A, LOWER_TRIANGULAR );
		assert(ierr == 0);
	      }
	      else if(c & SOLVES_SYMMETRIC_BOTH) {
	        int ierr = convert_matrix_symmetry( A, MC_STORE_BOTH );
		assert(ierr == 0);
	      }
	      else {
	        assert(0); // unhittable? TODO
	      }
	    }
	    break;
	  case LOWER_TRIANGULAR:
	    if(!(c & SOLVES_SYMMETRIC_LOWER_TRIANGULAR)) {
	      // choose upper triangular first: a cheaper conversion and faster calculations later
	      if(c & SOLVES_SYMMETRIC_UPPER_TRIANGULAR) {
	        int ierr = convert_matrix_symmetry( A, UPPER_TRIANGULAR );
		assert(ierr == 0);
	      }
	      else if(c & SOLVES_SYMMETRIC_BOTH) {
	        int ierr = convert_matrix_symmetry( A, MC_STORE_BOTH );
		assert(ierr == 0);
	      }
	      else {
	        assert(0); // unhittable? TODO
	      }
	    }
	    break;
	}
      }
      else if(c & SOLVES_UNSYMMETRIC) {
        // so we can convert it to an unsymmetric matrix
	int ierr = convert_matrix_symmetry( A, MC_STORE_BOTH );
	assert(ierr == 0);
	A->sym = SM_UNSYMMETRIC; // forget we were symmetric
      }
      else {
        assert(0); // TODO
      }
      break;
    case SM_SKEW_SYMMETRIC:
      assert(0); // TODO
      break;
    case SM_HERMITIAN:
      assert(0); // TODO
      break;
  }

  // format: don't convert unless we have to
  enum matrix_base_t base = A->base;
  switch(base) {
    case FIRST_INDEX_ZERO:
      if(!(c & SOLVES_BASE_ZERO)) {
        base = FIRST_INDEX_ONE; // need to convert
      }
      break;
    case FIRST_INDEX_ONE:
      if(!(c & SOLVES_BASE_ONE)) {
        base = FIRST_INDEX_ZERO; // need to convert
      }
      break;
  }

  // If you've given us a dense matrix for A, we need to convert it
  // to COO before going to some other format, so start there.
  if(A->format == DROW || A->format == DCOL) {
    int ierr = convert_matrix(A, SM_COO, base);
    assert(ierr == 0);
  }

  // TODO should we prefer CSC/CSR over COO (more efficient?)
  enum matrix_format_t format = A->format;
  switch(format) {
    case INVALID:
      assert(0); // fail!
      break;
    case DROW:
    case DCOL:
      assert(0); // TODO .. decide if we should be doing this in a sparse format?
    case SM_COO:
      if(!(c & SOLVES_FORMAT_COO)) { // can't handle COO format... figure out what to do next
        if(c & SOLVES_FORMAT_CSR) {
	  format = SM_CSR;
	}
	else if(c & SOLVES_FORMAT_CSC) {
	  format = SM_CSC;
	}
	else {
	  assert(0); // TODO
	}
      }
      break;
    case SM_CSC:
      if(!(c & SOLVES_FORMAT_CSC)) { // can't handle COO format... figure out what to do next
        if(c & SOLVES_FORMAT_COO) {
	  format = SM_COO;
	}
	else if(c & SOLVES_FORMAT_CSR) {
	  format = SM_CSR;
	}
	else {
	  assert(0); // TODO
	}
      }
      break;
    case SM_CSR:
      if(!(c & SOLVES_FORMAT_CSR)) { // can't handle COO format... figure out what to do next
        if(c & SOLVES_FORMAT_COO) {
	  format = SM_COO;
	}
	else if(c & SOLVES_FORMAT_CSC) {
	  format = SM_CSC;
	}
	else {
	  assert(0); // TODO
	}
      }
      break;
  }
  int ierr = convert_matrix( A, format, base );
  assert( ierr == 0 );
  return 0; // TODO return error code
}

static inline int _convert_matrix_b( const int solver, matrix_t* b );
static inline int _convert_matrix_b( const int solver, matrix_t* b ) {
  const unsigned int c = solver_lookup[solver].capabilities;
  // symmetry: don't convert unless we have to
  // Note: no solvers do expect to see a symmetric right-hand side, so convert to full matrix

  int ierr;
  if(b->sym != SM_UNSYMMETRIC) {
    ierr = convert_matrix_symmetry( b, MC_STORE_BOTH );
    assert(ierr == 0);
  }

  // format: don't convert unless we have to

  // should have at least one of these set in solver_lookup.h
  assert(c & (SOLVES_BASE_ZERO | SOLVES_BASE_ONE));
  enum matrix_base_t base = b->base;
  switch(base) {
    case FIRST_INDEX_ZERO:
      if(!(c & SOLVES_BASE_ZERO)) {
        base = FIRST_INDEX_ONE; // need to convert
      }
      break;
    case FIRST_INDEX_ONE:
      if(!(c & SOLVES_BASE_ONE)) {
        base = FIRST_INDEX_ZERO; // need to convert
      }
      break;
  }

  // should have at least one of these set in solver_lookup.h
  assert(c & (SOLVES_RHS_DROW | SOLVES_RHS_DCOL | SOLVES_RHS_COO | SOLVES_RHS_CSC | SOLVES_RHS_CSR));
  enum matrix_format_t format = b->format;
  switch(format) {
    case INVALID:
      assert(0); // fail!
      break;
    case DROW:
      if(!(c & SOLVES_RHS_DROW)) {
        if(c & SOLVES_RHS_DCOL) {
	  format = DCOL;
	}
	else if(c & SOLVES_RHS_COO) {
	  format = SM_COO;
	}
	else if(c & SOLVES_RHS_CSC) {
	  format = SM_CSC;
	}
	else if(c & SOLVES_RHS_CSR) {
	  format = SM_CSR;
	}
	else {
	  assert(0);
	}
      }
      break;
    case DCOL:
      if(!(c & SOLVES_RHS_DCOL)) {
        if(c & SOLVES_RHS_DROW) {
	  format = DROW;
	}
	else if(c & SOLVES_RHS_COO) {
	  format = SM_COO;
	}
	else if(c & SOLVES_RHS_CSC) {
	  format = SM_CSC;
	}
	else if(c & SOLVES_RHS_CSR) {
	  format = SM_CSR;
	}
	else {
	  assert(0);
	}
      }
      break;
    case SM_COO:
      if(!(c & SOLVES_RHS_COO)) {
	if(c & SOLVES_RHS_CSC) {
	  format = SM_CSC;
	}
	else if(c & SOLVES_RHS_CSR) {
	  format = SM_CSR;
	}
        else if(c & SOLVES_RHS_DROW) {
	  format = DROW;
	}
	else if(c & SOLVES_RHS_DCOL) {
	  format = DCOL;
	}
	else {
	  assert(0);
	}
      }
      break;
    case SM_CSC:
      if(!(c & SOLVES_RHS_CSC)) {
	if(c & SOLVES_RHS_COO) {
	  format = SM_COO;
	}
	else if(c & SOLVES_RHS_CSR) {
	  format = SM_CSR;
	}
        else if(c & SOLVES_RHS_DROW) {
	  format = DROW;
	}
	else if(c & SOLVES_RHS_DCOL) {
	  format = DCOL;
	}
	else {
	  assert(0);
	}
      }
      break;
    case SM_CSR:
      if(!(c & SOLVES_RHS_CSR)) {
	if(c & SOLVES_RHS_COO) {
	  format = SM_COO;
	}
	else if(c & SOLVES_RHS_CSC) {
	  format = SM_CSC;
	}
        else if(c & SOLVES_RHS_DROW) {
	  format = DROW;
	}
	else if(c & SOLVES_RHS_DCOL) {
	  format = DCOL;
	}
	else {
	  assert(0);
	}
      }
      break;
  }
  ierr = convert_matrix( b, format, base );
  assert( ierr == 0 );
  return 0; // TODO return error code
}

static inline int _valid_solver( const int solver );
static inline int _valid_solver( const int solver ) {
  return ( solver_lookup[solver].shortname != NULL );
}

// lookup functions
int lookup_solver_by_shortname( const char* shortname ) {
  assert( shortname != NULL );
  int i = 0;
  while ( solver_lookup[i].shortname != NULL ) {
    if ( strncmp( solver_lookup[i].shortname, shortname, SOLVER_SHORTNAME_MAX_LEN ) == 0 )
      return i;
    i++;
  }
  return -1;
}

const char* solver2str( const int solver ) {
  if ( _valid_solver( solver ) )
    return solver_lookup[solver].name;
  else
    return "<invalid>";
}

void printf_solvers( const unsigned int verbosity ) {
  printf( "Available solvers:\n" );
  int i, j;
  // determine max shortname length so everything can be aligned nicely
  size_t max_shortname_len = 0;

  // count solvers
  for ( i = 0; solver_lookup[i].shortname != NULL; i++ );
  size_t * const l = malloc( i * sizeof( size_t ) );

  for ( i = 0; solver_lookup[i].shortname != NULL; i++ ) {
    size_t ll = strlen( solver_lookup[i].shortname );
    l[i] = ll;
    if ( ll > max_shortname_len )
      max_shortname_len = ll;
  }
  for ( i = 0; solver_lookup[i].shortname != NULL; i++ ) {
    printf( "  %s", solver_lookup[i].shortname );
    for ( j = 0; j < max_shortname_len - l[i]; j++ )
      printf( " " );
    printf( "    %s %s (%s, %s)\n",
            solver_lookup[i].name,
            solver_lookup[i].version,
            solver_lookup[i].author,
            solver_lookup[i].license );
    if ( verbosity >= 1 ) { // organization, references
      for ( j = 0; j < max_shortname_len; j++ )
        printf( " " );
      printf( "      %s, %s\n",
              solver_lookup[i].organization,
              solver_lookup[i].url );
    }
    if ( verbosity >= 2 ) {
      printf( "    References:\n%s\n\n", solver_lookup[i].references );
      // TODO capabilities
    }
  }

  free( l );
}



// --------------------------------------------
// can the preferred solver solve this problem?
//   e.g. can the solver only handle Symmetric Postive Definite (SPD) matrices
// returns: 1 yes, 0 no
int solver_can_do( const int solver, matrix_t* A, matrix_t* b ) {
  return 1; // TODO something more clever, like a real answer
}

int solver_uses_mpi( const int solver ) {
  return (( solver_lookup[solver].multicore & SOLVER_CAN_USE_MPI ) != 0 );
}
int solver_requires_mpi( const int solver ) {
  return (( solver_lookup[solver].multicore & SOLVER_REQUIRES_MPI ) != 0 );
}
int solver_uses_omp( const int solver ) {
  return (( solver_lookup[solver].multicore & SOLVER_CAN_USE_OMP ) != 0 );
}
int solver_requires_omp( const int solver ) {
  return (( solver_lookup[solver].multicore & SOLVER_REQUIRES_OMP ) != 0 );
}

// select the most appropriate solver for this problem
//  - is it small and thus should be solved single-threaded (single processor)
//  - is it moderate and should be solved SMP (shared memory)
//  - is it huge and should be solved MPI (distributed memory)
int select_solver( matrix_t* A, matrix_t* b ) {
  return 0; // TODO something clever
}

// --------------------------------------------
int mc_mpi_omp_initialize(const int solver, int* is_mpi, int* is_omp) {
  *is_mpi = 0;
  *is_omp = 0;
 // if this solver is not thread/mpi safe, report an error if its been launched that way
  const int uses_mpi     = solver_uses_mpi( solver );
  const int requires_mpi = solver_requires_mpi( solver );
  const int uses_omp     = solver_uses_omp( solver );
  const int requires_omp = solver_requires_omp( solver );
  int c_mpi;
  {
    const char* mpi_world_size = getenv( "OMPI_COMM_WORLD_SIZE" );
    if ( mpi_world_size == NULL ) { // no env value configured
      c_mpi = 0;
      // printf( "no OMPI_COMM_WORLD_SIZE\n" );
    }
    else { // works #ifdef OPEN_MPI -- we're using openMPI rather than lam or mpich...
      int ret = sscanf( mpi_world_size, "%d", &c_mpi );
      assert( ret == 1 );
      // printf( "OMPI_COMM_WORLD_SIZE=%d\n", c_mpi );
    }
  }
  int c_omp;
  {
    const char* omp_threads = getenv( "OMP_NUM_THREADS" );
    if ( omp_threads == NULL ) { // no env value configured
      c_omp = 0;
      // printf( "no OMP_NUM_THREADS\n" );
    }
    else { // works #ifdef OPEN_MP -- we're using openMP
      int ret = sscanf( omp_threads, "%d", &c_omp );
      assert( ret == 1 );
      // printf( "OMP_NUM_THREADS=%d\n", c_omp );
    }
  }
  *is_mpi = ( requires_mpi || ( uses_mpi && ( c_mpi != 0 ) ) );
  *is_omp = ( requires_omp || ( uses_omp && ( c_omp != 0 ) ) );
  const int is_single_threaded = ( !*is_mpi && !*is_omp );

  // start up MPI, if we're using it
  // TODO FIXME safety check to see if MPI_Init has already been called...
  const int is_single_threaded_but_mpi = is_single_threaded && ( c_mpi > 1 );
  if ( *is_mpi || is_single_threaded_but_mpi ) {
    int ierr;
    int provided_threading;
    // argc and argv are null (MPI_Init doesn't use these anyways, right??)
    ierr = MPI_Init_thread( NULL, NULL, MPI_THREAD_MULTIPLE, &provided_threading );
    assert( ierr == 0 );
//    assert(provided_threading == MPI_THREAD_MULTIPLE);
    ierr = MPI_Comm_size( MPI_COMM_WORLD, &c_mpi );
    assert( ierr == 0 );
  }

  // wait till after firing up MPI, so we only put out an error on the rank=0 machine
  if ( is_single_threaded_but_mpi ) {
    int mpi_rank;
    int ierr = MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    assert( ierr == 0 );
    if ( mpi_rank == 0 )
      fprintf( stderr, "error: selected solver (%s) is single threaded but was launched with %d threads\n",
               solver2str( solver ), c_mpi ); // TODO report error through error code
    ierr = MPI_Finalize();
    assert( ierr == 0 );
    return 10;
  }

  return 0;
}

int mc_mpi_omp_finalize(const int solver) {
  // TODO refactor common code
  const int requires_mpi = solver_requires_mpi( solver );
  const int uses_mpi     = solver_uses_mpi( solver );

  int c_mpi;
  {
    const char* mpi_world_size = getenv( "OMPI_COMM_WORLD_SIZE" );
    if ( mpi_world_size == NULL ) { // no env value configured
      c_mpi = 0;
      // printf( "no OMPI_COMM_WORLD_SIZE\n" );
    }
    else { // works #ifdef OPEN_MPI -- we're using openMPI rather than lam or mpich...
      int ret = sscanf( mpi_world_size, "%d", &c_mpi );
      assert( ret == 1 );
      // printf( "OMPI_COMM_WORLD_SIZE=%d\n", c_mpi );
    }
  }

  const int is_mpi = ( requires_mpi || ( uses_mpi && ( c_mpi != 0 ) ) );
  if ( is_mpi ) {
    int ierr = MPI_Finalize();
    assert( ierr == 0 );
  }

  return 0; // TODO return errors rather than assert()
}




// --------------------------------------------
// wrapper function: solve 'A x = b' for 'x'
// calls initialize, analyze, factorize, evaluate, finalize
// returns x, the solution
void mc_solver( const int solver, const int verbosity, const int mpi_rank, matrix_t* A, matrix_t* b, matrix_t* x ) {
  solver_state_t* s = solver_init( solver, verbosity, mpi_rank, NULL );
  assert( s != NULL ); // malloc failure
  solver_analyze( s, A );
  solver_factorize( s, A );
  solver_evaluate( s, b, x );
  solver_finalize( s );
}

// wrapper function: solve 'A x = b' for 'x' w/o re-initializing solver
// calls analyze, factorize, evaluate
// must call initialize before and finalize after when all done
// returns x, the solution
void solver_solve( solver_state_t* s, matrix_t* A, matrix_t* b, matrix_t* x ) {
  solver_analyze( s, A );
  solver_factorize( s, A );
  solver_evaluate( s, b, x );
}

// --------------------------------------------
// initialize and finalize the solver state
solver_state_t* solver_init( const int solver, const int verbosity, const int mpi_rank, perftimer_t* timer ) {
  solver_state_t* s = malloc( sizeof( solver_state_t ) );
  assert( s != NULL );

  // configure state
  s->solver = solver;
  s->verbosity = verbosity;
  s->mpi_rank = mpi_rank;
  s->timer = timer;
  s->specific = NULL;
  if ( _valid_solver( solver ) && ( solver_lookup[solver].init != NULL ) )
    solver_lookup[solver].init( s );

  return s;
}

void solver_finalize( solver_state_t* s ) {
  assert( s != NULL );
  // clean up (and deallocate "s->specific" if required)
  const int solver = s->solver;
  if ( _valid_solver( solver ) && ( solver_lookup[solver].analyze != NULL ) )
    solver_lookup[solver].finalize( s );

  // make sure it won't get deallocated twice by mistake
  s->specific = NULL;

  free( s );
}

// evaluate the patterns in A, doesn't care about the actual values in the matrix (A->dd)
void solver_analyze( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  const int solver = s->solver;
  if(s->mpi_rank == 0) {
    assert( A != NULL );
    _convert_matrix_A(solver, A);
  }
  perftimer_inc( s->timer, "analyze", -1 );
  if ( _valid_solver( solver ) && ( solver_lookup[solver].analyze != NULL ) )
    solver_lookup[solver].analyze( s, A );
}


// factorize the matrix A, A must have the same pattern of non-zeros at that used in the solver_analyze stage
void solver_factorize( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  const int solver = s->solver;
  if(s->mpi_rank == 0) {
    assert( A != NULL );
    _convert_matrix_A(solver, A);
  }
  perftimer_inc( s->timer, "factorize", -1 );
  if ( _valid_solver( solver ) && ( solver_lookup[solver].factorize != NULL ) )
    solver_lookup[solver].factorize( s, A );
}


// solve the matrix 'A' for right-hand side 'b'
// returns 'x', the solution
void solver_evaluate( solver_state_t* s, matrix_t* b, matrix_t* x ) {
  assert( s != NULL );
  const int solver = s->solver;
  const unsigned int c = solver_lookup[solver].capabilities;

  // decide if we need to use MPI
  int is_mpi;
  int ierr = MPI_Initialized(&is_mpi);
  assert(ierr == 0);

  matrix_t bb = {0}; // Note that bb is only a shallow copy, it gets thrown away afterwards
  int loops;
  if(s->mpi_rank == 0) {
    assert( b != NULL );
    assert( x != NULL );

    if((b->n != 1) && (c & SOLVES_RHS_VECTOR_ONLY)) {
      // TODO FIXME: also need to handle non-DCOL format (sparse, etc)
      // Force to single column: We do one column at a time and
      // monkey with the data pointer. Note that if the solver's
      // wrapper does data conversion internally and reallocates
      // memory this method will break! For those solvers, we need to
      // support their desired format so they can't break this wrapper.
      // This special consideration is only necessary for
      // SOLVES_RHS_VECTOR_ONLY solvers.
      // ensure the RHS is in column major format
      assert(c & SOLVES_RHS_DCOL);
      convert_matrix(b, DCOL, FIRST_INDEX_ZERO);
      bb = *b; // shallow copy
      bb.n = 1; // pretend this right-hand side is only one column
      loops = b->n;
    }
    else { // don't need to do anything special
      _convert_matrix_b(s->solver, b);
      bb = *b; // shallow copy
      loops = 1;
    }
  }

  // share with all nodes, how many loops do we need to do when we can't handle more than a vector RHS
  if(is_mpi) {
    ierr = MPI_Bcast(&loops, 1, MPI_INT, 0, MPI_COMM_WORLD); // Bcast(var, n, datatype, root_rank, communicator)
    assert(ierr == 0);
    // TODO should set communicator in solver_initialize() instead of using COMM_WORLD by default
  }

  perftimer_inc( s->timer, "evaluate", -1 );
  if ( _valid_solver( solver ) && ( solver_lookup[solver].evaluate != NULL ) ) {
    // prepare to receive the results
    if(s->mpi_rank == 0) {
      clear_matrix(x);
      x->format = DCOL;
      assert(x->dd == NULL);
      assert(x->n == 0);
    }
    // Loops are for handling multiple RHS when the underlying solver can't.
    //
    // Note that this loop trickery will be utterly broken if the solver
    //   converts any of the data internally since we are holding shallow copies.
    int i;
    for(i=0; i<loops; i++) {
      matrix_t *const bb_p = (s->mpi_rank == 0)? &bb : NULL; // only rank-0 has a valid RHS matrix
      matrix_t xx = {0};
      matrix_t *const xx_p = (s->mpi_rank == 0)? &xx : NULL; // only rank-0 has a valid solution
      solver_lookup[solver].evaluate( s, bb_p, xx_p ); // TODO bb should be const const to be protected

      if(s->mpi_rank == 0) {
        // don't bother mucking with bb if its the last loop
        if(i != loops -1) {
	  // for the next loop, advance by a column
	  assert((bb.format == DCOL) || (bb.format == DROW));
	  assert(bb.n == 1);
	  convert_matrix(&bb, DCOL, FIRST_INDEX_ZERO);
	  // as a short-circuit conversion DROW->DCOL, this doesn't
	  // actually replace the memory pointers, so we're safe to
	  // continue with the monkey-ing with pointers (it converts
	  // to the correct index-ing and changes from DROW -> DCOL)

	  double* dd = bb.dd;
	  dd += bb.m; // advance by a column
	  bb.dd = dd;
        }

	// collect the results
	// TODO should probably abstract this into the matrix functions
	if(x->dd == NULL) {
	  // if this is the first row we can just steal the pointer
	  *x = xx;
	  xx.dd = NULL;
	}
	else { // otherwise we need to
	  x->m = xx.m; // rows
	  convert_matrix(&xx, DCOL, FIRST_INDEX_ZERO);

	  // resize to capture another column
	  x->dd = realloc(x->dd, (x->m * (x->n + xx.n))*sizeof(double));
	  assert(x->dd != NULL); // realloc failure
	  // and copy the data into the newly resized buffer
	  memcpy(x->dd + (x->m * x->n) * sizeof(double), xx.dd, x->m * xx.n * sizeof(double));
	  x->n += xx.n; // increment the column count
	}
	clear_matrix(&xx); // free any memory we accumulated from the solver
      }
    }
  }
  else {
    clear_matrix( x );
    x->format = INVALID;
    // TODO report an error
  }
  perftimer_inc( s->timer, "done", -1 );
}
