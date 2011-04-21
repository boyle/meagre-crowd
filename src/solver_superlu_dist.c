/* Meagre-Crowd: A sparse distributed matrix solver testbench for performance benchmarking.
 * Copyright (C) 2011 Alistair Boyle <alistair.js.boyle@gmail.com>
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
#include "solver_superlu_dist.h"
#include "solvers.h"
#include "matrix.h"
//#include <superlu_dist.h>

#include <stdlib.h> // malloc, free
#include <string.h> // memcpy
#include <assert.h>
#include <stdint.h> // int64_t
#include <unistd.h> // dup -> redirect stdout temporarily

#include <mpi.h>

typedef struct {
  int Arows;
  int Acols;
  int* Aii; // these are pointers to existing data (DON'T free)
  int* Ajj;
  double* Add;

  // Paradiso config
  int *iparm; // int * 64 - config
  double *dparm; // double * 64 - config
} solve_system_superlu_dist_t;

void solver_init_superlu_dist( solver_state_t* s ) {
  assert( s != NULL );
  solve_system_superlu_dist_t * const p = calloc( 1, sizeof( solve_system_superlu_dist_t ) );
  assert( p != NULL );
  s->specific = p;

  p->iparm = calloc( 64, sizeof( int ) );
  p->dparm = calloc( 64, sizeof( double ) );

  // TODO init
}

// TODO split analyze stage into ordering and symbolic factorization stages?
void solver_analyze_superlu_dist( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  solve_system_superlu_dist_t* const p = s->specific;
  assert( p != NULL );
//  p->IPARM(2) = WSMP_ANALYZE;
//  p->IPARM(3) = WSMP_ANALYZE;
/*
  int TEST = 1; // Note: WSMP needs this value or it dies trying to read a NULL value
  if(s->mpi_rank == 0) {
    assert( A != NULL );

    // TODO could do this smarter by giving WSMP the symmetric matrix to solve,
    // instead we're just going straight for the unsymmetric solver
    assert( A->sym == SM_UNSYMMETRIC );
    // TODO for SYMMETRIC matrices: superlu_dist wants all diagonal entries to exist in A, even if they are zero

    // prepare the matrix
    assert( A->format == SM_CSR );
    assert( A->base == FIRST_INDEX_ONE );
    assert( A->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO
    assert( A->m == A->n ); // TODO can only handle square matrices at present???

    // Compressed Row Format
    assert( A->ii[0] == 1 );
    assert( A->ii[A->n] == A->nz + 1 );

    int N = A->m; // rows of A
//    pwgsmp_( &N, ( int* ) A->ii, ( int* ) A->jj, A->dd, NULL, &TEST, NULL, NULL, p->iparm, p->dparm);
  }
  else {
    assert( A == NULL );
    int N = 0;
//    pwgsmp_( &N, NULL, NULL, NULL, NULL, &TEST, NULL, NULL, p->iparm, p->dparm);
  }
  const int error_code = p->IPARM(64);
  if ( error_code != NO_ERROR )
    fprintf( stderr, "error: superlu_dist analyze code %d\n", error_code ); // TODO decode
  assert( error_code == NO_ERROR );
*/
}

void solver_factorize_superlu_dist( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  solve_system_superlu_dist_t* const p = s->specific;
  assert( p != NULL );
//  p->IPARM(2) = WSMP_FACTORIZE;
//  p->IPARM(3) = WSMP_FACTORIZE;
/*
  int TEST = 1;
  if(s->mpi_rank == 0) {
    // TODO symmetric matrix handling
    assert( A != NULL );
    assert( A->sym == SM_UNSYMMETRIC );

    // prepare the matrix
    assert( A->format == SM_CSR );
    assert( A->base == FIRST_INDEX_ONE );
    assert( A->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO
    assert( A->m == A->n ); // TODO can only handle square matrices at present (UMFPACK?)


    // save A for the solve step.. needed for iterative refinement // TODO add A to solve stage to allow iterative refinement!, remove this (and from other solvers)
    p->Arows = A->m;
    p->Acols = A->n;
    p->Ajj = ( int* ) A->jj;
    p->Aii = ( int* ) A->ii;
    p->Add = A->dd;

    int N = A->m; // rows of A
//    pwgsmp_( &N, ( int* ) A->ii, ( int* ) A->jj, A->dd, NULL, &TEST, NULL, NULL, p->iparm, p->dparm);
  }
  else {
    int N = 0; // rows of A
//    pwgsmp_( &N, NULL, NULL, NULL, NULL, &TEST, NULL, NULL, p->iparm, p->dparm);
  }
  const int error_code = p->IPARM(64);
  if ( error_code != NO_ERROR )
    fprintf( stderr, "error: superlu_dist factorize code %d\n", error_code ); // TODO decode
  assert( error_code == NO_ERROR );
*/
}

// TODO b can be sparse... ??
// TODO can b be a matrix (vs a vector)?
void solver_evaluate_superlu_dist( solver_state_t* s, matrix_t* b, matrix_t* x ) {
  assert( s != NULL );
  solve_system_superlu_dist_t* const p = s->specific;
  assert( p != NULL );
//  p->IPARM(2) = WSMP_EVALUATE;
//  p->IPARM(3) = WSMP_ITERATIVE_REFINEMENT; // TODO split into seperate stage
/*
  if(s->mpi_rank == 0) {
    assert( b != NULL );
    assert( b != x ); // TODO allow this form

    // and we have a valid 'x' and 'b'
    assert( b->format == DCOL );
    assert( b->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO

    // allocate data space // TODO if required
    // TODO move this to master function
    if(b == x) {
      // go ahead, we're expecting b to be destroyed
      // need to put data from b into x, clear x's ptr
      clear_matrix( x );
      x->format = DCOL;
      x->sym = SM_UNSYMMETRIC;
      x->data_type = b->data_type;
      x->m = p->Arows;
      x->n = b->n;
      x->nz = x->m * x->n;
      if(b->nz > x->nz) {
	x->dd = malloc(b->nz * sizeof( double ));
      }
      else {
	x->dd = malloc(x->nz * sizeof( double ));
      }
      x->dd = b->dd;
      b->dd = NULL;
      clear_matrix(b);
    }
    else {
      // need to copy b since it's destroyed in the process
      clear_matrix(x);
      matrix_t* t = copy_matrix(b);
      // push copied t contents into x
      *x = *t;
      t->dd = NULL; // transfer dd pointer ownership to x
      x->m = p->Arows;
      x->n = b->n;
      x->nz = x->m * x->n;
      free(t);
      if(x->nz > b->nz) {
	x->dd = realloc(x->dd, x->nz * sizeof( double ));
	assert(t != NULL); // realloc failure -- TODO proper error code
      }
      assert(x->dd != NULL);
    }

    int N = p->Arows; // rows of A
    double* B = x->dd; // N x NRHS
    int LDB = b->m;  // rows of B, must be >= N
    int NRHS = b->n; // columns of B

    // have to share number of rhs with other workers
    MPI_Bcast(&NRHS, 1, MPI_INT, 0, MPI_COMM_WORLD); // TODO pass in MPI communicator

    assert(p->Aii != NULL);
    assert(p->Ajj != NULL);
    assert(p->Add != NULL);
    assert(B != NULL);
//    pwgsmp_( &N, ( int* ) p->Aii, ( int* ) p->Ajj, p->Add, B, &LDB, &NRHS, NULL, p->iparm, p->dparm);
  }
  else {
    int N = 0;
    int LDB = 1;
    int NRHS; // need to get NRHS from master
    MPI_Bcast(&NRHS, 1, MPI_INT, 0, MPI_COMM_WORLD); // TODO pass in MPI communicator
//    pwgsmp_( &N, NULL, NULL, NULL, NULL, &LDB, &NRHS, NULL, p->iparm, p->dparm);
  }
  const int error_code = p->IPARM(64);
  if ( error_code != NO_ERROR )
    fprintf( stderr, "error: superlu_dist evaluate code %d\n", error_code ); // TODO decode
  assert( error_code == NO_ERROR );
*/
}

void solver_finalize_superlu_dist( solver_state_t* s ) {
  if ( s == NULL )
    return;

  solve_system_superlu_dist_t* const p = s->specific;

  // release memory
  if ( p != NULL ) {
    int error_code = 0;

    // for symmetric matrices its wssmp/pwssmp
    // superlu_dist_clear(); // for SMP or
//    psuperlu_dist_clear(); // MPI
    if ( error_code != 0 )
      fprintf( stderr, "error: superlu_dist finalize code %d\n", error_code ); // TODO decode
    assert( error_code == 0 );

    free( p->iparm );
    free( p->dparm );
  }
  free( p );
  s->specific = NULL;
}
