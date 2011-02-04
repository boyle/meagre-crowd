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
#include "solver_paradiso.h"
#include "solvers.h"
#include "matrix.h"
//#include <paradiso.h>

#include <stdlib.h> // malloc, free
#include <string.h> // memcpy
#include <assert.h>

// NOTE: paradiso doesn't include a header file for its library, use these prototypes instead
void pardisoinit (void   * PT, int    * MTYPE,   int * SOLVER, int * IPARM, double * DPARM, int *ERROR);
void pardiso     (void   * PT, int    * MAXFCT,   int * MNUM, int * MTYPE,    int * PHASE, int * N,
                  double * A, int    * IA,    int * JA, int * PERM,   int * NRHS, int * IPARM,
                     int * MSGLVL, double * B, double * X, int * ERROR, double * DPARM);
void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
void pardiso_chkvec     (int *, int *, double *, int *);
void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);

enum paradiso_mtype { REAL_PATTERN_SYM=1, REAL_SYM_POSDEF = 2, REAL_SYM_INDEF = -2,
                      COMPLEX_PATTERN_SYM = 3, COMPLEX_HERMITIAN_POSDEF = 4,
                      COMPLEX_HERMITIAN_INDEF = -4,
                      COMPLEX_SYM = 6,
                      REAL_UNSYM = 11,
                      COMPLEX_UNSYM = 13 };
enum paradiso_solver { SPARSE_DIRECT_SOLVER = 0, MULTIRECURSIVE_ITERATIVE_SOLVER = 1};

enum paradiso_solve_error { NO_ERROR=0, INCONSISTENT_INPUT=-1,
                            MEM=-2, REORDERING=-3,
                            ZERO_PIVOT_NUMERICAL_FACT_OR_ITR_REFINEMENT=-4,
                            INTERNAL_ERR=-5, PREORDER_FAILED=-6, DIAG_MATRIX_PROB=-7,
                            INT32_OVERFLOW=-8,
                            // paradiso.lic problems
                            NO_LICENSE=-10, EXPIRED_LICENSE=-11, LICENSE_WRONG_HOST_USER=-12,
                            // Krylov subspace issues
                            MAX_ITER=-100, NO_CONVERGENCE=-101, // (w/in 25 iterations)
                            ITER_ERR=-102, ITER_BREAKDOWN=-103 };

#define IPARM(I) iparm[(I)-1]
#define DPARM(I) dparm[(I)-1]
typedef struct {
  int Arows;
  int Acols;
  int* Aii; // these are pointers to existing data (DON'T free)
  int* Ajj;
  double* Add;

  // Paradiso config
  int *PT; // int *64 (Paradiso private ptr)
  int MTYPE;
  int SOLVER; // 0 sparse direct, 1 multi-recursive iterative
  int *iparm; // int * 64 - config
  double *dparm; // double * 64 - config
  int error;
} solve_system_paradiso_t;

void solver_init_paradiso( solver_state_t* s ) {
  assert( s != NULL );
  solve_system_paradiso_t * const p = calloc( 1, sizeof( solve_system_paradiso_t ) );
  assert( p != NULL );
  s->specific = p;

  p->PT = malloc(64*sizeof(int));
  assert(p->PT != NULL);
  p->MTYPE = 11; // default to real/unsym
  p->SOLVER = 0; // default to direct solver
  p->iparm = calloc(64,sizeof(int));
  p->dparm = calloc(64,sizeof(double));

  p->IPARM(3) = 1; // must be set to match number of processors TODO c_omp * c_mpi or just c_omp???
  p->IPARM(52) = 1; // number of compute nodes (MPI)

  // launch paradiso
  // TODO do we need to know the matrix type when we start here?? might need to move to analyze stage...
  pardisoinit(p->PT, &(p->MTYPE), &(p->SOLVER), p->iparm, p->dparm, &(p->error) );
  if(p->error != NO_ERROR)
    fprintf(stderr, "error: paradiso initialization code %d\n", p->error); // TODO decode
  assert(p->error == NO_ERROR);
}

void solver_analyze_paradiso( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  assert( A != NULL );
  solve_system_paradiso_t* const p = s->specific;
  assert( p != NULL );

  // TODO could do this smarter by giving MUMPS the symmetric matrix to solve,
  // instead we're just going straight for the unsymmetric solver
  // TODO this logic should really be moved to the wrapper solver() functions using the capabilities masks
  if ( A->sym == SM_SYMMETRIC )
    convert_matrix_symmetry( A, BOTH );

  // prepare the matrix
  int ierr = convert_matrix( A, SM_CSC, FIRST_INDEX_ZERO );
  assert( ierr == 0 );
  assert(( A->sym == SM_UNSYMMETRIC ) || ( A->sym == SM_SYMMETRIC ) );
  assert( A->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO
  assert( A->m == A->n ); // TODO can only handle square matrices at present (UMFPACK?)

  // Compressed Column Format
  assert( A->jj[0] == 0 );
  assert( A->jj[A->n] == A->nz );

  // umfpack_di_symbolic( A->m, A->n, ( int* ) A->jj, ( int* ) A->ii, A->dd, &( p->Symbolic ), NULL, NULL );
}

void solver_factorize_paradiso( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  assert( A != NULL );
  solve_system_paradiso_t* const p = s->specific;
  assert( p != NULL );

  // prepare the matrix
  int ierr = convert_matrix( A, SM_CSC, FIRST_INDEX_ZERO );
  assert( ierr == 0 );
  assert(( A->sym == SM_UNSYMMETRIC ) || ( A->sym == SM_SYMMETRIC ) );
  assert( A->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO
  assert( A->m == A->n ); // TODO can only handle square matrices at present (UMFPACK?)

  // saved for evaluation phase
  p->Arows = A->m;
  p->Acols = A->n;
  p->Ajj = ( int* ) A->jj;
  p->Aii = ( int* ) A->ii;
  p->Add = A->dd;
  // umfpack_di_numeric( p->Ajj, p->Aii, p->Add, p->Symbolic, &( p->Numeric ), NULL, NULL );
}

// TODO b can be sparse... ??
// TODO can b be a matrix (vs a vector)?
void solver_evaluate_paradiso( solver_state_t* s, matrix_t* b, matrix_t* x ) {
  assert( s != NULL );
  assert( b != NULL );
  assert( b != x ); // TODO allow this form
  solve_system_paradiso_t* const p = s->specific;
  assert( p != NULL );

  // and we have a valid 'x' and 'b'
  int ierr = convert_matrix( b, DCOL, FIRST_INDEX_ZERO );
  assert( ierr == 0 );
  assert( b->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO
  assert( b->n == 1 );
  assert( b->m == p->Acols ); // TODO move to wrapper level check?
  assert( b->data_type == REAL_DOUBLE );

  //umfpack_di_solve( UMFPACK_A, p->Ajj, p->Aii, p->Add, x->dd, b->dd, p->Numeric, NULL, NULL ) ;
}

void solver_finalize_paradiso( solver_state_t* s ) {
  if ( s == NULL )
    return;

  solve_system_paradiso_t* const p = s->specific;

  // release memory
  if ( p != NULL ) {
    free(p->PT);
    free(p->iparm);
    free(p->dparm);
  }
  free( p );
}
