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
#include "solver_umfpack.h"
#include "solvers.h"
#include "matrix.h"
#include <umfpack.h>

#include <stdlib.h> // malloc, free
#include <string.h> // memcpy
#include <assert.h>

typedef struct {
  int Arows;
  int Acols;
  int* Aii; // these are pointers to existing data (DON'T free)
  int* Ajj;
  double* Add;

  void* Symbolic;
  void* Numeric;
} solve_system_umfpack_t;

void solver_init_umfpack( solver_state_t* s ) {
  assert( s != NULL );
  s->specific = calloc( 1, sizeof( solve_system_umfpack_t ) );
  assert( s->specific != NULL );
}

void solver_analyze_umfpack( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  assert( A != NULL );
  solve_system_umfpack_t* const p = s->specific;
  assert( p != NULL );


  // TODO could do this smarter by giving UMFPACK the symmetric matrix to solve,
  // instead we're just going straight for the unsymmetric solver
  assert( A->sym == SM_UNSYMMETRIC );

  // prepare the matrix
  assert( A->format == SM_CSC );
  assert( A->base == FIRST_INDEX_ZERO );
  assert( A->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO
  assert( A->m == A->n ); // TODO can only handle square matrices at present (UMFPACK?)

  // Compressed Column Format
  assert( A->jj[0] == 0 );
  assert( A->jj[A->n] == A->nz );

  if(s->verbosity >= 4) {
    double Control [UMFPACK_CONTROL];
    umfpack_di_defaults(Control); // set defaults
    Control[UMFPACK_PRL] = 6;
    umfpack_di_report_matrix( A->m, A->n, ( int* ) A->jj, ( int* ) A->ii, A->dd, 0, Control);
  }

  int status = umfpack_di_symbolic( A->m, A->n, ( int* ) A->jj, ( int* ) A->ii, A->dd, &( p->Symbolic ), NULL, NULL );
  if(status != UMFPACK_OK)
    umfpack_di_report_status(NULL, status);
  assert(status == UMFPACK_OK);
}

void solver_factorize_umfpack( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  assert( A != NULL );
  solve_system_umfpack_t* const p = s->specific;
  assert( p != NULL );

  // prepare the matrix
  assert( A->format == SM_CSC );
  assert( A->base == FIRST_INDEX_ZERO );
  assert( A->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO
  assert( A->m == A->n ); // TODO can only handle square matrices at present (UMFPACK?)

  // saved for evaluation phase
  p->Arows = A->m;
  p->Acols = A->n;
  p->Ajj = ( int* ) A->jj;
  p->Aii = ( int* ) A->ii;
  p->Add = A->dd;

  int status = umfpack_di_numeric( p->Ajj, p->Aii, p->Add, p->Symbolic, &( p->Numeric ), NULL, NULL );
  if(status != UMFPACK_OK)
    umfpack_di_report_status(NULL, status);
  assert(status == UMFPACK_OK);
}

// TODO b can be sparse... ??
// TODO can b be a matrix (vs a vector)?
void solver_evaluate_umfpack( solver_state_t* s, matrix_t* b, matrix_t* x ) {
  assert( s != NULL );
  assert( b != NULL );
  assert( b != x ); // TODO allow this form
  solve_system_umfpack_t* const p = s->specific;
  assert( p != NULL );

  // and we have a valid 'x' and 'b'
  int ierr = convert_matrix( b, DCOL, FIRST_INDEX_ZERO );
  assert( ierr == 0 );
  assert( b->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO
  assert( b->n == 1 );
  assert( b->m == p->Acols ); // TODO move to wrapper level check?
  assert( b->data_type == REAL_DOUBLE );

  // allocate x, if required
  if (( x->format != DCOL ) || ( x->m != p->Acols ) || ( x->n != b->n ) ) {
    clear_matrix( x );
    x->format = DCOL;
    x->data_type = b->data_type;
    x->m = p->Acols;
    x->n = b->n;
    x->nz = x->m * x->n;
    x->dd = calloc(( x->m ) * ( x->n ), sizeof( double ) ); // TODO this shouldn't need to be a calloc!
    assert(x->dd != NULL);
  }

  int status = umfpack_di_solve( UMFPACK_A, p->Ajj, p->Aii, p->Add, x->dd, b->dd, p->Numeric, NULL, NULL );
  if(status != UMFPACK_OK)
    umfpack_di_report_status(NULL, status);
  assert(status == UMFPACK_OK);
}

void solver_finalize_umfpack( solver_state_t* s ) {
  if ( s == NULL )
    return;

  solve_system_umfpack_t* const p = s->specific;

  // release memory
  if ( p != NULL ) {
//    umfpack_di_free_numeric(p->Numeric);
//    umfpack_di_free_symbolic(p->Symbolic);
  }
  free( p );
}
