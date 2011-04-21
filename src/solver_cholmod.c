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
#include "solver_cholmod.h"
#include "solvers.h"
#include "matrix.h"
#include <cholmod.h>

#include <stdlib.h> // malloc, free
#include <string.h> // memcpy
#include <assert.h>

typedef struct {
  cholmod_common common;
  cholmod_factor* factor;
} solve_system_cholmod_t;

void solver_init_cholmod( solver_state_t* s ) {
  assert( s != NULL );

  solve_system_cholmod_t* const p = calloc( 1, sizeof( solve_system_cholmod_t ) );
  s->specific = p;
  assert( s->specific != NULL );

  assert( cholmod_start( &p->common ) == 1 );
  // TODO handle errors nicely (avoid asserting on malloc failures...)
}


static inline void _matrix2cholmodsparse( matrix_t* A, cholmod_sparse* B );
static inline void _matrix2cholmodsparse( matrix_t* A, cholmod_sparse* B )  {
  // check we're starting with a sane matrix
  assert( validate_matrix( A ) == 0 );

  assert( A->sym == SM_SYMMETRIC );
  assert( A->location == UPPER_TRIANGULAR );
  assert( A->format == SM_CSC );
  assert( A->base == FIRST_INDEX_ZERO );

  // Note: must have symmetric and positive definite
  // TODO check its positive definite

  B->nrow = A->m;
  B->ncol = A->n;
  B->nzmax = A->nz;
  B->p = A->jj; // column pointers
  B->i = A->ii; // row indices
  B->nz = NULL; // we use packed matrices
  B->x = A->dd; // data
  // z is NULL unless complex and in MATLAB format
  // (complex is split into real and imag components)
  B->z = NULL;

  if ( A->location == MC_STORE_BOTH )
    B->stype = 0; // TODO optimization: could be +1 or -1, then other half would be *ignored*
  else if ( A->location == UPPER_TRIANGULAR )
    B->stype = + 1;
  else // LOWER_TRIANGULAR
    B->stype = -1;

  B->itype = CHOLMOD_INT; // TODO or CHOLMOD_INTLONG or CHOLMOD_LONG

  // TODO refactor
  switch ( A->data_type ) {
      // TODO handle CHOLMOD_ZOMPLEX (matlab split real/imag vectors)
    case REAL_DOUBLE:
      B->xtype = CHOLMOD_REAL;
      B->dtype = CHOLMOD_DOUBLE;
      break;
    case REAL_SINGLE:
      B->xtype = CHOLMOD_REAL;
      B->dtype = CHOLMOD_SINGLE;
      break;
    case COMPLEX_DOUBLE:
      B->xtype = CHOLMOD_COMPLEX;
      B->dtype = CHOLMOD_DOUBLE;
      break;
    case COMPLEX_SINGLE:
      B->xtype = CHOLMOD_COMPLEX;
      B->dtype = CHOLMOD_SINGLE;
      break;
    case SM_PATTERN:
      B->xtype = CHOLMOD_PATTERN;
      B->dtype = 0; // don't care
      break;
  }

  B->sorted = 0; // TODO TRUE if columns are sorted, FALSE otherwise
  B->packed = 1; // TODO nz is ignored
}
static inline void _cholmodsparse2matrix( cholmod_sparse* B, matrix_t* A );
static inline void _cholmodsparse2matrix( cholmod_sparse* B, matrix_t* A )  {

  clear_matrix( A );
  A->format = SM_CSC;
  A->base = FIRST_INDEX_ZERO;

  // TODO check there aren't any undefined fields in A?

  A->m = B->nrow;
  A->n = B->ncol;
  A->nz = (( int* ) B->p )[B->ncol];
  // TODO B->nzmax;

  A->jj = B->p; // column pointers
  A->ii = B->i; // row indices
  assert( B->nz == NULL ); // we use packed matrices
  A->dd = B->x; // data
  // z is NULL unless complex and in MATLAB format
  // (complex is split into real and imag components)
  assert( B->z == NULL );

  // Note: can't handle other types of symmetry
  A->sym = SM_SYMMETRIC;
  if ( B->stype == 0 )
    A->location = MC_STORE_BOTH;
  else if ( B->stype > 0 )
    A->location = UPPER_TRIANGULAR;
  else  // B->stype < 0
    A->location = LOWER_TRIANGULAR;


  assert( B->itype == CHOLMOD_INT ); // TODO or CHOLMOD_INTLONG or CHOLMOD_LONG

  // TODO refactor
  if ( B->xtype == CHOLMOD_REAL ) {
    if ( B->dtype == CHOLMOD_DOUBLE )
      A->data_type = REAL_DOUBLE;
    else if ( B->dtype == CHOLMOD_SINGLE )
      A->data_type = REAL_SINGLE;
    else
      assert( 0 ); // unknown type
  }
  else if ( B->xtype == CHOLMOD_COMPLEX ) {
    if ( B->dtype == CHOLMOD_DOUBLE )
      A->data_type = COMPLEX_DOUBLE;
    else if ( B->dtype == CHOLMOD_SINGLE )
      A->data_type = COMPLEX_SINGLE;
    else
      assert( 0 ); // unknown type
  }
  else if ( B->xtype == CHOLMOD_PATTERN ) {
    A->data_type = SM_PATTERN;
  }
  else {
    assert( 0 ); // can't handle this type
  }

  // TODO currently, don't care re: B->sorted == FALSE/TRUE
  assert( B->packed == 1 ); // TODO nz is ignored

  assert( validate_matrix( A ) == 0 );
}

static inline void _matrix2cholmoddense( matrix_t* A, cholmod_dense* B );
static inline void _matrix2cholmoddense( matrix_t* A, cholmod_dense* B ) {
  // and we have a valid 'x' and 'b'
  int ierr = convert_matrix( A, DCOL, FIRST_INDEX_ZERO );
  assert( ierr == 0 );

  B->nrow = A->m;
  B->ncol = A->n;
  B->nzmax = A->m * A->n;
  B->d = B->nrow;
  B->x = A->dd;
  B->z = NULL; // TODO zomplex (matlab split real/imag)

  switch ( A->data_type ) {
      // TODO handle CHOLMOD_ZOMPLEX (matlab split real/imag vectors)
    case REAL_DOUBLE:
      B->xtype = CHOLMOD_REAL;
      B->dtype = CHOLMOD_DOUBLE;
      break;
    case REAL_SINGLE:
      B->xtype = CHOLMOD_REAL;
      B->dtype = CHOLMOD_SINGLE;
      break;
    case COMPLEX_DOUBLE:
      B->xtype = CHOLMOD_COMPLEX;
      B->dtype = CHOLMOD_DOUBLE;
      break;
    case COMPLEX_SINGLE:
      B->xtype = CHOLMOD_COMPLEX;
      B->dtype = CHOLMOD_SINGLE;
      break;
    case SM_PATTERN:
      B->xtype = CHOLMOD_PATTERN;
      B->dtype = 0; // don't care
      break;
  }
}

static inline void _cholmoddense2matrix( cholmod_dense* B, matrix_t* A );
static inline void _cholmoddense2matrix( cholmod_dense* B, matrix_t* A ) {
  clear_matrix( A );
  A->format = DCOL;

  A->m = B->nrow;
  A->n = B->ncol;
  A->nz = A->m * A->n;
  // TODO use B->nzmax?
  // TODO? B->d = B->nrow;
  A->dd = B->x;
  assert( B->z == NULL ); // TODO zomplex (matlab split real/imag)

  // TODO refactor
  if ( B->xtype == CHOLMOD_REAL ) {
    if ( B->dtype == CHOLMOD_DOUBLE )
      A->data_type = REAL_DOUBLE;
    else if ( B->dtype == CHOLMOD_SINGLE )
      A->data_type = REAL_SINGLE;
    else
      assert( 0 ); // unknown type
  }
  else if ( B->xtype == CHOLMOD_COMPLEX ) {
    if ( B->dtype == CHOLMOD_DOUBLE )
      A->data_type = COMPLEX_DOUBLE;
    else if ( B->dtype == CHOLMOD_SINGLE )
      A->data_type = COMPLEX_SINGLE;
    else
      assert( 0 ); // unknown type
  }
  else if ( B->xtype == CHOLMOD_PATTERN ) {
    A->data_type = SM_PATTERN;
  }
  else {
    assert( 0 ); // can't handle this type
  }
}


void solver_analyze_cholmod( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  assert( A != NULL );
  assert( A->sym == SM_SYMMETRIC );
  solve_system_cholmod_t* const p = s->specific;
  assert( p != NULL );

  cholmod_sparse B;
  _matrix2cholmodsparse( A, &B );

  p->factor = cholmod_analyze( &B, &p->common );
  assert( p->factor != NULL );
  // TODO handle errors nicely
}

void solver_factorize_cholmod( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  assert( A != NULL );
  assert( A->sym == SM_SYMMETRIC );
  solve_system_cholmod_t* const p = s->specific;
  assert( p != NULL );

  cholmod_sparse B;
  _matrix2cholmodsparse( A, &B );

  assert( cholmod_factorize( &B, p->factor, &p->common ) == 1 );
  // TODO handle errors nicely
}

// b can be sparse...
static void _solver_evaluate_cholmod_dense( solver_state_t* s, matrix_t* b, matrix_t* x );
static void _solver_evaluate_cholmod_sparse( solver_state_t* s, matrix_t* b, matrix_t* x );

// TODO can b be a matrix (vs a vector)?
void solver_evaluate_cholmod( solver_state_t* s, matrix_t* b, matrix_t* x ) {
  assert( s != NULL );
  assert( b != NULL );
  assert( b != x ); // TODO allow this form
  convert_matrix( b, DCOL, 0 ); // TODO rm (see below)
  if ( b->format == DCOL || b->format == DROW )
    _solver_evaluate_cholmod_dense( s, b, x );
  else
    _solver_evaluate_cholmod_sparse( s, b, x );
//
//    TODO sparse solver requires _matrix2cholmodsparse be able to handle
//    unsymmetric matrices properly, while checking for symmetric positive
//    definite A elsewhere
}

static void _solver_evaluate_cholmod_dense( solver_state_t* s, matrix_t* b, matrix_t* x ) {
  solve_system_cholmod_t* const p = s->specific;
  assert( p != NULL );

  // solve
  cholmod_dense B;
  _matrix2cholmoddense( b, &B );
  cholmod_dense* xc = cholmod_solve( CHOLMOD_A, p->factor, &B, &p->common );
  assert( xc != NULL );

  // xc -> x
  _cholmoddense2matrix( xc, x );
  // TODO assert( cholmod_free_dense( &xc, &p->common ) == 1 );
  // TODO assert( xc == NULL );
}

static void _solver_evaluate_cholmod_sparse( solver_state_t* s, matrix_t* b, matrix_t* x ) {
  solve_system_cholmod_t* const p = s->specific;
  assert( p != NULL );

  // solve
  cholmod_sparse B;
  _matrix2cholmodsparse( b, &B );
  cholmod_sparse* xc = cholmod_spsolve( CHOLMOD_A, p->factor, &B, &p->common );
  assert( xc != NULL );

  // TODO xc -> x
  _cholmodsparse2matrix( xc, x );
  // TODO assert( cholmod_free_sparse( &xc, &p->common ) == 1 );
  // TODO assert( xc == NULL );
}

void solver_finalize_cholmod( solver_state_t* s ) {
  if ( s == NULL )
    return;

  solve_system_cholmod_t* const p = s->specific;

  // release memory
  if ( p != NULL ) {
    if ( p->factor != NULL )
      cholmod_free_factor( &p->factor, &p->common );
    cholmod_finish( &p->common );
  }
  free( p );

}
