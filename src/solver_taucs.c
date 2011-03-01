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
#include "solver_taucs.h"
#include "solvers.h"
#include "matrix.h"
#include <taucs.h>

#include <stdlib.h> // malloc, free
#include <string.h> // memcpy
#include <assert.h>
#include <stdint.h> // int64_t
#include <unistd.h> // dup -> redirect stdout temporarily

typedef struct {
  int Arows;
  int Acols;
  int* Aii; // these are pointers to existing data (DON'T free)
  int* Ajj;
  double* Add;

  void* F; // saved factorization for solve stage
  taucs_io_handle* LU;
  int* colperm;
} solve_system_taucs_t;

static const char* taucs_errorcode_to_string_(int err);
static const char* taucs_errorcode_to_string_(int err) {
  switch(err) {
    case TAUCS_SUCCESS:          return "no error: SUCCESS!";
    case TAUCS_ERROR:            return "general error";
    case TAUCS_ERROR_NOMEM:      return "memory allocation failure";
    case TAUCS_ERROR_BADARGS:    return "bad arguments";
    case TAUCS_ERROR_MAXDEPTH:   return "recursion limit (stack overflow)";
    case TAUCS_ERROR_INDEFINITE: return "expecting posdef matrix";
    default: return "unknown error code";
  }
}

void solver_init_taucs( solver_state_t* s ) {
  assert( s != NULL );
  solve_system_taucs_t * const p = calloc( 1, sizeof( solve_system_taucs_t ) );
  assert( p != NULL );
  s->specific = p;

  // if set to be verbose set this to stdout TODO
  if(s->verbosity >= 3)
    taucs_logfile("stdout");
}

void solver_analyze_taucs( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  assert( A != NULL );
  solve_system_taucs_t* const p = s->specific;
  assert( p != NULL );

  // TODO could do this smarter by giving MUMPS the symmetric matrix to solve,
  // instead we're just going straight for the unsymmetric solver
  // TODO this logic should really be moved to the wrapper solver() functions using the capabilities masks
  if ( A->sym == SM_SYMMETRIC )
    convert_matrix_symmetry( A, BOTH );
  // TODO for SYMMETRIC matrices: taucs wants all diagonal entries to exist in A, even if they are zero
  // TODO taucs can handle pattern symmetric matrices that have unsymmetric data

  // TODO speial handling if the matrix has a symmetric pattern but unsymmetric data .. copy, delete data, set to pattern and test

  // prepare the matrix
  int ierr = convert_matrix( A, SM_CSC, FIRST_INDEX_ZERO );
  assert( ierr == 0 );
  assert(( A->sym == SM_UNSYMMETRIC ) || ( A->sym == SM_SYMMETRIC ) );
  assert( A->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO
  assert( A->m == A->n ); // TODO can only handle square matrices at present???

  // Compressed Column Format
  assert( A->jj[0] == 0 );
  assert( A->jj[A->m] == A->nz );

  // TODO orderings? taucs_ccs_order

// TODO not working for unsym matrices...
/*
  static char* options[] =
    { "taucs.factor.numeric=false",
      "taucs.factor.LU=true",
      // TODO taucs.factor.mf or ll for multifrontal or left-looking
      // taucs.factor.ordering chooses the ordering method
      NULL }; // NULL terminated string list

  const int nrhs = 0;
  int ret = taucs_linsolve( &m, &(p->F), nrhs, NULL, NULL, options, NULL ); // factorize
  if(ret != TAUCS_SUCCESS)
    fprintf(stderr,"taucs analysis error: %s\n",taucs_errorcode_to_string_(ret));
  assert(ret == TAUCS_SUCCESS);
*/

  taucs_ccs_matrix m;
  m.m = A->m; // rows
  m.n = A->m; // columns
  m.flags = TAUCS_DOUBLE;
  m.colptr = (int*) A->jj; // colptr
  m.rowind = (int*) A->ii; // row indices
  m.values.d = (taucs_double*) A->dd;

  int* invperm = NULL;
  taucs_ccs_order(&m, &(p->colperm), &invperm, "colamd"); // return void
  free(invperm); // don't care about the inverse permutation
}

void solver_factorize_taucs( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  assert( A != NULL );
  solve_system_taucs_t* const p = s->specific;
  assert( p != NULL );

  // TODO symmetric matrix handling
  if ( A->sym == SM_SYMMETRIC )
    convert_matrix_symmetry( A, BOTH );

  // prepare the matrix
  int ierr = convert_matrix( A, SM_CSC, FIRST_INDEX_ZERO );
  assert( ierr == 0 );
  assert(( A->sym == SM_UNSYMMETRIC ) || ( A->sym == SM_SYMMETRIC ) );
  assert( A->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO
  assert( A->m == A->n ); // TODO can only handle square matrices at present (UMFPACK?)


  // save A for the solve step.. needed for iterative refinement // TODO add A to solve stage to allow iterative refinement!, remove this (and from other solvers)
  p->Arows = A->m;
  p->Acols = A->n;
  p->Ajj = ( int* ) A->jj;
  p->Aii = ( int* ) A->ii;
  p->Add = A->dd;

// TODO currently not working for unsym matrices
/*
  static char* options[] =
    { "taucs.factor.LU=true",
      "taucs.factor.symbolic=false", // use the analyze results
      NULL }; // NULL terminated string list

  const int nrhs = 0;
  int ret = taucs_linsolve( &m, &(p->F), nrhs, NULL, NULL, options, NULL ); // factorize
  if(ret != TAUCS_SUCCESS)
    fprintf(stderr,"taucs factorization error: %s\n",taucs_errorcode_to_string_(ret));
  assert(ret == TAUCS_SUCCESS);

  // TODO currently this re-runs the ordering again, instead of reusing the old one from _anaylze()
*/
  taucs_ccs_matrix m;
  m.m = A->m; // rows
  m.n = A->n; // columns
  m.flags = TAUCS_DOUBLE;
  m.colptr = (int*) A->jj; // colptr
  m.rowind = (int*) A->ii; // row indices
  m.values.d = (taucs_double*) A->dd;

  // use taucs out-of-core solver
  if(p->LU == NULL)
    p->LU = taucs_io_create_multifile("/tmp/taucs-");
  const double memory = 1e6;// TODO
  int ret = taucs_ooc_factor_lu(&m, p->colperm, p->LU, memory);
  if(ret != TAUCS_SUCCESS)
    fprintf(stderr,"taucs evaluation error: %s\n",taucs_errorcode_to_string_(ret));
  assert(ret == TAUCS_SUCCESS);
  
}

// TODO b can be sparse... ??
// TODO can b be a matrix (vs a vector)?
void solver_evaluate_taucs( solver_state_t* s, matrix_t* b, matrix_t* x ) {
  assert( s != NULL );
  assert( b != NULL );
  assert( b != x ); // TODO allow this form
  solve_system_taucs_t* const p = s->specific;
  assert( p != NULL );

  // and we have a valid 'x' and 'b'
  int ierr = convert_matrix( b, DCOL, FIRST_INDEX_ONE ); // TODO support for sparse rhs too
  assert( ierr == 0 );
  assert( b->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO

  // TODO if solver only handles single rhs, loops solver and collect answers...

  // allocate data space // TODO if required
  // TODO move this to master function
  // TODO handle sparse answers??
  if(x != b)
    clear_matrix( x );
  double *const xdd = malloc(( p->Arows * b->n ) * sizeof( double ) ); // TODO or handle complex!

// TODO currently not working for unsym matrices
/*
  static char* options[] = 
    { "taucs.factor.LU=true",
      "taucs.factor.symbolic=false", // use the analyze results
      //"taucs.factor.numeric=false", // use the factorize results
      "taucs.ooc=true", // TAUCS only supports out-of-core solver for unsym matrices -- there is an in-core solver for sym matrices TODO use this when we can
      "taucs.ooc.basename=/tmp/taucs-",
      NULL }; // NULL terminated string list

  int nrhs = b->n;

  taucs_ccs_matrix m;
  m.m = p->Arows; // rows
  m.n = p->Acols; // columns
  m.flags = TAUCS_DOUBLE;
  m.colptr = p->Ajj; // colptr
  m.rowind = p->Aii; // row indices
  m.values.d = (taucs_double*) p->Add;
  int ret = taucs_linsolve( &m, &(p->F), nrhs, xdd, b->dd, options, NULL ); // solve
  if(ret != TAUCS_SUCCESS)
    fprintf(stderr,"taucs evaluation error: %s\n",taucs_errorcode_to_string_(ret));
  assert(ret == TAUCS_SUCCESS);
*/

  // use taucs out-of-core solver
  int ret = taucs_ooc_solve_lu(p->LU, xdd, b->dd);
  if(ret != TAUCS_SUCCESS)
    fprintf(stderr,"taucs evaluation error: %s\n",taucs_errorcode_to_string_(ret));
  assert(ret == TAUCS_SUCCESS);


  // return result
  x->format = DCOL;
  x->sym = SM_UNSYMMETRIC;
  x->data_type = b->data_type;
  x->m = p->Arows;
  x->n = b->n;
  x->nz = x->m * x->n;
  x->dd = xdd;
}

void solver_finalize_taucs( solver_state_t* s ) {
  if ( s == NULL )
    return;

  solve_system_taucs_t* const p = s->specific;

  // release memory
  if ( p != NULL ) {
    if(p->F != NULL) {
      int ret = taucs_linsolve(NULL,&(p->F),0,NULL,NULL,NULL,NULL); // free the factorization
      if(ret != TAUCS_SUCCESS)
        fprintf(stderr,"taucs finalize error: %s\n",taucs_errorcode_to_string_(ret));
      assert(ret == TAUCS_SUCCESS);
    }
    if(p->LU != NULL) {
      int ret = taucs_io_delete(p->LU);
      if(ret != TAUCS_SUCCESS)
        fprintf(stderr,"taucs finalize error: %s\n",taucs_errorcode_to_string_(ret));
      assert(ret == TAUCS_SUCCESS);
    }
    p->F = NULL;
    p->LU = NULL;
  }

  free( p );
  s->specific = NULL; // TODO check this is cleared in all finalize functions...
}
