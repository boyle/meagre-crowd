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
#include "solver_umfpack.h"


#include <stdlib.h> // malloc, free
#include <string.h> // memcpy
#include <assert.h>

#include <bebop/smc/sparse_matrix_ops.h>
#include <bebop/smc/csc_matrix.h>
#include <bebop/util/enumerations.h>

#include <umfpack.h>

// Initialize a MUMPS instance. Use MPI_COMM_WORLD.
// Note: only using A to determine matrix type
solve_system_umfpack_t* solver_init_umfpack(struct parse_args* args, perftimer_t* timer, matrix_t* A) {
  return calloc(1,sizeof(solve_system_umfpack_t));
}

// if NULL, do nothing to A or b
// TODO b can be sparse...
void solver_data_prep_umfpack(solve_system_umfpack_t* p, matrix_t* A, matrix_t* b) {
  if(A != NULL) {
    // prepare the matrix
    int ierr = convert_matrix(A, SM_CSC, FIRST_INDEX_ZERO);
    assert(ierr == 0);
    assert(A->sym == SM_UNSYMMETRIC);
    assert(A->data_type == REAL_DOUBLE); // don't handle complex... yet TODO

    assert(A->m == A->n); // TODO can only handle square matrices at present (UMFPACK?)
    p->m = p->n = A->m;

    // Compressed Column Format
    p->Ap = (int*) A->jj; // start index for column n, column 0 = 0, column N+1 = nz
    p->Ai = (int*) A->ii; // row indices
    p->Ax = A->dd;
    assert(p->Ap[0] == 0);
    assert(p->Ap[p->n] == A->nz);

    // allocate x
    p->x = malloc(p->m * sizeof(double));
  }

  if(b != NULL) {
    free(p->b); // no-op if NULL
    p->b = malloc(p->n * sizeof(double));
    assert(p->b != NULL); // malloc failure
    int ret = convert_matrix(b, DROW, FIRST_INDEX_ZERO);
    assert(ret == 0);
    if(A != NULL)
      assert(b->m == A->m);
    assert(b->n == 1); // TOOD MUMPS can handle vectors...
    assert(b->data_type == REAL_DOUBLE);
    memcpy(p->b, b->dd, p->n * sizeof(double));
  }
}

void solver_solve_umfpack(solve_system_umfpack_t* p, struct parse_args* args, perftimer_t* timer, matrix_t* ans) {
  void *Symbolic, *Numeric;

  perftimer_inc(timer,"analyze",-1);
  umfpack_di_symbolic(p->m, p->n, p->Ap, p->Ai, p->Ax, &Symbolic, NULL, NULL);

  perftimer_inc(timer,"factorize",-1);
  umfpack_di_numeric(p->Ap, p->Ai, p->Ax, Symbolic, &Numeric, NULL, NULL);
  umfpack_di_free_symbolic(&Symbolic); // done with Symbolic data

  perftimer_inc(timer,"solve",-1);
  umfpack_di_solve(UMFPACK_A, p->Ap, p->Ai, p->Ax, p->x, p->b, Numeric, NULL, NULL) ;
  umfpack_di_free_numeric(&Numeric);

  perftimer_inc(timer,"done",-1);

  clear_matrix(ans);
  ans->m = p->n;
  ans->n = 1;
  ans->format = DROW;
  ans->data_type = REAL_DOUBLE;
  ans->dd = malloc(p->n * sizeof(double));
  assert(ans->dd != NULL);
  memcpy(ans->dd, p->x, p->n * sizeof(double));
}


void solver_finalize_umfpack(solve_system_umfpack_t* p) {
  free(p->x);
  free(p->b);
  free(p);
}
