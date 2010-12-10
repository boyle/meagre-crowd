#include "config.h"
#include "solver_umfpack.h"


#include <stdlib.h> // malloc, free
#include <string.h> // memcpy
#include <assert.h>

#include <bebop/smc/csc_matrix.h>
#include <bebop/util/enumerations.h>

#include <umfpack.h>

// Initialize a MUMPS instance. Use MPI_COMM_WORLD.
// Note: only using A to determine matrix type
solve_system_umfpack_t* solver_init_umfpack(struct parse_args* args, perftimer_t* timer, struct sparse_matrix_t* A) {
  return calloc(1,sizeof(solve_system_umfpack_t));
}

// if NULL, do nothing to A or b
// TODO b can be sparse...
void solver_data_prep_umfpack(solve_system_umfpack_t* p, struct sparse_matrix_t* A, double* b) {
  if(A != NULL) {
    // prepare the matrix
    int ierr = sparse_matrix_convert(A, CSC); assert(ierr == 0);
    struct csc_matrix_t* Acsc = A->repr;
    assert(Acsc != NULL);
    assert(Acsc->symmetry_type == UNSYMMETRIC);
    assert(Acsc->value_type == REAL); // don't handle complex... yet TODO
    // TODO do soemthing with A.ownership, so we can tell bebop to clean itself up, but not have to copy the elements
    // TODO really we should just copy this to be CORRECT/TYPESAFE (not worth being clever...)

    p->m = p->n = Acsc->m;
    assert(Acsc->m == Acsc->n); // TODO can only handle square matrices at present (UMFPACK?)

    // Compressed Column Format
    p->Ap = Acsc->colptr; // start index for column n, column 0 = 0, column N+1 = nz
    p->Ai = Acsc->rowidx; // row indices
    p->Ax = Acsc->values;
    assert(p->Ap[0] == 0);
    assert(p->Ap[p->n] == Acsc->nnz);

    // TODO allocate x
  }

  if(b != NULL) {
    free(p->b); // no-op if NULL
    p->b = malloc(p->n * sizeof(double));
    assert(p->b != NULL); // malloc failure
    memcpy(p->b, b, p->n);
  }
}

void solver_solve_umfpack(solve_system_umfpack_t* p, struct parse_args* args, perftimer_t* timer) {
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
}


void solver_finalize_umfpack(solve_system_umfpack_t* p) {
  free(p->x);
  free(p->b);
  free(p);
}
