#ifndef _SOLVER_UMFPACK_H_
#define _SOLVER_UMFPACK_H_

#include "config.h"
#include <bebop/smc/sparse_matrix.h>
#include "args.h"
#include "perftimer.h"


typedef struct {
  int m;
  int n;
  int* Ap;
  int* Ai;
  double* Ax;
  double* b;
  double* x;
} solve_system_umfpack_t;

solve_system_umfpack_t* solver_init_umfpack(struct parse_args* args, perftimer_t* timer, struct sparse_matrix_t* A);

void solver_data_prep_umfpack(solve_system_umfpack_t* p, struct sparse_matrix_t* A, double* b);

double* solver_solve_umfpack(solve_system_umfpack_t* p, struct parse_args* args, perftimer_t* timer);

void solver_finalize_umfpack(solve_system_umfpack_t* p);

#endif
