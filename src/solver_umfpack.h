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
