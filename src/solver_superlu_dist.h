/* Meagre-Crowd: A sparse distributed matrix solver testbench for performance benchmarking.
 * Copyright (C) 2011 Alistair Boyle <alistair.js.boyle@gmail.com>
 *
 *     This file is part of Meagre-Crowd.
 *
 *     Meagre-Crowd program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _SOLVER_SUPERLU_DIST_H_
#define _SOLVER_SUPERLU_DIST_H_

#include "config.h"
#include "solvers.h"
#include "matrix.h"

//const char* solver_version_superlu_dist();
void solver_init_superlu_dist( solver_state_t* s );
void solver_analyze_superlu_dist( solver_state_t* s, matrix_t* A );
void solver_factorize_superlu_dist( solver_state_t* s, matrix_t* A );
void solver_evaluate_superlu_dist( solver_state_t* s, matrix_t* b, matrix_t* x );
// TODO void solver_refine_superlu_dist( solver_state_t* s, matrix_t* A, matrix_t* b, matrix_t* x );
void solver_finalize_superlu_dist( solver_state_t* s );

#endif
