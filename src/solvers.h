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
#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "config.h"
#include "perftimer.h"
#include "matrix.h"

enum solvers_t { UMFPACK, MUMPS }; // TODO rm (obsolete)

// --------------------------------------------
// structures and enums
typedef struct {
  int                solver;
  int                mpi_rank;
  int                verbosity;
  perftimer_t*       timer;
  void*              specific; // further solver-specific state
} solver_state_t;



// utility functions for user interface
int lookup_solver_by_shortname( const char* s );
const char* solver2str( const int solver );
void printf_solvers( const unsigned int verbosity );

// --------------------------------------------
// can the preferred solver solve this problem?
//   e.g. can the solver only handle Symmetric Postive Definite (SPD) matrices
// returns: 1 yes, 0 no
int solver_can_do( const int solver, matrix_t* A, matrix_t* b );

// select the most appropriate solver for this problem
//  - is it small and thus should be solved single-threaded (single processor)
//  - is it moderate and should be solved SMP (shared memory)
//  - is it huge and should be solved MPI (distributed memory)
int select_solver( matrix_t* A, matrix_t* b );

// --------------------------------------------
// wrapper function: solve 'A x = b' for 'x'
// calls initialize, analyze, factorize, evaluate, finalize
// returns x, the solution
// TODO cleaner way of passing in MPI info, if required?
void solver( const int solver, const int verbosity, const int mpi_rank, matrix_t* A, matrix_t* b, matrix_t* x );

// wrapper function: solve 'A x = b' for 'x' w/o re-initializing solver
// calls analyze, factorize, evaluate
// must call initialize before and finalize after when all done
// returns x, the solution
void solver_solve( solver_state_t* state, matrix_t* A, matrix_t* b, matrix_t* x );

// --------------------------------------------
// initialize and finalize the solver state
solver_state_t* solver_init( const int solver, const int verbosity, const int mpi_rank, perftimer_t* timer );
void solver_finalize( solver_state_t* p );

// evaluate the patterns in A, doesn't care about the actual values in the matrix (A->dd)
void solver_analyze( solver_state_t* p, matrix_t* A );
// factorize the matrix A, A must have the same pattern of non-zeros at that used in the solver_analyze stage
void solver_factorize( solver_state_t* p, matrix_t* A );
// solve the matrix 'A' for right-hand side 'b'
// returns 'x', the solution
void solver_evaluate( solver_state_t* p, matrix_t* b, matrix_t* x );

#endif
