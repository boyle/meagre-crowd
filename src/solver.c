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
#include "solver.h"
#include "perftimer.h"
#include "matrix.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "solver_lookup.h"

static inline int _valid_solver(const int solver);
static inline int _valid_solver(const int solver) {
  if((solver >= 0) && (solver < SOLVER_INDEX_COUNT))
    return 1;
  else
    return 0;
}

// lookup functions
int lookup_solver_by_shortname(const char* shortname) {
  assert(shortname != NULL);
  int i;
  for(i=0; i < SOLVER_INDEX_COUNT; i++) {
    if(strncmp(solver_lookup[i].shortname, shortname, SOLVER_SHORTNAME_MAX_LEN) == 0)
      return i;
  }
  return -1;
}

const char* solver2str(const int solver) {
  if(_valid_solver(solver))
    return solver_lookup[solver].name;
  else
    return "<invalid>";
}

void printf_solvers(const unsigned int verbosity) {
  printf("Available solvers:\n");
  int i, j;
  // determine max shortname length so everything can be aligned nicely
  size_t max_shortname_len = 0;
  size_t l [SOLVER_INDEX_COUNT];
  for(i=0; i < SOLVER_INDEX_COUNT; i++) {
    size_t ll = strlen(solver_lookup[i].shortname);
    l[i] = ll;
    if(ll > max_shortname_len)
      max_shortname_len = ll;
  }
  for(i=0; i < SOLVER_INDEX_COUNT; i++) {
    printf("  %s", solver_lookup[i].shortname);
    for(j=0; j < max_shortname_len - l[i]; j++)
      printf(" ");
    printf("    %s %s (%s, %s)\n",
      solver_lookup[i].name,
      solver_lookup[i].version,
      solver_lookup[i].author,
      solver_lookup[i].license);
    if(verbosity >= 1) { // organization, references
      for(j=0; j < max_shortname_len; j++)
        printf(" ");
      printf("      %s, %s\n",
        solver_lookup[i].organization,
        solver_lookup[i].url);
    }
    if(verbosity >= 2) {
      printf("    references:\n%s\n\n",solver_lookup[i].references);
      // TODO capabilities
    }
  }
}



// --------------------------------------------
// can the preferred solver solve this problem?
//   e.g. can the solver only handle Symmetric Postive Definite (SPD) matrices
// returns: 1 yes, 0 no
int solver_can_do(const int solver, matrix_t* A, matrix_t* b) {
  return 1; // TODO something more clever, like a real answer
}

// select the most appropriate solver for this problem
//  - is it small and thus should be solved single-threaded (single processor)
//  - is it moderate and should be solved SMP (shared memory)
//  - is it huge and should be solved MPI (distributed memory)
int select_solver(matrix_t* A, matrix_t* b) {
  return UMFPACK; // TODO something clever
}

// --------------------------------------------
// wrapper function: solve 'A x = b' for 'x'
// calls initialize, analyze, factorize, evaluate, finalize
// returns x, the solution
void solver(const int solver, const int verbosity, const int mpi_rank, matrix_t* A, matrix_t* b, matrix_t* x) {
  solver_state_t* s = solver_init(solver, verbosity, mpi_rank, NULL);
  assert(s != NULL); // malloc failure
  solver_analyze(s, A);
  solver_factorize(s, A);
  solver_evaluate(s, b, x);
  solver_finalize(s);
}

// wrapper function: solve 'A x = b' for 'x' w/o re-initializing solver
// calls analyze, factorize, evaluate
// must call initialize before and finalize after when all done
// returns x, the solution
void solver_solve(solver_state_t* s, matrix_t* A, matrix_t* b, matrix_t* x) {
  solver_analyze(s, A);
  solver_factorize(s, A);
  solver_evaluate(s, b, x);
}

// --------------------------------------------
// initialize and finalize the solver state
solver_state_t* solver_init(const int solver, const int verbosity, const int mpi_rank, perftimer_t* timer) {
  solver_state_t* s = malloc(sizeof(solver_state_t));
  assert(s != NULL);
  
  // configure state
  s->solver = solver;
  s->verbosity = verbosity;
  s->mpi_rank = mpi_rank;
  s->timer = timer;
  s->specific = NULL;
  if(_valid_solver(solver) && (solver_lookup[solver].init != NULL))
    solver_lookup[solver].init(s);

  return s;
}

void solver_finalize(solver_state_t* s) {
  assert(s != NULL);
  // clean up (and deallocate "s->specific" if required)
  const int solver = s->solver;
  if(_valid_solver(solver) && (solver_lookup[solver].analyze != NULL))
    solver_lookup[solver].finalize(s);

  // make sure it won't get deallocated twice by mistake
  s->specific = NULL;

  free(s);
}

// evaluate the patterns in A, doesn't care about the actual values in the matrix (A->dd)
void solver_analyze(solver_state_t* s, matrix_t* A) {
  assert(s != NULL);
  assert(A != NULL);
  perftimer_inc(s->timer,"analyze",-1);
  const int solver = s->solver;
  if(_valid_solver(solver) && (solver_lookup[solver].analyze != NULL))
    solver_lookup[solver].analyze(s, A);
}


// factorize the matrix A, A must have the same pattern of non-zeros at that used in the solver_analyze stage
void solver_factorize(solver_state_t* s, matrix_t* A) {
  assert(s != NULL);
  assert(A != NULL);
  perftimer_inc(s->timer,"factorize",-1);
  const int solver = s->solver;
  if(_valid_solver(solver) && (solver_lookup[solver].factorize != NULL))
    solver_lookup[solver].factorize(s, A);
}


// solve the matrix 'A' for right-hand side 'b'
// returns 'x', the solution
void solver_evaluate(solver_state_t* s, matrix_t* b, matrix_t* x) {
  assert(s != NULL);
  assert(b != NULL);
  assert(x != NULL);
  perftimer_inc(s->timer,"evaluate",-1);
  const int solver = s->solver;
  if(_valid_solver(solver) && (solver_lookup[solver].evaluate != NULL)) {
    solver_lookup[solver].evaluate(s, b, x);
  }
  else {
    clear_matrix(x);
    x->format = INVALID;
  }
  perftimer_inc(s->timer,"done",-1);
}
