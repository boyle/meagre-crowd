/* Meagre-Crowd: A sparse distributed matrix solver testbench for performance benchmarking.
 * Copyright (C) 2011 Alistair Boyle <alistair.js.boyle@gmail.com>
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
#include "solver_superlu_dist.h"
#include "solvers.h"
#include "matrix.h"
#include "matrix_share.h"
#include <superlu_ddefs.h>

#include <stdlib.h> // malloc, free
#include <string.h> // memcpy
#include <assert.h>
#include <stdint.h> // int64_t
#include <unistd.h> // dup -> redirect stdout temporarily

#include <mpi.h>
#include <math.h>

typedef struct {
  superlu_options_t options;
  gridinfo_t grid;
  int active; // is this node active in the superlu grid?
  int rank0; // who is the old rank0 in the new communicator
  SuperMatrix A;
  ScalePermstruct_t scale_permute;
  LUstruct_t lu;
} solve_system_superlu_dist_t;

void solver_init_superlu_dist( solver_state_t* s ) {
  assert( s != NULL );
  solve_system_superlu_dist_t * const p = calloc( 1, sizeof( solve_system_superlu_dist_t ) );
  assert( p != NULL );
  s->specific = p;

  // set default options
  set_default_options_dist(&(p->options));
  if((s->mpi_rank == 0) && (s->verbosity >= 3)) {
    p->options.PrintStat = YES;
    print_options_dist(&(p->options));
  }
  else {
    p->options.PrintStat = NO;
  }

  // determine the size of the MPI communicator
  const MPI_Comm initial_comm = MPI_COMM_WORLD; // TODO pass in communicator instead of assuming MPI_COMM_WORLD
  int size;
  int ret = MPI_Comm_size(initial_comm, &size);
  assert(ret == MPI_SUCCESS);

  // determine size: a 2D grid
  int nprow = floor(sqrt(size));
  int npcol = nprow;
  superlu_gridinit(initial_comm, nprow, npcol, &(p->grid));
  p->active = (p->grid.iam < nprow * npcol);
  if(s->mpi_rank == 0) {
    assert(p->active); // must have the rank=0 node active or we're broken (A is only loaded on rank=0)
    p->rank0 = p->grid.iam;
  }
  ret = MPI_Bcast(&(p->rank0), 1, MPI_INT, p->rank0, p->grid.comm);
  assert(ret == MPI_SUCCESS);
}

// TODO split analyze stage into ordering and symbolic factorization stages?
void solver_analyze_superlu_dist( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  solve_system_superlu_dist_t* const p = s->specific;
  assert( p != NULL );
  if( !p->active ) // check if this grid node is active
    return;

  // setup the A matrix, shared globally
  matrix_t* AA;
  if(s->mpi_rank == 0) {
    AA = copy_matrix(A);
  }
  else {
    AA = malloc_matrix();
  }
  matrix_bcast(AA, p->rank0, p->grid.comm);
  dCreate_CompCol_Matrix_dist(&(p->A), AA->m, AA->n, AA->nz,
                              AA->dd, (int*) AA->ii, (int*) AA->jj, SLU_NC, SLU_D, SLU_GE);
  // last 3 enums are: stype=column-wise(no super-nodes), dtype=double, mtype=general);

  // clear the data pointers since these are now held by p->A and
  // release the rest of the AA matrix pointer
  AA->ii = NULL;
  AA->jj = NULL;
  AA->dd = NULL;
  free_matrix(AA);
  AA = NULL;

  // TODO analyze
}

void solver_factorize_superlu_dist( solver_state_t* s, matrix_t* A ) {
  // TODO factorize
}

void solver_evaluate_superlu_dist( solver_state_t* s, matrix_t* b, matrix_t* x ) {
  assert( s != NULL );
  solve_system_superlu_dist_t* const p = s->specific;
  assert( p != NULL );
  if( !p->active ) // check if this grid node is active
    return;

  assert(b->format == DCOL);
  assert(b->data_type == REAL_DOUBLE);

  // initialize structures
  SuperLUStat_t stat;
  ScalePermstructInit(p->A.nrow, p->A.ncol, &(p->scale_permute));
  LUstructInit(p->A.nrow, p->A.ncol, &(p->lu));
  PStatInit(&stat);

  // setup for solver
  int ldb;
  int nrhs;
  if(s->mpi_rank == 0) {
    nrhs = b->n;
    ldb = b->m;
  }
  int ret;
  ret = MPI_Bcast(&nrhs, 1, MPI_INT, p->rank0, p->grid.comm);
  assert(ret == MPI_SUCCESS);
  ret = MPI_Bcast(&ldb, 1, MPI_INT, p->rank0, p->grid.comm);
  assert(ret == MPI_SUCCESS);
  // now we know the size of the rhs on all nodes, allocate space
  double *bb = doubleMalloc_dist(nrhs*ldb);
  double* berr = doubleMalloc_dist(nrhs*ldb);
  assert(bb != NULL);
  assert(berr != NULL);
  // share the rhs to all nodes
  if(s->mpi_rank == 0) {
    memcpy(bb, b->dd, sizeof(double)*nrhs*ldb);
  }
  ret = MPI_Bcast(&bb, nrhs*ldb, MPI_DOUBLE, p->rank0, p->grid.comm);
  assert(ret == MPI_SUCCESS);

  // call solver
  int info;
  pdgssvx_ABglobal(&(p->options), &(p->A), &(p->scale_permute), bb, ldb, nrhs, &(p->grid),
                   &(p->lu), berr, &stat, &info);

  // TODO calculate inf_norm

  // print statistics
  if(s->verbosity >= 3)
    PStatPrint(&(p->options), &stat, &(p->grid));

  // release structures
  ScalePermstructFree(&(p->scale_permute));
  Destroy_LU(p->A.nrow, &(p->grid), &(p->lu));

  // copy the anser into x
  if(s->mpi_rank == 0) {
    clear_matrix(x);
    x->format = DCOL;
    x->m = ldb; // since the A matrix is square, the rows in b match the columns in A, which match the rows in x
    x->n = nrhs; // matches the b's columns
    x->nz = ldb * nrhs;
    x->dd = malloc(sizeof(double)*(x->m)*(x->n));
    assert(x->dd != NULL);
    memcpy(x->dd, bb, sizeof(double)*(x->m)*(x->n));
    assert(validate_matrix(x) == 0);
  }

  // release data
  PStatFree(&stat);
  SUPERLU_FREE(bb);
  SUPERLU_FREE(berr);
}

void solver_finalize_superlu_dist( solver_state_t* s ) {
  if ( s == NULL )
    return;

  solve_system_superlu_dist_t* const p = s->specific;

  // release memory
  if ( p != NULL ) {

    // release the A matrix
    if( p->active )
      Destroy_CompCol_Matrix_dist(&(p->A));

    // shutdown the MPI grid for superlu
    superlu_gridexit(&(p->grid));
  }
  free( p );
  s->specific = NULL;
}
