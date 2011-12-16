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

#include "config.h"
#include "matrix_share.h"

#include <assert.h>

// broadcast a matrix A to all nodes in the MPI communicator 'comm'
// MPI node 'root' intially holds the original matrix
// afterwards, all nodes hold the whole matrix 'A'
// return: MPI_SUCCESS on success, MPI_ERR_* on failure (see Mpi_Bcast)
int matrix_bcast(matrix_t* A, int root, MPI_Comm comm) {
  assert(A != NULL);

  // get communicator info
  int ret;
  int myrank, size;
  ret = MPI_Comm_rank(comm, &myrank);
  assert(ret == MPI_SUCCESS);
  ret = MPI_Comm_size(comm, &size);
  assert(ret == MPI_SUCCESS);
  if(size <= 1)
    return ret; // nothing to do

  // check our root node is within the communicator
  assert(root < size); // root to big for communicator
  if(myrank == root) {
    assert(A->format == SM_CSC); // can only handle CSC format for now
    assert(A->data_type == REAL_DOUBLE); // con only handle double
  }

  // make sure we aren't going to lose any memory
  if(myrank != root)
    clear_matrix(A);

  // send round one: sizes for malloc
  // TODO this could be more efficient if sending all three, then waiting for completion
  // TODO return an error rather than aborting
  ret = MPI_Bcast(&(A->m),         1, MPI_INT, root, comm);
  assert(ret == MPI_SUCCESS);
  ret = MPI_Bcast(&(A->n),         1, MPI_INT, root, comm);
  assert(ret == MPI_SUCCESS);
  ret = MPI_Bcast(&(A->nz),        1, MPI_INT, root, comm);
  assert(ret == MPI_SUCCESS);
  ret = MPI_Bcast(&(A->base),      1, MPI_INT, root, comm);
  assert(ret == MPI_SUCCESS);
  ret = MPI_Bcast(&(A->format),    1, MPI_INT, root, comm);
  assert(ret == MPI_SUCCESS);
  ret = MPI_Bcast(&(A->sym),       1, MPI_INT, root, comm);
  assert(ret == MPI_SUCCESS);
  ret = MPI_Bcast(&(A->location),  1, MPI_INT, root, comm);
  assert(ret == MPI_SUCCESS);
  ret = MPI_Bcast(&(A->data_type), 1, MPI_INT, root, comm);
  assert(ret == MPI_SUCCESS);

  // allocate memory where required
  if(myrank != root) {
    // TODO refactor: this could be a generic matrix_realloc(A, m_new, n_new, nz_new)
    // assuming CSC format
    A->ii = malloc(A->nz * sizeof(unsigned int)); // row indices
    assert(A->ii != NULL); // TODO malloc error
    A->jj = malloc(((A->n)+1) * sizeof(unsigned int)); // col ptrs
    assert(A->jj != NULL); // TODO malloc error
    A->dd = malloc(A->nz * sizeof(double)); // data
    assert(A->dd != NULL); // TODO malloc error
  }

  // send round two: data
  ret = MPI_Bcast(A->ii, A->nz, MPI_INT, root, comm);
  assert(ret == MPI_SUCCESS);
  ret = MPI_Bcast(A->jj, (A->n) +1, MPI_INT, root, comm);
  assert(ret == MPI_SUCCESS);
  ret = MPI_Bcast(A->dd, A->nz, MPI_DOUBLE, root, comm);
  assert(ret == MPI_SUCCESS);

  return ret;
}

