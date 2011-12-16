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
#ifndef _MATRIX_SHARE_H_
#define _MATRIX_SHARE_H_

#include "config.h"
#include "matrix.h"
#include <mpi.h>

// broadcast a matrix A to all nodes in the MPI communicator 'comm'
// MPI node 'root' intially holds the original matrix
// afterwards, all nodes hold the whole matrix 'A'
// return: MPI_SUCCESS on success, MPI_ERR_* on failure (see Mpi_Bcast)
int matrix_bcast(matrix_t* A, int root, MPI_Comm comm);

#endif
