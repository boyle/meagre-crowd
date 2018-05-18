/* Meagre-Crowd: A sparse distributed matrix solver testbench for performance benchmarking.
 * Copyright (C) 2010 Alistair Boyle <alistair.js.boyle@gmail.com>
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
#ifndef _FILE_H_
#define _FILE_H_

#include "config.h"
#include "matrix.h"

// load a matrix from file "n" into matrix A
// returns 0: success, <0 failure
int load_matrix( char* n, matrix_t* A );

// save a matrix into file "n" from matrix A
// returns 0: success, 1: failure
int save_matrix( matrix_t* A, char* n );

// returns number of rows in matrix A
unsigned int matrix_rows( const matrix_t* const A );

#endif
