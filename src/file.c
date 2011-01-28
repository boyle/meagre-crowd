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
#include "file.h"
#include <stdio.h> // fprintf
#include <string.h> // strnlen, strcmp
#include <assert.h>

#include "matrix.h"

#include <bebop/util/init.h>
#include <bebop/util/enumerations.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h> // load_sparse_matrix
#include <bebop/smc/coo_matrix.h> // convert to coo

static int _identify_format_from_extension(char* n, enum sparse_matrix_file_format_t* ext, int is_input);

// load a matrix from file "n" into matrix A
// returns 0: success, 1: failure
int load_matrix(char* n, matrix_t* AA) {
  assert(AA != NULL);
  if(n == NULL) {
    fprintf(stderr,"input error: No input specified (-i)\n");
    return 1; // failure
  }

  int retval;
  enum sparse_matrix_file_format_t ext;
  if( (retval = _identify_format_from_extension(n, &ext, 1)) != 0)
    return retval;

  struct sparse_matrix_t* A = load_sparse_matrix(ext, n);
  if(A == NULL) {
    fprintf(stderr,"input error: Failed to load matrix\n");
    return 1;
  }
  int ierr = sparse_matrix_convert(A, COO);
  assert(ierr == 0);
  assert(A->format == COO);
  struct coo_matrix_t* Acoo = A->repr;
  Acoo->ownership = USER_DEALLOCATES;
  assert(Acoo->index_base == ZERO);
  assert(Acoo->symmetry_type == UNSYMMETRIC);
  assert(Acoo->value_type == REAL);

  // store data
  clear_matrix(AA);
  AA->m = Acoo->m;
  AA->n = Acoo->n;
  AA->nz = Acoo->nnz;
  AA->base = FIRST_INDEX_ZERO;
  AA->format = SM_COO;
  AA->sym = SM_UNSYMMETRIC;
  AA->location = BOTH;
  AA->data_type = REAL_DOUBLE;
  AA->dd = Acoo->val;
  AA->ii = (unsigned int*) Acoo->II;
  AA->jj = (unsigned int*) Acoo->JJ;

  // destroy BeBOP's matrix, guts (excepting the data which we now own)
  free(A->repr);

  return 0; // success
}

// save a matrix into file "n" from matrix A
// returns 0: success, 1: failure
int save_matrix(matrix_t* AA, char* n) {
  if(n == NULL) {
    fprintf(stderr,"output error: No output specified (-o)\n");
    return 1; // failure
  }

  int ret;
  enum sparse_matrix_file_format_t ext;
  if( (ret = _identify_format_from_extension(n, &ext, 0)) != 0)
    return ret;

  // make sure we're in the right format (COO) and base 0
  if( (ret = convert_matrix(AA, SM_COO, FIRST_INDEX_ZERO)) != 0)
    return ret;
  assert(AA->format == SM_COO);
  assert(AA->base == FIRST_INDEX_ZERO);

  assert(AA->sym == SM_UNSYMMETRIC);
  assert(AA->data_type == REAL_DOUBLE);

  // copy data into the BeBOP format
  struct sparse_matrix_t A;
  struct coo_matrix_t Acoo;
  A.format = COO;
  A.repr = &Acoo;
  Acoo.m = AA->m;
  Acoo.n = AA->n;
  Acoo.nnz = AA->nz;
  Acoo.II = (signed int*) AA->ii;
  Acoo.JJ = (signed int*) AA->jj;
  Acoo.val = AA->dd;
  Acoo.index_base = ZERO;
  Acoo.symmetry_type = UNSYMMETRIC;
  Acoo.value_type = REAL;
  Acoo.ownership = USER_DEALLOCATES; // don't let BeBOP blow away our data

  save_sparse_matrix(n, &A, ext);
  return 0; // success
}


// returns number of rows in matrix A
inline unsigned int matrix_rows(const matrix_t* const A) {
  if(A == NULL)
    return 0;
  return A->m;
}

// identify the file format from the extension
// return 0: success, 1: failure
// Note: static -- only visible w/in this file
static int _identify_format_from_extension(char* n, enum sparse_matrix_file_format_t* ext, int is_input) {

  size_t s = strnlen(n,100);
  char *e = n + s - 3;
  // strcmp returned match
  if((s > 3) && (strncmp(e,".mm",100) == 0)) {
    *ext = MATRIX_MARKET;
    return 0; // success
  }
  else if((s>3) && (strncmp(e,".hb",100) == 0)) {
    *ext = HARWELL_BOEING;
    if(is_input)
      fprintf(stderr,"input error: Sorry Harwell-Boeing reader is broken\n");
    else
      fprintf(stderr,"output error: Sorry Harwell-Boeing writer is broken\n");
    return 1; // failure
  }
  else if((s > 3) && (strncmp(e,".rb",100) == 0)) {
    *ext = HARWELL_BOEING;
    if(is_input)
      fprintf(stderr,"input error: Sorry Rutherford-Boeing reader is broken\n");
    else
      fprintf(stderr,"output error: Sorry Rutherford-Boeing writer is broken\n");
    return 1; // failure
  }
  else if((s > 4) && (strncmp(e-1,".mat",100) == 0)) {
    *ext = MATLAB;
    if(is_input)
      fprintf(stderr,"input error: Sorry Matlab reader is broken\n");
    else
      fprintf(stderr,"error: Sorry Matlab writer is broken\n");
    return 1; // failure
  } // TODO test if the matlab reader is actually busted
  else{
    if(is_input)
      fprintf(stderr,"input error: Unrecognized file extension\n");
    else
      fprintf(stderr,"output error: Unrecognized file extension\n");
    return 1; // failure
  }
}
