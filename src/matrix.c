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
#include "matrix.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

// TODO rm these
#include <bebop/util/init.h>
#include <bebop/util/enumerations.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h>
#include <bebop/smc/coo_matrix.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/csc_matrix.h>

inline matrix_t* malloc_matrix() {
  return malloc(sizeof(matrix_t));
}

void inline free_matrix(matrix_t* m) {
  assert(m != NULL);
  free(m->dd);
  free(m->ii);
  free(m->jj);
  free(m);
}

void inline clear_matrix(matrix_t* m) {
  assert(m != NULL);
  free(m->dd);
  free(m->ii);
  free(m->jj);
  *m = (matrix_t){0}; // assign all zeros
}


// deep copy
// TODO const correctness
matrix_t* copy_matrix(matrix_t* m) {
  assert(m != NULL);

  matrix_t* ret = malloc(sizeof(matrix_t));
  if(ret == NULL) // malloc failed
    return NULL;

  *ret = *m; // shallow copy
  size_t data_bytes = 0;
  switch(m->data_type) {
    case REAL_SINGLE: data_bytes = sizeof(float); break;
    case REAL_DOUBLE: data_bytes = sizeof(double); break;
    case COMPLEX_SINGLE: data_bytes = 2*sizeof(float); break;
    case COMPLEX_DOUBLE: data_bytes = 2*sizeof(double); break;
    case SM_PATTERN: data_bytes = 0; break;
  }
  // if its only a pattern, or an empty matrix, then there is nothing in 'dd', it's a NULL ptr
  if((m->data_type != SM_PATTERN) && (m->format != INVALID)) {
    ret->dd = malloc((m->nz)*data_bytes);
    if(ret->dd == NULL) { // malloc failed
      free(ret);
      return NULL;
    }
    memcpy(ret->dd, m->dd, (m->nz)*data_bytes); // memcpy(*dest,*src,n)
  }

  switch(m->format) {
    case SM_COO:
      ret->ii = malloc((m->m)*sizeof(unsigned int));
      if(ret->ii == NULL) { // malloc failed
        free(ret->dd);
	free(ret);
	return NULL;
      }
      memcpy(ret->ii, m->ii, (m->m)*sizeof(unsigned int)); // memcpy(*dest,*src,n)
      ret->jj = malloc((m->n)*sizeof(unsigned int));
      if(ret->jj == NULL) { // malloc failed
        free(ret->ii);
        free(ret->dd);
	free(ret);
	return NULL;
      }
      memcpy(ret->jj, m->jj, (m->n)*sizeof(unsigned int)); // memcpy(*dest,*src,n)
      break;

    case SM_CSR:
      ret->ii = malloc((m->m+1)*sizeof(unsigned int));
      if(ret->ii == NULL) { // malloc failed
	free(ret->dd);
	free(ret);
	return NULL;
      }
      memcpy(ret->ii, m->ii, (m->m+1)*sizeof(unsigned int)); // memcpy(*dest,*src,n)
      ret->jj = malloc((m->nz)*sizeof(unsigned int));
      if(ret->jj == NULL) { // malloc failed
	free(ret->ii);
	free(ret->dd);
	free(ret);
	return NULL;
      }
      memcpy(ret->jj, m->jj, (m->nz)*sizeof(unsigned int)); // memcpy(*dest,*src,n)
      break;

    case SM_CSC:
      ret->jj = malloc((m->n+1)*sizeof(unsigned int));
      if(ret->jj == NULL) { // malloc failed
	free(ret->dd);
	free(ret);
	return NULL;
      }
      memcpy(ret->jj, m->jj, (m->n+1)*sizeof(unsigned int)); // memcpy(*dest,*src,n)
      ret->ii = malloc((m->nz)*sizeof(unsigned int));
      if(ret->ii == NULL) { // malloc failed
	free(ret->jj);
	free(ret->dd);
	free(ret);
	return NULL;
      }
      memcpy(ret->ii, m->ii, (m->nz)*sizeof(unsigned int)); // memcpy(*dest,*src,n)
      break;

    case DROW: case DCOL: // dense matrix
      ret->ii = malloc((m->m)*sizeof(unsigned int));
      if(ret->ii == NULL) { // malloc failed
	free(ret->dd);
	free(ret);
	return NULL;
      }
      memcpy(ret->ii, m->ii, (m->m)*sizeof(unsigned int)); // memcpy(*dest,*src,n)
      ret->jj = malloc((m->n)*sizeof(unsigned int));
      if(ret->jj == NULL) { // malloc failed
	free(ret->ii);
	free(ret->dd);
	free(ret);
	return NULL;
      }
      memcpy(ret->jj, m->jj, (m->n)*sizeof(unsigned int)); // memcpy(*dest,*src,n)
      break;
    case INVALID:
      break;
  }

  return ret;
}

// compare matrices
// returns: zero on match
// TODO compare with-in a given precision (floating pt. data a->dd)
// TODO const correctness
int cmp_matrix(matrix_t* a, matrix_t* b) {
  assert(a != NULL);
  assert(b != NULL);

  // if empty, then don't compare further
  if((a->format == INVALID) && (b->format == INVALID))
    return 0; // matching
  else if((a->format == INVALID) && (b->format == INVALID))
    return -1; // not matching

  // do sizes/types match?
  if((a->m != b->m) || (a->n != b->n) || (a->nz != b->nz) ||
     (a->data_type != b->data_type))
    return -2;

  // TODO deal with symmetry issues (sym, location)

  // convert them to matching formats if required
  matrix_t* bb;
  int copied;
  if((a->format != b->format) || (a->base != b->base)) {
    copied = 1;
    bb = copy_matrix(b);
    assert(bb != NULL);
    int ret = convert_matrix(bb, a->format, a->base);
    if(ret != 0) {
      free_matrix(bb);
      return -3;
    }
  }
  else {
    copied = 0;
    bb = b;
  }

  // decide how many bytes of each ptr to compare
  size_t iilen, jjlen, ddlen, dwidth;
  switch(a->data_type) { // TODO refactor (duplicated code!)
    case REAL_SINGLE:    dwidth = sizeof(float);
    case REAL_DOUBLE:    dwidth = sizeof(double);
    case COMPLEX_SINGLE: dwidth = 2*sizeof(float);
    case COMPLEX_DOUBLE: dwidth = 2*sizeof(double);
    case SM_PATTERN:     dwidth = 0;
  }
  iilen = jjlen = 0; // make compiler happy: hide warning
  switch(a->format) {
    case INVALID:
      assert(0); // shouldn't be able to get here (checked for above)
    case DROW:
    case DCOL:
      iilen = jjlen = 0;
      ddlen = dwidth * a->m * a->n;
      break;
    case SM_COO:
      iilen = jjlen = a->nz;
      ddlen = dwidth * a->nz;
      break;
    case SM_CSC:
      iilen = a->nz;
      jjlen = a->n+1;
      ddlen = dwidth * a->nz;
      break;
    case SM_CSR:
      iilen = a->m+1;
      jjlen = a->nz;
      ddlen = dwidth * a->nz;
      break;
  }

  // now compare data
  if(iilen != 0) {
    assert(a->ii != NULL);
    assert(bb->ii != NULL);
    int ret = memcmp(a->ii, bb->ii, iilen);
    if(ret != 0) {
      if(copied)
        free_matrix(bb);
      return -3;
    }
  }

  if(jjlen != 0) {
    assert(a->jj != NULL);
    assert(bb->jj != NULL);
    int ret = memcmp(a->jj, bb->jj, jjlen);
    if(ret != 0) {
      if(copied)
        free_matrix(bb);
      return -4;
    }
  }

  if(ddlen != 0) {
    assert(a->dd != NULL);
    assert(bb->dd != NULL);
    int ret = memcmp(a->dd, bb->dd, ddlen);
    if(ret != 0) {
      if(copied)
        free_matrix(bb);
      return -5;
    }
  }

  if(copied)
    free_matrix(bb);

  return 0;
}

// local conversion functions
// NOTE: these have no checking built in!
// TODO these currently use BeBOP sparse matrix stuff inside, to be replaced eventually!

// --------------------------
struct sparse_matrix_t* _bebop_input(matrix_t* m, enum sparse_matrix_storage_format_t f);
struct sparse_matrix_t* _bebop_input(matrix_t* m, enum sparse_matrix_storage_format_t f) {
  struct sparse_matrix_t* A = malloc(sizeof(struct sparse_matrix_t));

  enum value_type_t v = REAL;
  switch(m->data_type) {
    case COMPLEX_SINGLE:
    case REAL_SINGLE:
      assert(0); // TODO bebop can't handle this!
      break;
    case REAL_DOUBLE:
      v = REAL;
      break;
    case COMPLEX_DOUBLE:
      v = COMPLEX;
      break;
    case SM_PATTERN:
      v = PATTERN;
      break;
  }

  switch(f) {
    case CSC: {
      A->format = CSC;
      struct csc_matrix_t* p = calloc(1,sizeof(struct csc_matrix_t));
      A->repr = p;
      p->m = m->m;
      p->n = m->n;
      p->nnz = m->nz;
      p->values = m->dd;
      p->rowidx = (int*) m->ii;
      p->colptr = (int*) m->jj;
      p->symmetry_type = UNSYMMETRIC; // TODO
      p->value_type = v;
      p->ownership = USER_DEALLOCATES;
      break;
    }
    case CSR: {
      A->format = CSR;
      struct csr_matrix_t* p = calloc(1,sizeof(struct csr_matrix_t));
      A->repr = p;
      p->m = m->m;
      p->n = m->n;
      p->nnz = m->nz;
      p->values = m->dd;
      p->rowptr = (int*) m->ii;
      p->colidx = (int*) m->jj;
      p->symmetry_type = UNSYMMETRIC; // TODO
      p->value_type = v;
      p->ownership = USER_DEALLOCATES;
      break;
    }
    case COO: {
      A->format = COO;
      struct coo_matrix_t* p = calloc(1,sizeof(struct coo_matrix_t));
      A->repr = p;
      p->m = m->m;
      p->n = m->n;
      p->nnz = m->nz;
      p->val = m->dd;
      p->II = (int*) m->ii;
      p->JJ = (int*) m->jj;
      p->index_base = ZERO; // TODO don't lose this!
      p->symmetry_type = UNSYMMETRIC; // TODO
      p->value_type = v;
      p->ownership = USER_DEALLOCATES;
      break;
    }
    default:
      assert(0); // we don't handle the rest
  }
  return A;
}

// NOTE: the DATA isn't destroyed, just the containers
void _bebop_destroy(struct sparse_matrix_t* A);
void _bebop_destroy(struct sparse_matrix_t* A) {
  free(A->repr);
  free(A);
}
// --------------------------


//
// CSC -> COO
int _csc2coo(matrix_t* m);
int _csc2coo(matrix_t* m) {
  struct sparse_matrix_t* A =_bebop_input(m, CSC);

  int ierr = sparse_matrix_convert(A, COO); assert(ierr == 0);

  struct coo_matrix_t* p = A->repr;
  m->format = COO;
  m->ii = (unsigned int*) p->II;
  m->jj = (unsigned int*) p->JJ;
  m->dd = p->val;

  _bebop_destroy(A);
  return 0;
}

// COO -> CSC
int _coo2csc(matrix_t* m);
int _coo2csc(matrix_t* m) {
  struct sparse_matrix_t* A =_bebop_input(m, COO);

  int ierr = sparse_matrix_convert(A, CSC); assert(ierr == 0);

  struct csc_matrix_t* p = A->repr;
  m->format = CSC;
  m->ii = (unsigned int*) p->rowidx;
  m->jj = (unsigned int*) p->colptr;
  m->dd = p->values;

  _bebop_destroy(A);
  return 0;
}

// CSR -> COO
int _csr2coo(matrix_t* m);
int _csr2coo(matrix_t* m) {
  struct sparse_matrix_t* A =_bebop_input(m, CSR);

  int ierr = sparse_matrix_convert(A, COO); assert(ierr == 0);

  struct coo_matrix_t* p = A->repr;
  m->format = COO;
  m->ii = (unsigned int*) p->II;
  m->jj = (unsigned int*) p->JJ;
  m->dd = p->val;

  _bebop_destroy(A);
  return 0;
}

// COO -> CSR
int _coo2csr(matrix_t* m);
int _coo2csr(matrix_t* m) {
  struct sparse_matrix_t* A =_bebop_input(m, COO);

  int ierr = sparse_matrix_convert(A, CSR); assert(ierr == 0);

  struct csr_matrix_t* p = A->repr;
  m->format = CSR;
  m->ii = (unsigned int*) p->rowptr;
  m->jj = (unsigned int*) p->colidx;
  m->dd = p->values;

  _bebop_destroy(A);
  return 0;
}

// COO -> DROW
int _coo2drow(matrix_t* m);
int _coo2drow(matrix_t* m) {
  assert(0); // TODO implement me
  return 0;
}

// DROW -> COO
int _drow2coo(matrix_t* m);
int _drow2coo(matrix_t* m) {
  assert(0); // TODO implement me
  return 0;
}

// DROW -> DCOL
int _drow2dcol(matrix_t* m);
int _drow2dcol(matrix_t* m) {
  assert(0); // TODO implement me
  return 0;
}

// DCOL -> DROW
int _dcol2drow(matrix_t* m);
int _dcol2drow(matrix_t* m) {
  assert(0); // TODO implement me
  return 0;
}


// convert between formats: some conversions might take more than one step
// non-zero means failure: -1 to/from INVALID, +1 malloc/realloc failed
int convert_matrix(matrix_t* m, enum matrix_format_t f, enum matrix_base_t b) {
  if(m->base != f) {
    int i;
    switch(m->format) {
      case INVALID:
        return -1;
      case DROW:
      case DCOL:
        break; // nothing needs doing
      case SM_COO: // adjust row and col
        if(b == FIRST_INDEX_ZERO) {
          for(i=0;i<m->nz;i++) {
            m->ii[i]--;
            m->jj[i]--;
          }
        }
        else { // FIRST_INDEX_ONE
          for(i=0;i<m->nz;i++) {
            m->ii[i]++;
            m->jj[i]++;
          }
        }
        break;
      case SM_CSC: // adjust col
        if(b == FIRST_INDEX_ZERO) {
          for(i=0;i<m->nz;i++) {
            m->ii[i]--;
          }
        }
        else { // FIRST_INDEX_ONE
          for(i=0;i<m->nz;i++) {
            m->ii[i]++;
          }
        }
        break;
      case SM_CSR: // adjust row
        if(b == FIRST_INDEX_ZERO) {
          for(i=0;i<m->nz;i++) {
            m->jj[i]--;
          }
        }
        else { // FIRST_INDEX_ONE
          for(i=0;i<m->nz;i++) {
            m->jj[i]++;
          }
        }
        break;
    }
  }

  int ret1, ret2, ret3;
  switch(m->format) {
    case INVALID:
      return -1;
    case DROW:
      switch(f) {
        case INVALID:
	  return -1;
	case DROW:
	  return 0; // nothing to do
	case DCOL:
	  ret1 = _drow2dcol(m);
	  return ret1;
	case SM_COO:
	  ret1 = _drow2coo(m);
	  return ret1;
	case SM_CSC:
	  ret1 = _drow2coo(m);
	  ret2 = _coo2csc(m);
	  return (ret1 || ret2);
	case SM_CSR:
	  ret1 = _drow2coo(m);
	  ret2 = _coo2csr(m);
	  return (ret1 || ret2);
      }
    case DCOL:
      switch(f) {
        case INVALID:
	  return -1;
	case DROW:
	  ret1 = _dcol2drow(m);
	  return ret1;
	case DCOL:
	  return 0; // nothing to do
	case SM_COO:
	  ret1 = _dcol2drow(m);
	  ret2 = _drow2coo(m);
	  return (ret1 || ret2);
	case SM_CSC:
	  ret1 = _dcol2drow(m);
	  ret2 = _drow2coo(m);
	  ret3 = _coo2csc(m);
	  return (ret1 || ret2 || ret3);
	case SM_CSR:
	  ret1 = _dcol2drow(m);
	  ret2 = _drow2coo(m);
	  ret3 = _coo2csr(m);
	  return (ret1 || ret2 || ret3);
      }
    case SM_COO:
      switch(f) {
        case INVALID:
	  return -1;
	case DROW:
	  ret1 = _coo2drow(m);
	  return ret1;
	case DCOL:
	  ret1 = _coo2drow(m);
	  ret2 = _drow2dcol(m);
	  return (ret1 || ret2);
	case SM_COO:
	  return 0; // nothing to do
	case SM_CSC:
	  ret1 = _coo2csc(m);
	  return ret1;
	case SM_CSR:
	  ret1 = _coo2csr(m);
	  return ret1;
      }
    case SM_CSC:
      switch(f) {
        case INVALID:
	  return -1;
	case DROW:
	  ret1 = _csc2coo(m);
	  ret2 = _coo2drow(m);
	  return (ret1 || ret2);
	case DCOL:
	  ret1 = _csc2coo(m);
	  ret2 = _coo2drow(m);
	  ret3 = _drow2dcol(m);
	  return (ret1 || ret2 || ret3);
	case SM_COO:
	  ret1 = _csc2coo(m);
	  return ret1;
	case SM_CSC:
	  return 0; // nothing to do
	case SM_CSR:
	  ret1 = _csc2coo(m);
	  ret2 = _coo2csr(m);
	  return (ret1 || ret2);
      }
    case SM_CSR:
      switch(f) {
        case INVALID:
	  return -1;
	case DROW:
	  ret1 = _csr2coo(m);
	  ret2 = _coo2drow(m);
	  return (ret1 || ret2);
	case DCOL:
	  ret1 = _csr2coo(m);
	  ret2 = _coo2drow(m);
	  ret3 = _drow2dcol(m);
	  return (ret1 || ret2 || ret3);
	case SM_COO:
	  ret1 = _csr2coo(m);
	  return ret1;
	case SM_CSC:
	  ret1 = _csr2coo(m);
	  ret2 = _coo2csc(m);
	  return (ret1 || ret2);
	case SM_CSR:
	  return 0; // nothing to do
      }
  }
  assert(0); // shouldn't be able to get here due to returns
  return -2;
}

// test matrix
// returns zero on pass
//   -1: bad size
//   -2: bad ptrs
//   -3: CSR/CSC wrong size in ptr array
//   -4: CSR/CSC bad first ptr entry (expect 0)
int validate_matrix(matrix_t* m) {
  assert(m != NULL);

  // must be invalid if zero size matrix
  if(((m->m == 0) || (m->n == 0)) && (m->format != INVALID)) // empty matrix
    return -1;

  // empty matrices don't need to be validated further
  // but if its empty it aught to match the descriptor: no ptrs, no size
  if(m->format == INVALID) { // empty matrix
    if((m->m != 0) || (m->n != 0) || (m->nz != 0))
      return -1;
    else if((m->ii != NULL) || (m->jj != NULL) || (m->dd != NULL))
      return -2;
    else
      return 0;
  }

  // if there is no data, ptrs should be null
  if((m->nz == 0) && ((m->ii != NULL) || (m->jj != NULL) || (m->dd != NULL)))
    return -2;

  // pattern type matrices can't hold data, only indices
  if((m->data_type == SM_PATTERN) && (m->dd != NULL))
    return -2;

  // if dense, ii, jj = NULL
  if(((m->format == DROW) || (m->format == DCOL)) && ((m->ii != NULL) || (m->jj != NULL)))
    return -2;

  // if CSC/CSR, then nz must match expected value in first & last element of m->ii/jj
  if((m->format == SM_CSR) && ((m->ii)[m->m+1] != m->nz))
    return -3;
  if((m->format == SM_CSC) && ((m->jj)[m->n+1] != m->nz))
    return -3;
  if((m->format == SM_CSR) && ((m->ii)[0] != 0))
    return -4;
  if((m->format == SM_CSC) && ((m->jj)[0] != 0))
    return -4;

  // TODO for CSC/CSR/COO check for duplicate entrys (should be summed)
  // TODO if symmetric check there aren't any extra values in the other triangle
  // TODO check that no matrix indices are outside the declared matrix rows/cols

  // all is okay
  return 0;
}
