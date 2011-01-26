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

  // and now the deep copy (ii, jj, dd)
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
  else { // otherwise, make sure the shallow copy didn't have cruft in it
    ret->dd = NULL;
  }

  switch(m->format) {
    case SM_COO:
      ret->ii = malloc((m->nz)*sizeof(unsigned int));
      if(ret->ii == NULL) { // malloc failed
        free(ret->dd);
	free(ret);
	return NULL;
      }
      memcpy(ret->ii, m->ii, (m->nz)*sizeof(unsigned int)); // memcpy(*dest,*src,n)
      ret->jj = malloc((m->n)*sizeof(unsigned int));
      if(ret->jj == NULL) { // malloc failed
        free(ret->ii);
        free(ret->dd);
	free(ret);
	return NULL;
      }
      memcpy(ret->jj, m->jj, (m->nz)*sizeof(unsigned int)); // memcpy(*dest,*src,n)
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
    case INVALID: // invalid matrix
      ret->ii = NULL;
      ret->jj = NULL;
      break;
  }

  return ret;
}

// from enum, returns width of ea. value in the matrix in bytes
inline size_t _data_width(const enum matrix_data_type_t t);
inline size_t _data_width(const enum matrix_data_type_t t) {
  switch(t) {
    case REAL_SINGLE:    return sizeof(float);
    case REAL_DOUBLE:    return sizeof(double);
    case COMPLEX_SINGLE: return 2*sizeof(float);
    case COMPLEX_DOUBLE: return 2*sizeof(double);
    case SM_PATTERN:     return 0;
  }
  assert(0); // unhittable
  return 0;
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
  size_t iilen, jjlen, ddlen;
  iilen = jjlen = ddlen = 0;
  const size_t dwidth = _data_width(a->data_type);
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
  const size_t dwidth = _data_width(m->data_type);
  void* d_new = calloc((m->m)*(m->n), dwidth);
  if(d_new == NULL)
    return -1; // malloc failure

  // convert from COO to DROW
  int i;
  const unsigned int cols = m->n;
  for(i=0; i < m->nz; i++) {
    // find index in row-major order, given dwidth size entries
    // copy to the appropriate location in the dense array
    const void* src  = (char*) m->dd + i*dwidth;
    void*       dest = (char*) d_new + ((m->ii[i] * cols) + m->jj[i])*dwidth;
    memcpy(dest, src, dwidth);
  }
  // rest of the entries in the array are zero from calloc()

  // now clean up
  free(m->dd);
  free(m->ii);
  free(m->jj);
  m->dd = d_new;
  m->ii = NULL;
  m->jj = NULL;
  m->format = DROW;

  return 0;
}

int _cmp_matrix_entry(void* a, void* b, const enum data_type_t t, const double tol);
int _cmp_matrix_entry(void* a, void* b, const enum data_type_t t, const double tol) {
  switch(t) {
    case REAL_DOUBLE:
      double* aa = (double*) a;
      double* bb = (double*) b;
      if( aa > bb + tol )
        return 1;
      else if( aa < bb - tol )
        return -1;
      else
        return 0;
      break;
    case REAL_SINGLE: // TODO
    case COMPLEX_SINGLE: // TODO
    case COMPLEX_DOUBLE: // TODO
      break;
  }
  assert(0);
  return 0;
}

// DROW -> COO
int _drow2coo(matrix_t* m);
int _drow2coo(matrix_t* m) {
  const double tol = 1e-15; // tolerance: what to approximate as zero when converting // TODO use machine epsilon*2?
  const size_t dwidth = _data_width(m->data_type);
  const void* zero = calloc(2,sizeof(double)); // largest data size: COMPLEX_DOUBLE=2*double

  // allocate maximum size, then realloc later to reduce to the appropriate size ptr
  // data (dd) is already maximum size
  m->ii = malloc((m->m)*(m->n)*sizeof(unsigned int));
  m->jj = malloc((m->m)*(m->n)*sizeof(unsigned int));
  m->nz = 0;

  // in-place compression of data, row-by-row
  int i;
  const unsigned int rows = m->m;
  const unsigned int cols = m->n;
  unsigned int  index = 0;
  void*         d_old = m->dd;
  void*         d_new = m->dd;
  unsigned int* i_new = m->ii;
  unsigned int* j_new = m->jj;
  for(i=0; i<rows; i++) {
    for(j=0; j<cols; j++) {
      // index = (i*cols + j); // row-major indexing
      assert(index == (i*cols + j)); // check indexing is correct

      // TODO this cmp is inefficient, since there's a switch on data type within the loop, its probably better to code 4 loops, one for each data type
      if(_cmp_matrix_entry(d_old, zero, m->data_type, tol) != 0) { // store, otherwise its zero so skip
        // store, if not the same address
	if(d_old != d_new) { // addresses can't overlap
	  memcpy(d_new, d_old, dwidth); // memcpy(dest,src,size)
	}

        d_new += dwidth;
        i_new++;
        j_new++;
      }
      d_old += dwidth;
      index++; // TODO rm unused
    }
  }

  // resize ptr arrays to the correct size, now that
  // we know exactly how many non-zero entries there are
  void* p = realloc(m->ii, m->nz*sizeof(unsigned int));
  if(p != NULL)
    m->ii = p;
  p = realloc(m->jj, m->nz*sizeof(unsigned int));
  if(p != NULL)
    m->jj = p;
  p = realloc(m->dd, m->nz*dwidth);
  if(p != NULL)
    m->dd = p;

  free(zero);
  return 0;
}


// DROW -> DCOL
int _drow2dcol(matrix_t* m);
int _drow2dcol(matrix_t* m) {
  void* d_new = malloc((m->m)*(m->n)*_data_width(m->data_type));
  if(d_new == NULL)
    return -1; // malloc failure

  // swaps rows and columns
  unsigned int i, j;
  void* d_old = m->dd;
  const unsigned int rows = m->m;
  const unsigned int cols = m->n;
  const size_t dwidth = _data_width(m->data_type);
  for(i=0; i < rows; i++) {
    for(j=0; j < cols; j++) {
      const void* src =  (char*)d_old + (i*cols + j)*dwidth;
      void*       dest = (char*)d_new + (j*rows + i)*dwidth;
      memcpy(dest, src, dwidth);
    }
  }

  // swap ptrs
  free(m->dd);
  m->dd = d_new;
  m->format = DCOL;

  return 0;
}

// DCOL -> DROW
int _dcol2drow(matrix_t* m);
int _dcol2drow(matrix_t* m) {
  void* d_new = malloc((m->m)*(m->n)*_data_width(m->data_type));
  if(d_new == NULL)
    return -1; // malloc failure

  // swaps rows and columns
  unsigned int i, j;
  void* d_old = m->dd;
  const unsigned int rows = m->m;
  const unsigned int cols = m->n;
  const size_t dwidth = _data_width(m->data_type);
  // TODO common code: only swapping roll of i,j and indexing between _dcol2drow and _drow2dcol
  for(i=0; i < cols; i++) {
    for(j=0; j < rows; j++) {
      const void* src = (char*)d_old + (i*rows + j)*dwidth;
      void*       dest =(char*)d_new + (j*cols + i)*dwidth;
      memcpy(dest, src, dwidth);
    }
  }

  // swap ptrs
  free(m->dd);
  m->dd = d_new;
  m->format = DROW;

  return 0;
}


// convert between formats: some conversions might take more than one step
// non-zero means failure: -1 to/from INVALID, +1 malloc/realloc failed
int convert_matrix(matrix_t* m, enum matrix_format_t f, enum matrix_base_t b) {
  if(m->base != b) {
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
      return -2;
    case DROW:
      switch(f) {
        case INVALID:
	  return -3;
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
	  return -4;
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
	  return -5;
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
	  return -6;
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
	  return -7;
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
  return -8;
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
