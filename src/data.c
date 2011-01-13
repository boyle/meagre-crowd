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
#include "data.h"

#include <stdlib.h>
#include <string.h>

void inline free_matrix(matrix_t* m) {
  free(m->dd);
  free(m->ii);
  free(m->jj);
  free(m);
}

void inline clear_matrix(matrix_t* m) {
  free(m->dd);
  free(m->ii);
  free(m->jj);
  *m = (matrix_t){0}; // assign all zeros
}

void inline _normalize_matrix(matrix_t* m);
void inline _normalize_matrix(matrix_t* m) {
  if((m->m == 0) || (m->n == 0)) {
    clear_matrix(m);
    return;
  }
    
  switch(m->format) {
    case CSC:
      m->format = COMPRESSED;
      m->order = COLUMN;
      break;
    case CSR:
      m->format = COMPRESSED;
      m->order = ROW;
      break;
    default: break; // do nothing
  }
}

// deep copy
matrix_t* copy_matrix(matrix_t* m) {
  matrix_t* ret = malloc(sizeof(matrix_t));
  if(ret == NULL) // malloc failed
    return NULL;

  _normalize_matrix(m);
  *ret = *m; // shallow copy
  size_t data_bytes = 0;
  switch(m->data_type) {
    case REAL_SINGLE: data_bytes = sizeof(float); break;
    case REAL_DOUBLE: data_bytes = sizeof(double); break;
    case COMPLEX_SINGLE: data_bytes = 2*sizeof(float); break;
    case COMPLEX_DOUBLE: data_bytes = 2*sizeof(double); break;
    case PATTERN: data_bytes = 0; break;
  }
  // if its only a pattern, or an empty matrix, then there is nothing in 'dd', it's a NULL ptr
  if((m->data_type != PATTERN) && (m->m != 0)) {
    ret->dd = malloc((m->nz)*data_bytes);
    if(ret->dd == NULL) { // malloc failed
      free(ret);
      return NULL;
    }
    memcpy(ret->dd, m->dd, (m->nz)*data_bytes); // memcpy(*dest,*src,n)
  }

  switch(m->format) {
    case COO:
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

    // TODO remove redundancy - remove COMPRESSED and 'order' field in favour of CSR/CSC -- will need to change to DENSE_ROW/DENSE_COLUMN from DENSE
    case CSC: case CSR: case COMPRESSED: 
      switch(m->order) {
        case ROW:

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

	case COLUMN:

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
      }
      break;

      case DENSE:

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
  }

  return ret;
}

int validate_matrix(matrix_t* m) {
  // both row and column must be zero
  if((m->m == 0) || (m->n == 0)) // empty matrix
    return -1;

  // empty matrices don't need to be validated further
  if((m->m == 0) || (m->n == 0)) // empty matrix
    return 0;

  // pattern type matrices can't hold data, only indices
  if((m->data_type == PATTERN) && (m->dd != NULL))
    return -2;

  // if dense, ii, jj = NULL
  if(m->format == DENSE) {
    if(m->ii != NULL)
      return -3;
    if(m->jj != NULL)
      return -4;
  }

  // if CSC/CSR, then nz must match expected value in first & last element of m->ii/jj
  if(m->format == COMPRESSED) {
    if((m->order == ROW) && ((m->ii)[m->m+1] != m->nz))
      return -5;
    if((m->order == COLUMN) && ((m->jj)[m->n+1] != m->nz))
      return -6;
    if((m->order == ROW) && ((m->ii)[0] != 0))
      return -7;
    if((m->order == COLUMN) && ((m->jj)[0] != 0))
      return -8;
  }

  // TODO if symmetric check there aren't any extra values

  // all is okay
  return 0;
}
