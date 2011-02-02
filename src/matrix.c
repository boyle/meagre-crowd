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

static inline int _realloc_arrays( matrix_t* m, size_t nz );

// sort by i=0: row, i=1: column
struct _qsort_coo_data_double {
  unsigned int ii;
  unsigned int jj;
  double       dd;
} _qsort_coo_data_double;

static void _qsort_coo( matrix_t* m, int i );
static int _qsort_coo_cmp_rows_cols( const void* a, const void* b );
static int _qsort_coo_cmp_cols_rows( const void* a, const void* b );

// TODO this is quick, dirty and cheap.. can be much improved and maybe avoid the two extra copies for SPEED
static void _qsort_coo( matrix_t* m, int x ) {
  assert( m->format = SM_COO );
  assert( m->data_type == REAL_DOUBLE );

  struct _qsort_coo_data_double* data = malloc(( m->nz ) * sizeof( _qsort_coo_data_double ) );
  assert( data != NULL );

  int i;
  for ( i = 0;i < m->nz;i++ ) {
    ( data[i] ).ii = m->ii[i];
    ( data[i] ).jj = m->jj[i];
    ( data[i] ).dd = (( double* )m->dd )[i];
  }

  if ( x == 0 )
    qsort( data, m->nz, sizeof( struct _qsort_coo_data_double ), &_qsort_coo_cmp_rows_cols );
  else
    qsort( data, m->nz, sizeof( struct _qsort_coo_data_double ), &_qsort_coo_cmp_cols_rows );

  for ( i = 0;i < m->nz;i++ ) {
    m->ii[i] = data[i].ii;
    m->jj[i] = data[i].jj;
    (( double* )m->dd )[i] = data[i].dd;
  }

  free( data );
}

static int _qsort_coo_cmp_rows_cols( const void* a, const void* b ) {
  const struct _qsort_coo_data_double* aa = a;
  const struct _qsort_coo_data_double* bb = b;
  if (( *aa ).ii < ( *bb ).ii ) {
    return -1;
  }
  else if (( *aa ).ii > ( *bb ).ii ) {
    return + 1;
  }
  else {
    if (( *aa ).jj < ( *bb ).jj )
      return -1;
    else if (( *aa ).jj > ( *bb ).jj )
      return + 1;
    else
      return 0;
  }
}

static int _qsort_coo_cmp_cols_rows( const void* a, const void* b ) {
  const struct _qsort_coo_data_double* aa = a;
  const struct _qsort_coo_data_double* bb = b;
  if (( *aa ).jj < ( *bb ).jj ) {
    return -1;
  }
  else if (( *aa ).jj > ( *bb ).jj ) {
    return + 1;
  }
  else {
    if (( *aa ).ii < ( *bb ).ii )
      return -1;
    else if (( *aa ).ii > ( *bb ).ii )
      return + 1;
    else
      return 0;
  }
}



inline matrix_t* malloc_matrix() {
  matrix_t* m = malloc( sizeof( matrix_t ) );
  // make this safe to free_matrix() whatever comes out of this
  if ( m != NULL ) {
    m->format = INVALID;
    m->ii = NULL;
    m->jj = NULL;
    m->dd = NULL;
  }
  return m;
}

void inline free_matrix( matrix_t* m ) {
  if ( m != NULL ) {
    free( m->dd );
    free( m->ii );
    free( m->jj );
    free( m );
  }
}

void inline clear_matrix( matrix_t* m ) {
  assert( m != NULL );
  free( m->dd );
  free( m->ii );
  free( m->jj );
  *m = ( matrix_t ) {
    0
  }
  ; // assign all zeros
}

// from enum, returns width of ea. value in the matrix in bytes
inline size_t _data_width( const enum matrix_data_type_t t );
inline size_t _data_width( const enum matrix_data_type_t t ) {
  switch ( t ) {
    case REAL_SINGLE:
      return sizeof( float );
    case REAL_DOUBLE:
      return sizeof( double );
    case COMPLEX_SINGLE:
      return 2*sizeof( float );
    case COMPLEX_DOUBLE:
      return 2*sizeof( double );
    case SM_PATTERN:
      return 0;
  }
  assert( 0 ); // unhittable
  return 0;
}

// deep copy
// TODO const correctness
matrix_t* copy_matrix( matrix_t* m ) {
  assert( m != NULL );

  matrix_t* ret = malloc( sizeof( matrix_t ) );
  if ( ret == NULL ) // malloc failed
    return NULL;

  *ret = *m; // shallow copy

  // and now the deep copy (ii, jj, dd)
  const size_t dwidth = _data_width( m->data_type );
  // if its only a pattern, or an empty matrix, then there is nothing in 'dd', it's a NULL ptr
  if (( m->data_type != SM_PATTERN ) && ( m->format != INVALID ) && ( m->nz != 0 ) ) {

    size_t n; // entries to copy
    if (( m->format == DROW ) || ( m->format == DCOL ) )
      n = ( m->m ) * ( m->n ); // rows*cols entries
    else
      n = m->nz;

    ret->dd = malloc( n * dwidth ); // nz entries

    // malloc failed
    if ( ret->dd == NULL ) {
      free( ret );
      return NULL;
    }

    memcpy( ret->dd, m->dd, n*dwidth ); // memcpy(*dest,*src,n)
  }
  else { // otherwise, make sure the shallow copy didn't have cruft in it
    ret->dd = NULL;
  }

  switch ( m->format ) {
    case SM_COO:
      ret->ii = malloc(( m->nz ) * sizeof( unsigned int ) );
      if ( ret->ii == NULL ) { // malloc failed
        free( ret->dd );
        free( ret );
        return NULL;
      }
      memcpy( ret->ii, m->ii, ( m->nz )*sizeof( unsigned int ) ); // memcpy(*dest,*src,n)
      ret->jj = malloc(( m->n ) * sizeof( unsigned int ) );
      if ( ret->jj == NULL ) { // malloc failed
        free( ret->ii );
        free( ret->dd );
        free( ret );
        return NULL;
      }
      memcpy( ret->jj, m->jj, ( m->nz )*sizeof( unsigned int ) ); // memcpy(*dest,*src,n)
      break;

    case SM_CSR:
      ret->ii = malloc(( m->m + 1 ) * sizeof( unsigned int ) );
      if ( ret->ii == NULL ) { // malloc failed
        free( ret->dd );
        free( ret );
        return NULL;
      }
      memcpy( ret->ii, m->ii, ( m->m + 1 )*sizeof( unsigned int ) ); // memcpy(*dest,*src,n)
      ret->jj = malloc(( m->nz ) * sizeof( unsigned int ) );
      if ( ret->jj == NULL ) { // malloc failed
        free( ret->ii );
        free( ret->dd );
        free( ret );
        return NULL;
      }
      memcpy( ret->jj, m->jj, ( m->nz )*sizeof( unsigned int ) ); // memcpy(*dest,*src,n)
      break;

    case SM_CSC:
      ret->jj = malloc(( m->n + 1 ) * sizeof( unsigned int ) );
      if ( ret->jj == NULL ) { // malloc failed
        free( ret->dd );
        free( ret );
        return NULL;
      }
      memcpy( ret->jj, m->jj, ( m->n + 1 )*sizeof( unsigned int ) ); // memcpy(*dest,*src,n)
      ret->ii = malloc(( m->nz ) * sizeof( unsigned int ) );
      if ( ret->ii == NULL ) { // malloc failed
        free( ret->jj );
        free( ret->dd );
        free( ret );
        return NULL;
      }
      memcpy( ret->ii, m->ii, ( m->nz )*sizeof( unsigned int ) ); // memcpy(*dest,*src,n)
      break;

    case DROW:
    case DCOL: // dense matrix
    case INVALID: // invalid matrix
      ret->ii = NULL;
      ret->jj = NULL;
      break;
  }

  return ret;
}

// compare matrices
// returns: zero on match
// TODO compare with-in a given precision (floating pt. data a->dd)
// TODO const correctness
int cmp_matrix( matrix_t* a, matrix_t* b ) {
  assert( a != NULL );
  assert( b != NULL );

  // if empty, then don't compare further
  if (( a->format == INVALID ) && ( b->format == INVALID ) )
    return 0; // matching
  else if (( a->format == INVALID ) || ( b->format == INVALID ) )
    return -1; // not matching

  // do sizes/types match?
  if (( a->m != b->m ) || ( a->n != b->n ) ||
      ( a->data_type != b->data_type ) )
    return -2;

  // TODO deal with symmetry issues (sym, location)
  // can't match if symmetry type doesn't match...
  // unless its an undetected symmetric matrix (unsymmetric only)
  if (( a->sym != b->sym ) && ( b->sym != SM_UNSYMMETRIC ) )
    return -3;

  // convert them to matching formats if required
  matrix_t* bb;
  int copied;
  if (( a->format != b->format ) || ( a->base != b->base ) ||
      (( a->sym != b->sym ) && ( b->sym == SM_UNSYMMETRIC ) ) || // maybe b is symmetric?
      (( a->sym == b->sym ) && ( b->sym != SM_UNSYMMETRIC ) && ( a->location != b->location ) ) ) {
    copied = 1;
    bb = copy_matrix( b );
    assert( bb != NULL );
    int ret = convert_matrix( bb, a->format, a->base );
    if ( ret != 0 ) {
      free_matrix( bb );
      return -5;
    }
  }
  else {
    copied = 0;
    bb = b;
  }

  // try detecting symmetry
  if (( a->sym != bb->sym ) && ( bb->sym == SM_UNSYMMETRIC ) )
    assert( detect_matrix_symmetry( bb ) == 0 );

  // if there is symmetry, make sure its in the same format
  if (( a->sym == bb->sym ) && ( bb->sym != SM_UNSYMMETRIC ) )
    assert( convert_matrix_symmetry( bb, a->location ) == 0 );


  if ( a->nz != bb->nz ) {
    if ( copied )
      free_matrix( bb );
    return -7;
  }

  // decide how many bytes of each ptr to compare
  size_t iilen, jjlen, ddlen;
  iilen = jjlen = ddlen = 0;
  const size_t dwidth = _data_width( a->data_type );
  switch ( a->format ) {
    case INVALID:
      assert( 0 ); // shouldn't be able to get here (checked for above)
    case DROW:
    case DCOL:
      iilen = jjlen = 0;
      ddlen = dwidth * a->m * a->n;
      break;
    case SM_COO:
      iilen = jjlen = a->nz;
      ddlen = dwidth * a->nz;
      // TODO this is a hack... really need to sort both matrices for all data types except DROW/DCOL! (need to copy both a and b...)
      _qsort_coo( a, 0 ); // TODO rm
      _qsort_coo( bb, 0 ); // TODO rm
      break;
    case SM_CSC:
      iilen = a->nz;
      jjlen = a->n + 1;
      ddlen = dwidth * a->nz;
      break;
    case SM_CSR:
      iilen = a->m + 1;
      jjlen = a->nz;
      ddlen = dwidth * a->nz;
      break;
  }

  // now compare data
  if ( iilen != 0 ) {
    assert( a->ii != NULL );
    assert( bb->ii != NULL );
    int ret = memcmp( a->ii, bb->ii, iilen );
    if ( ret != 0 ) {
      if ( copied )
        free_matrix( bb );
      return -10;
    }
  }

  if ( jjlen != 0 ) {
    assert( a->jj != NULL );
    assert( bb->jj != NULL );
    int ret = memcmp( a->jj, bb->jj, jjlen );
    if ( ret != 0 ) {
      if ( copied )
        free_matrix( bb );
      return -11;
    }
  }

  if ( ddlen != 0 ) {
    assert( a->dd != NULL );
    assert( bb->dd != NULL );


    int ret;
//    ret = memcmp(a->dd, bb->dd, ddlen); // TODO compare with some tolerance
    switch ( a->data_type ) {
      case REAL_DOUBLE: {
        const double tol = 1e-15; // tolerance: what to approximate as zero when converting // TODO use machine epsilon*2?
        double* d_old = a->dd;
        double* d_new = bb->dd;
        int i;
        for ( i = 0; i < a->nz; i++ ) {
          if (( *d_old < *d_new - tol ) || ( *d_old > *d_new + tol ) ) { // if !zero store, otherwise skip
            ret = -1;
            break;
          }
          d_old++;
          d_new++;
        }
        ret = 0;
        break;
      }
      case REAL_SINGLE:
      case COMPLEX_DOUBLE:
      case COMPLEX_SINGLE:
      case SM_PATTERN:
        assert( 0 ); // TODO
        break;
    }

    if ( ret != 0 ) {
      if ( copied )
        free_matrix( bb );
      return -12;
    }
  }

  if ( copied )
    free_matrix( bb );

  return 0;
}

// local conversion functions
// NOTE: these have no checking built in!
// TODO these currently use BeBOP sparse matrix stuff inside, to be replaced eventually!

// --------------------------
// NOTE: there is memory loss here -- I think its within BeBOP but its hard to see w/o removing it

struct sparse_matrix_t* _bebop_input( matrix_t* m, enum sparse_matrix_storage_format_t f );
struct sparse_matrix_t* _bebop_input( matrix_t* m, enum sparse_matrix_storage_format_t f ) {
  struct sparse_matrix_t* A = malloc( sizeof( struct sparse_matrix_t ) );

  enum value_type_t v = REAL;
  switch ( m->data_type ) {
    case COMPLEX_SINGLE:
    case REAL_SINGLE:
      assert( 0 ); // TODO bebop can't handle this!
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

  switch ( f ) {
    case CSC: {
      A->format = CSC;
      struct csc_matrix_t* p = calloc( 1, sizeof( struct csc_matrix_t ) );
      A->repr = p;
      p->m = m->m;
      p->n = m->n;
      p->nnz = m->nz;
      p->values = m->dd;
      p->rowidx = ( int* ) m->ii;
      p->colptr = ( int* ) m->jj;
      p->symmetry_type = UNSYMMETRIC; // TODO
      p->value_type = v;
      p->ownership = USER_DEALLOCATES;
      break;
    }
    case CSR: {
      A->format = CSR;
      struct csr_matrix_t* p = calloc( 1, sizeof( struct csr_matrix_t ) );
      A->repr = p;
      p->m = m->m;
      p->n = m->n;
      p->nnz = m->nz;
      p->values = m->dd;
      p->rowptr = ( int* ) m->ii;
      p->colidx = ( int* ) m->jj;
      p->symmetry_type = UNSYMMETRIC; // TODO
      p->value_type = v;
      p->ownership = USER_DEALLOCATES;
      break;
    }
    case COO: {
      A->format = COO;
      struct coo_matrix_t* p = calloc( 1, sizeof( struct coo_matrix_t ) );
      A->repr = p;
      p->m = m->m;
      p->n = m->n;
      p->nnz = m->nz;
      p->val = m->dd;
      p->II = ( int* ) m->ii;
      p->JJ = ( int* ) m->jj;
      if ( m->base == FIRST_INDEX_ZERO )
        p->index_base = ZERO;
      else
        p->index_base = ONE;
      p->symmetry_type = UNSYMMETRIC; // TODO
      p->value_type = v;
      p->ownership = USER_DEALLOCATES;
      break;
    }
    default:
      assert( 0 ); // we don't handle the rest
  }
  return A;
}

// NOTE: the DATA isn't destroyed, just the containers
void _bebop_destroy( struct sparse_matrix_t* A );
void _bebop_destroy( struct sparse_matrix_t* A ) {
  free( A->repr );
  free( A );
}
// --------------------------


//
// CSC -> COO
int _csc2coo( matrix_t* m );
int _csc2coo( matrix_t* m ) {
  assert( m->format == SM_CSC );

  // BeBOP only handles zero-based CSC matrices
  const enum matrix_base_t old_base = m->base;
  int ret;
  ret = convert_matrix( m, m->format, FIRST_INDEX_ZERO );
  assert( ret == 0 ); // should never fail

  // store, ready for BeBOP
  struct sparse_matrix_t* A = _bebop_input( m, CSC );

  int ierr = sparse_matrix_convert( A, COO );
  assert( ierr == 0 ); // might fail on malloc?

  // store new matrix
  struct coo_matrix_t* p = A->repr;
  m->format = SM_COO;
  m->ii = ( unsigned int* ) p->II;
  m->jj = ( unsigned int* ) p->JJ;
  m->dd = p->val;

  // clean up after BeBOP
  _bebop_destroy( A );

  // now convert back to appropriate base
  ret = convert_matrix( m, m->format, old_base );
  assert( ret == 0 ); // should never fail

  return 0;
}

// COO -> CSC
int _coo2csc( matrix_t* m );
int _coo2csc( matrix_t* m ) {
  assert( m->format == SM_COO );

  // BeBOP only handles zero-based CSC matrices
  const enum matrix_base_t old_base = m->base;
  int ret;
  ret = convert_matrix( m, m->format, FIRST_INDEX_ZERO );
  assert( ret == 0 ); // should never fail

  // store, ready for BeBOP
  struct sparse_matrix_t* A = _bebop_input( m, COO );

  int ierr = sparse_matrix_convert( A, CSC );
  assert( ierr == 0 ); // could fail on malloc?

  // store new matrix
  struct csc_matrix_t* p = A->repr;
  m->format = SM_CSC;
  m->ii = ( unsigned int* ) p->rowidx;
  m->jj = ( unsigned int* ) p->colptr;
  m->dd = p->values;

  // stupid BeBOP doesn't allocate the last value for CSC
  assert( m->jj[m->n] == m->nz );
  assert( m->jj[0] == 0 );

  // clean up after BeBOP
  _bebop_destroy( A );

  // now convert back to appropriate base
  ret = convert_matrix( m, m->format, old_base );
  assert( ret == 0 ); // should never fail

  return 0;
}

// CSR -> COO
int _csr2coo( matrix_t* m );
int _csr2coo( matrix_t* m ) {
  assert( m->format == SM_CSR );

  // BeBOP only handles zero-based CSC matrices
  const enum matrix_base_t old_base = m->base;
  int ret;
  ret = convert_matrix( m, m->format, FIRST_INDEX_ZERO );
  assert( ret == 0 ); // should never fail

  // store, ready for BeBOP
  struct sparse_matrix_t* A = _bebop_input( m, CSR );

  int ierr = sparse_matrix_convert( A, COO );
  assert( ierr == 0 ); // could fail on malloc?

  // store new matrix
  struct coo_matrix_t* p = A->repr;
  m->format = SM_COO;
  m->ii = ( unsigned int* ) p->II;
  m->jj = ( unsigned int* ) p->JJ;
  m->dd = p->val;

  // clean up after BeBOP
  _bebop_destroy( A );

  // now convert back to appropriate base
  ret = convert_matrix( m, m->format, old_base );
  assert( ret == 0 ); // should never fail

  return 0;
}

// COO -> CSR
int _coo2csr( matrix_t* m );
int _coo2csr( matrix_t* m ) {
  assert( m->format == SM_COO );

  // BeBOP only handles zero-based CSC matrices
  const enum matrix_base_t old_base = m->base;
  int ret;
  ret = convert_matrix( m, m->format, FIRST_INDEX_ZERO );
  assert( ret == 0 ); // should never fail

  // store, ready for BeBOP
  struct sparse_matrix_t* A = _bebop_input( m, COO );

  int ierr = sparse_matrix_convert( A, CSR );
  assert( ierr == 0 ); // could fail on malloc?

  struct csr_matrix_t* p = A->repr;
  m->format = SM_CSR;
  m->ii = ( unsigned int* ) p->rowptr;
  m->jj = ( unsigned int* ) p->colidx;
  m->dd = p->values;

  // clean up after BeBOP
  _bebop_destroy( A );

  // now convert back to appropriate base
  ret = convert_matrix( m, m->format, old_base );
  assert( ret == 0 ); // should never fail

  return 0;
}

// COO -> DROW
int _coo2drow( matrix_t* m );
int _coo2drow( matrix_t* m ) {
  assert( m->format == SM_COO );
  assert( m->base == FIRST_INDEX_ZERO );
  const size_t dwidth = _data_width( m->data_type );
  void* d_new = calloc(( m->m ) * ( m->n ), dwidth );
  if ( d_new == NULL )
    return -1; // malloc failure

  // convert from COO to DROW
  int i;
  const unsigned int cols = m->n;
  for ( i = 0; i < m->nz; i++ ) {
    // find index in row-major order, given dwidth size entries
    // copy to the appropriate location in the dense array
    const void* src  = ( char* ) m->dd + i * dwidth;
    void*       dest = ( char* ) d_new + (( m->ii[i] * cols ) + m->jj[i] ) * dwidth;
    memcpy( dest, src, dwidth );
  }
  // rest of the entries in the array are zero from calloc()

  // now clean up
  free( m->dd );
  free( m->ii );
  free( m->jj );
  m->dd = d_new;
  m->ii = NULL;
  m->jj = NULL;
  m->format = DROW;
  m->nz = m->m * m->n; // not really valid, but might as well set it to a sane value

  return 0;
}


// DROW -> COO
int _drow2coo( matrix_t* m, const enum matrix_base_t b );
int _drow2coo( matrix_t* m, const enum matrix_base_t b ) {
  assert( m->format == DROW );
  const double tol = 1e-15; // tolerance: what to approximate as zero when converting // TODO use machine epsilon*2?
  const size_t dwidth = _data_width( m->data_type );

  // allocate maximum size, then realloc later to reduce to the appropriate size ptr
  // data (dd) is already maximum size
  m->ii = malloc(( m->m ) * ( m->n ) * sizeof( unsigned int ) );
  if ( m->ii == NULL )
    return -1;
  m->jj = malloc(( m->m ) * ( m->n ) * sizeof( unsigned int ) );
  if ( m->jj == NULL ) {
    free( m->ii );
    return -1;
  }

  // in-place compression of data, row-by-row
  const unsigned int rows = m->m;
  const unsigned int cols = m->n;
  unsigned int* i_new = m->ii;
  unsigned int* j_new = m->jj;
  int i, j;

  int nz = 0;
  switch ( m->data_type ) {
    case REAL_DOUBLE: {
      double* d_old = m->dd;
      double* d_new = m->dd;
      for ( i = 0; i < rows; i++ ) { // (i initialized already)
        for ( j = 0; j < cols; j++ ) { // (j initialized already)
          // index = (i*cols + j); // row-major indexing
          if (( *d_old < 0.0 - tol ) || ( *d_old > 0.0 + tol ) ) { // if !zero store, otherwise skip
            // store, if not the same address
            if ( d_old != d_new ) // addresses can't overlap for memcpy
              memcpy( d_new, d_old, dwidth ); // memcpy(dest,src,size)

            if ( b == FIRST_INDEX_ONE ) {// then loop variables have to be +1
              *i_new = i + 1;
              *j_new = j + 1;
            }
            else { // FIRST_INDEX_ZERO
              *i_new = i;
              *j_new = j;
            }

            d_new++;
            i_new++;
            j_new++;
            nz++; // update the entry count
          }
          d_old++;
        }
      }
    }
    break;
    case REAL_SINGLE:
    case COMPLEX_DOUBLE:
    case COMPLEX_SINGLE:
    case SM_PATTERN:
      assert( 0 );
      break;
  }

  // optimization!
  // resize ptr arrays to the correct size, now that
  // we know exactly how many non-zero entries there are
  // -- if these fail we still have the original ptr
  _realloc_arrays( m, nz );
  m->format = SM_COO;
  m->base = b;

  return 0;
}


// DROW -> DCOL
int _drow2dcol( matrix_t* m );
int _drow2dcol( matrix_t* m ) {
  assert( m->format == DROW );
  void* d_new = malloc(( m->m ) * ( m->n ) * _data_width( m->data_type ) );
  if ( d_new == NULL )
    return -1; // malloc failure

  // TODO dcol == drow if rows or cols == 1, short-circuit in that case and change format
  // TODO do this going the other way too dcol -> drow

  // swaps rows and columns
  unsigned int i, j;
  void* d_old = m->dd;
  const unsigned int rows = m->m;
  const unsigned int cols = m->n;
  const size_t dwidth = _data_width( m->data_type );
  for ( i = 0; i < rows; i++ ) {
    for ( j = 0; j < cols; j++ ) {
      const void* src = ( char* )d_old + ( i * cols + j ) * dwidth;
      void*       dest = ( char* )d_new + ( j * rows + i ) * dwidth;
      memcpy( dest, src, dwidth );
    }
  }

  // swap ptrs
  free( m->dd );
  m->dd = d_new;
  m->format = DCOL;

  return 0;
}

// DCOL -> DROW
int _dcol2drow( matrix_t* m );
int _dcol2drow( matrix_t* m ) {
  assert( m->format == DCOL );
  void* d_new = malloc(( m->m ) * ( m->n ) * _data_width( m->data_type ) );
  if ( d_new == NULL )
    return -1; // malloc failure

  // swaps rows and columns
  unsigned int i, j;
  void* d_old = m->dd;
  const unsigned int rows = m->m;
  const unsigned int cols = m->n;
  const size_t dwidth = _data_width( m->data_type );
  // TODO common code: only swapping roll of i,j and indexing between _dcol2drow and _drow2dcol
  for ( i = 0; i < cols; i++ ) {
    for ( j = 0; j < rows; j++ ) {
      const void* src = ( char* )d_old + ( i * rows + j ) * dwidth;
      void*       dest = ( char* )d_new + ( j * cols + i ) * dwidth;
      memcpy( dest, src, dwidth );
    }
  }

  // swap ptrs
  free( m->dd );
  m->dd = d_new;
  m->format = DROW;

  return 0;
}


// convert between formats: some conversions might take more than one step
// non-zero means failure: -1 to/from INVALID, +1 malloc/realloc failed
int convert_matrix( matrix_t* m, enum matrix_format_t f, enum matrix_base_t b ) {
  // convert to base 0 if going to dense format
  if (( f == DROW ) || ( f == DCOL ) )
    b = FIRST_INDEX_ZERO;

  // do base conversion
  if ( m->base != b ) {
    int i;
    switch ( m->format ) {
      case INVALID:
        return -1;
      case DROW:
      case DCOL:
        break; // nothing needs doing
      case SM_COO: // adjust row and col
        if ( b == FIRST_INDEX_ZERO ) {
          for ( i = 0;i < m->nz;i++ ) {
            m->ii[i]--;
            m->jj[i]--;
          }
        }
        else { // FIRST_INDEX_ONE
          for ( i = 0;i < m->nz;i++ ) {
            m->ii[i]++;
            m->jj[i]++;
          }
        }
        break;
      case SM_CSC: // adjust col
        if ( b == FIRST_INDEX_ZERO ) {
          for ( i = 0;i < m->nz;i++ ) {
            m->ii[i]--;
          }
        }
        else { // FIRST_INDEX_ONE
          for ( i = 0;i < m->nz;i++ ) {
            m->ii[i]++;
          }
        }
        break;
      case SM_CSR: // adjust row
        if ( b == FIRST_INDEX_ZERO ) {
          for ( i = 0;i < m->nz;i++ ) {
            m->jj[i]--;
          }
        }
        else { // FIRST_INDEX_ONE
          for ( i = 0;i < m->nz;i++ ) {
            m->jj[i]++;
          }
        }
        break;
    }
    m->base = b;
  }

  int ret1, ret2, ret3;
  switch ( m->format ) {
    case INVALID:
      return -2;
    case DROW:
      switch ( f ) {
        case INVALID:
          return -3;
        case DROW:
          return 0; // nothing to do
        case DCOL:
          ret1 = _drow2dcol( m );
          return ret1;
        case SM_COO:
          ret1 = _drow2coo( m, b );
          return ret1;
        case SM_CSC:
          ret1 = _drow2coo( m, b );
          ret2 = _coo2csc( m );
          return ( ret1 || ret2 );
        case SM_CSR:
          ret1 = _drow2coo( m, b );
          ret2 = _coo2csr( m );
          return ( ret1 || ret2 );
      }
    case DCOL:
      switch ( f ) {
        case INVALID:
          return -4;
        case DROW:
          ret1 = _dcol2drow( m );
          return ret1;
        case DCOL:
          return 0; // nothing to do
        case SM_COO:
          ret1 = _dcol2drow( m );
          ret2 = _drow2coo( m, b );
          return ( ret1 || ret2 );
        case SM_CSC:
          ret1 = _dcol2drow( m );
          ret2 = _drow2coo( m, b );
          ret3 = _coo2csc( m );
          return ( ret1 || ret2 || ret3 );
        case SM_CSR:
          ret1 = _dcol2drow( m );
          ret2 = _drow2coo( m, b );
          ret3 = _coo2csr( m );
          return ( ret1 || ret2 || ret3 );
      }
    case SM_COO:
      switch ( f ) {
        case INVALID:
          return -5;
        case DROW:
          ret1 = _coo2drow( m );
          return ret1;
        case DCOL:
          ret1 = _coo2drow( m );
          ret2 = _drow2dcol( m );
          return ( ret1 || ret2 );
        case SM_COO:
          return 0; // nothing to do
        case SM_CSC:
          ret1 = _coo2csc( m );
          return ret1;
        case SM_CSR:
          ret1 = _coo2csr( m );
          return ret1;
      }
    case SM_CSC:
      switch ( f ) {
        case INVALID:
          return -6;
        case DROW:
          ret1 = _csc2coo( m );
          ret2 = _coo2drow( m );
          return ( ret1 || ret2 );
        case DCOL:
          ret1 = _csc2coo( m );
          ret2 = _coo2drow( m );
          ret3 = _drow2dcol( m );
          return ( ret1 || ret2 || ret3 );
        case SM_COO:
          ret1 = _csc2coo( m );
          return ret1;
        case SM_CSC:
          return 0; // nothing to do
        case SM_CSR:
          ret1 = _csc2coo( m );
          ret2 = _coo2csr( m );
          return ( ret1 || ret2 );
      }
    case SM_CSR:
      switch ( f ) {
        case INVALID:
          return -7;
        case DROW:
          ret1 = _csr2coo( m );
          ret2 = _coo2drow( m );
          return ( ret1 || ret2 );
        case DCOL:
          ret1 = _csr2coo( m );
          ret2 = _coo2drow( m );
          ret3 = _drow2dcol( m );
          return ( ret1 || ret2 || ret3 );
        case SM_COO:
          ret1 = _csr2coo( m );
          return ret1;
        case SM_CSC:
          ret1 = _csr2coo( m );
          ret2 = _coo2csc( m );
          return ( ret1 || ret2 );
        case SM_CSR:
          return 0; // nothing to do
      }
  }
  assert( 0 ); // shouldn't be able to get here due to returns
  return -8;
}

// swap upper-to-lower triangular and vice-versa
static inline void _symmetry_swap( matrix_t* m );
static inline void _symmetry_swap( matrix_t* m ) {
  assert( m->format == SM_COO );
  assert( m->sym = SM_SYMMETRIC );
  assert( m->location != BOTH );

  // swap ii and jj (rows and column indices)
  // an item at (1,2) is moved to (2,1) -- the other side of the triangle
  void *const t = m->ii;
  m->ii = m->jj;
  m->jj = t;

  // update matrix info
  if ( m->location == UPPER_TRIANGULAR )
    m->location = LOWER_TRIANGULAR;
  else // LOWER_TRIANGULAR
    m->location = UPPER_TRIANGULAR;
}


// duplicate data
// returns non-zero on malloc failure
static inline int _symmetry_both( matrix_t* m );
static inline int _symmetry_both( matrix_t* m ) {
  assert( m->format == SM_COO );
  assert( m->sym = SM_SYMMETRIC );
  assert( m->location != BOTH );

  const size_t dwidth = _data_width( m->data_type );
  const size_t nz_old = m->nz;
  if ( _realloc_arrays( m, nz_old*2 ) != 0 )
    return -1;

  assert( m->nz == nz_old*2 );

  // duplicate data
  // append the converse indices, duplicate the data
  // reverse ii and jj in the copied portions, appended to the old data
  memcpy( m->ii + nz_old, m->jj, nz_old*sizeof( unsigned int ) );
  memcpy( m->jj + nz_old, m->ii, nz_old*sizeof( unsigned int ) );
  memcpy( m->dd + nz_old * dwidth, m->dd, nz_old*dwidth );


  // update ptrs
  m->location = BOTH;
  // TODO clean this up!
  // TODO -- better to strip out entries on the diagonal as we go! (doing it the current way means copying the data multiple times
  // but now we might have duplicate entries if there were any on the diagonal
  // so remove anything on the diagonal that occurs in the second set of arrays
  int del = 0;
  int i;
  for ( i = nz_old;i < m->nz; i++ ) {
    if ( m->ii[i] == m->jj[i] ) { // on the diagonal
      del++; // deleting the i-th entry
      if ( i + 1 != m->nz ) { // unless its the last entry
        memcpy( m->ii + i, m->ii + ( i + 1 ), ( m->nz - i - 1 )*sizeof( unsigned int ) );
        memcpy( m->jj + i, m->jj + ( i + 1 ), ( m->nz - i - 1 )*sizeof( unsigned int ) );
        memcpy( m->dd + i * dwidth, m->dd + ( i + 1 ) * dwidth, ( m->nz - i - 1 )*dwidth );
      }
    }
  }
  _realloc_arrays( m, m->nz - del );
  // ignore return value .. we're shrinking the arrays and if the realloc fails we can continue on
  return 0;
}


static inline int _realloc_arrays( matrix_t* m, size_t nz ) {
  if ( nz == m->nz )
    return 0;

  const size_t dwidth = _data_width( m->data_type );
  // resize arrays
  unsigned int* ii_new = realloc( m->ii, nz * sizeof( unsigned int ) );
  unsigned int* jj_new = realloc( m->jj, nz * sizeof( unsigned int ) );
  void* dd_new         = realloc( m->dd, nz * dwidth );

  // udpate ptrs
  if (( ii_new != NULL ) && ( jj_new != NULL ) && ( dd_new != NULL ) ) {
    m->ii = ii_new;
    m->jj = jj_new;
    m->dd = dd_new;
    m->nz = nz;
    return 0;
  }
  else {
    // if the update shrank the arrays, its safe to ignore the failed realloc
    if ( m->nz > nz )
      m->nz = nz;
    return -1;
  }
}


// search for duplicate matrix entries
// combine the duplicates by adding
static inline void _coo_merge_duplicate_entries( matrix_t* m );
static inline void _coo_merge_duplicate_entries( matrix_t* m ) {
  assert( m->format == SM_COO );
  const size_t dwidth = _data_width( m->data_type );
  int del = 0;
  // expensive: this is O(log(n))?
  int i, j;
  assert( m->data_type == REAL_DOUBLE ); // TODO other data types
  double* d = ( double* ) m->dd;
  for ( i = 0;i < m->nz; i++ ) {
    for ( j = i + 1;i < m->nz; j++ ) {
      if (( m->ii[i] == m->ii[j] ) && ( m->jj[i] == m->jj[j] ) ) {
        // combine the duplicate entries (addition)
        d[i] = d[i] + d[j];

        del++; // deleting the j-th entry
        if ( j + 1 != m->nz ) { // unless its the last entry
          memcpy( m->ii + j, m->ii + ( j + 1 ), ( m->nz - j - 1 )*sizeof( unsigned int ) );
          memcpy( m->jj + j, m->jj + ( j + 1 ), ( m->nz - j - 1 )*sizeof( unsigned int ) );
          memcpy( m->dd + j * dwidth, m->dd + ( j + 1 ) * dwidth, ( m->nz - j - 1 )*dwidth );
        }
      }
    }
  }
  _realloc_arrays( m, m->nz - del );
}


// convert from symmetry BOTH -> LOWER_TRIANGULAR
// note that realloc might fail but we can carry on: no failure cases
static inline void _symmetry_lower( matrix_t* m );
static inline void _symmetry_lower( matrix_t* m ) {
  assert( m->format == SM_COO );
  assert( m->sym == SM_SYMMETRIC );
  assert( m->location == BOTH );

  // remove redundant entries in the upper triangle
  const size_t dwidth = _data_width( m->data_type );
  size_t nz = 0;
  int i;
  for ( i = 0;i < m->nz; i++ ) {
    if ( m->ii[i] >= m->jj[i] ) { // if row is >= column (lower triangle), keep this entry
      if ( i != nz ) { // don't need to copy if its staying at the current location
        memcpy( m->ii + nz , m->ii + i, sizeof( unsigned int ) );
        memcpy( m->jj + nz , m->jj + i, sizeof( unsigned int ) );
        memcpy( m->dd + nz * dwidth, m->dd + i * dwidth, dwidth );
      }
      nz++;
    }
  }
  _realloc_arrays( m, nz );
  m->location = LOWER_TRIANGULAR;
}

int convert_matrix_symmetry( matrix_t* m, enum matrix_symmetric_storage_t loc ) {
  assert( m->sym == SM_SYMMETRIC );

  // short circuit if no work to do
  //if(m->location == loc)
  //  return 0;

  int ret;
  const enum matrix_format_t old_format = m->format;

  // convert to COO format
  // TODO handle other formats directly (changing formats is expensive)?
  if (( ret = convert_matrix( m, SM_COO, m->base ) ) != 0 )
    return ret;

  int ret1 = 0;
  switch ( m->location ) {
    case BOTH:
      switch ( loc ) {
        case BOTH: // nothing to do
          break;
        case UPPER_TRIANGULAR:
          _symmetry_lower( m );
          _symmetry_swap( m );
          break;
        case LOWER_TRIANGULAR:
          _symmetry_lower( m );
          break;
      }
      break;
    case UPPER_TRIANGULAR:
      switch ( loc ) {
        case BOTH:
          ret1 = _symmetry_both( m );
          break;
        case UPPER_TRIANGULAR:
          break; // nothing to do
        case LOWER_TRIANGULAR:
          _symmetry_swap( m );
          break;
      }
      break;
    case LOWER_TRIANGULAR:
      switch ( loc ) {
        case BOTH:
          ret1 = _symmetry_both( m );
          break;
        case UPPER_TRIANGULAR:
          _symmetry_swap( m );
          break;
        case LOWER_TRIANGULAR:
          break; // nothing to do
      }
      break;
  }

  // return to original format
  ret = convert_matrix( m, old_format, m->base );
  return ( ret || ret1 );
}

// non-zero indicates failure (malloc)
int detect_matrix_symmetry( matrix_t* m ) {
  // if the matrix isn't currently unsymmetric, then
  // its already been decided that the matrix is symmetric
  if ( m->sym != SM_UNSYMMETRIC )
    return 0;

  const enum matrix_format_t old_format = m->format;
  int ret;

  // convert to COO format
  // TODO handle other formats directly (changing formats is expensive)
  if (( ret = convert_matrix( m, SM_COO, m->base ) ) != 0 )
    return ret;

  const double tol = 1e-15; // TODO use limits.h machine precision -- different tolerances for patterns
  int pattern_is_symmetric = 1;
  int i, j;
  // TODO handle other data formats
  assert( m->data_type == REAL_DOUBLE );
  double * const dd = m->dd;
  // TODO if entries were sorted this could be made more efficient
  // TODO need to ensure there are no entries for the same location? or does this matter? (probably, since these should be added together first)
  // TODO need to check for hermitian, skew symetric too
  // TODO only need to check half the entries?
  // loop over all entries, stop if we find a single unsymmetric entry
  // yuck: O(n^2) currently
  // TODO _coo_merge_duplicate_entries(m);
  for ( i = 0;i < m->nz && pattern_is_symmetric; i++ ) {
    // ignore entries on the diagonal
    if ( m->ii[i] != m->jj[i] ) {
      int match = 0;
      // then look for a corresponding entry on the other side of the matrix
      for ( j = 0; j < m->nz && !match; j++ ) { // stop when we find a match
        if (( i != j ) && // not the same entry
            ( m->ii[i] == m->jj[j] ) && // indices are flipped (diagonal symmetry
            ( m->jj[i] == m->ii[j] ) ) {
          // found a matching entry: potential match
          if (( dd[i] < dd[j] + tol ) && // and data matches within machine precision
              ( dd[i] > dd[j] - tol ) ) {
            match = 1;
          }
          else { // data doesn't match.. we can quit now!
            pattern_is_symmetric = 0; // done searching!
            match = 1; // just so we can break out of here
          }
        }
      }
      if ( !match )
        pattern_is_symmetric = 0;
    }
  }

  if ( pattern_is_symmetric ) {
    m->sym = SM_SYMMETRIC;
    m->location = BOTH;
  }

  // return to original format
  ret = convert_matrix( m, old_format, m->base );
  return ret;
}


// test matrix
// returns zero on pass
//   -1: bad size
//   -2: bad ptrs
//   -3: CSR/CSC wrong size in ptr array
//   -4: CSR/CSC bad first ptr entry (expect 0)
int validate_matrix( matrix_t* m ) {
  assert( m != NULL );

  // must be invalid if zero size matrix
  if ((( m->m == 0 ) || ( m->n == 0 ) ) && ( m->format != INVALID ) ) // empty matrix
    return -1;

  if ((( m->format == DROW ) || ( m->format == DCOL ) ) && ( m->nz != m->m * m->n ) )
    return -1;

  // empty matrices don't need to be validated further
  // but if its empty it aught to match the descriptor: no ptrs, no size
  if ( m->format == INVALID ) { // empty matrix
    if (( m->m != 0 ) || ( m->n != 0 ) || ( m->nz != 0 ) )
      return -1;
    else if (( m->ii != NULL ) || ( m->jj != NULL ) || ( m->dd != NULL ) )
      return -2;
    else
      return 0;
  }

  // if there is no data, ptrs should be null
  if (( m->nz == 0 ) && (( m->ii != NULL ) || ( m->jj != NULL ) || ( m->dd != NULL ) ) )
    return -2;

  // pattern type matrices can't hold data, only indices
  if (( m->data_type == SM_PATTERN ) && ( m->dd != NULL ) )
    return -2;

  // if dense, ii, jj = NULL
  if ((( m->format == DROW ) || ( m->format == DCOL ) ) && (( m->ii != NULL ) || ( m->jj != NULL ) ) )
    return -2;

  // if CSC/CSR, then nz must match expected value in first & last element of m->ii/jj
  if (( m->format == SM_CSR ) && (( m->ii )[m->m] != m->nz ) )
    return -3;
  if (( m->format == SM_CSC ) && (( m->jj )[m->n] != m->nz ) )
    return -3;
  if (( m->format == SM_CSR ) && (( m->ii )[0] != 0 ) )
    return -4;
  if (( m->format == SM_CSC ) && (( m->jj )[0] != 0 ) )
    return -4;

  // TODO for CSC/CSR/COO check for duplicate entrys (should be summed)
  // TODO if symmetric check there aren't any extra values in the other triangle
  // TODO check that no matrix indices are outside the declared matrix rows/cols

  // all is okay
  return 0;
}

void printf_matrix( char const *const pre, matrix_t* m ) {
  assert( m != NULL );

  matrix_t *const  c = m; // TODO copy_matrix( m );
  assert( c != NULL );
  convert_matrix( c, SM_COO, FIRST_INDEX_ZERO ); // TODO return value?

  int i;
  assert( c->data_type == REAL_DOUBLE );
  double* v = c->dd;
  if ( c->n == 1 ) { // if its a vector
    for ( i = 0; i < c->nz; i++ )
      printf( "%s(%i)=%.2f\n", pre, c->ii[i], v[i] );
  }
  else { // its a matrix
    for ( i = 0; i < c->nz; i++ )
      printf( "%s(%i,%i)=%.2f\n", pre, c->ii[i], c->jj[i], v[i] );
  }

  //TODO free_matrix( c );
}
