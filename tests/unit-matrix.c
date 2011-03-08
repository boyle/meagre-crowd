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
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "matrix.h"


// prototypes
struct enum2format_t {
  const enum matrix_format_t f;
  const char* s;
};

// how many matrix_format_t entries are there
#define matrix_format_t__MAX 6
static const struct enum2format_t enum2format[] = { {INVALID, "invalid"},
  {DROW,    "dense (rows)"},
  {DCOL,    "dense (columns)"},
  {SM_COO,  "COO"},
  {SM_CSC,  "CSC"},
  {SM_CSR,  "CSR"}
};

struct enum2base_t {
  const enum matrix_base_t b;
  const char* s;
};

const struct enum2base_t enum2base[] = { {FIRST_INDEX_ZERO, "zero"},
  {FIRST_INDEX_ONE, "one"}
};

void print_matrix( matrix_t* a );
void test_formats( matrix_t* a );
void test_symmetry( matrix_t* a );
void test_copy( matrix_t* a );
void test_basic();
void build_test_matrix( matrix_t** m, int i );


void print_matrix( matrix_t* a ) {
  printf( "    %zux%zu (nz:%zu) %s%s%s%s %s %s%d\n",
          a->m, a->n, a->nz,
          enum2format[a->format].s,
          ( a->dd == NULL ) ? " !dd" : "",
          ( a->ii == NULL ) ? " !ii" : "",
          ( a->jj == NULL ) ? " !jj" : "",
          ( a->base == FIRST_INDEX_ZERO ) ? "base0" : "base1",
          ( a->sym == SM_UNSYMMETRIC ) ? "unsym" : "sym",
          ( a->location ) );
  if ( a->format == SM_COO && a->data_type == REAL_DOUBLE ) {
    printf( "      dd:" );
    int i;
    for ( i = 0; i < a->nz; i++ )
      printf( " %e", ( double )(( double* )( a->dd ) )[i] );
    printf( "\n" );
    printf( "      ii:" );
    for ( i = 0; i < a->nz; i++ )
      printf( " %d", ( unsigned int )(( unsigned int* )( a->ii ) )[i] );
    printf( "\n" );
    printf( "      jj:" );
    for ( i = 0; i < a->nz; i++ )
      printf( " %d", ( unsigned int )(( unsigned int* )( a->jj ) )[i] );
    printf( "\n" );
  }
}

void test_symmetry( matrix_t* a ) {
  printf( "symmetry test\n" );
  assert( a->sym == SM_SYMMETRIC );
  matrix_t* b = copy_matrix( a );
  assert( b != NULL );
  assert( validate_matrix( a ) == 0 );
  assert( validate_matrix( b ) == 0 );
  print_matrix( b );

  // de-symmetrize the matrix and check it gets detected correctly
  assert( a->sym == SM_SYMMETRIC );
  assert( convert_matrix_symmetry( a, BOTH ) == 0 );
  a->sym = SM_UNSYMMETRIC;
  print_matrix( a );
  assert( detect_matrix_symmetry( a ) == 0 );
  print_matrix( a );
  printf( "cmp_matrix(b[orig],a(just-detected))=%d\n", cmp_matrix( b, a ) );
  assert( cmp_matrix( b, a ) == 0 );
  print_matrix( a );

  // try all conversions
  int i, j;
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0;j < 3; j++ ) {

      printf( "  %d -> %d\n", i, j );

      assert( convert_matrix_symmetry( a, i ) == 0 );
      assert( a->sym == SM_SYMMETRIC );
      assert( a->location == i );

      assert( convert_matrix_symmetry( a, j ) == 0 );
      assert( a->sym == SM_SYMMETRIC );
      assert( a->location == j );

      print_matrix( a );

      assert( cmp_matrix( b, a ) == 0 );
    }
  }

  free_matrix( b );
}

void test_copy( matrix_t* a ) {
  printf( "copy test\n" );
  int j;
  for ( j = 1; j < matrix_format_t__MAX; j++ ) { // storage formats, skipping INVALID
    matrix_t* b = copy_matrix( a );
    assert( b != NULL );
    assert( validate_matrix( a ) == 0 );
    assert( validate_matrix( b ) == 0 );

    int ret = convert_matrix( b, ( enum matrix_format_t ) j, FIRST_INDEX_ZERO );
    assert(ret == 0);

    matrix_t* c = copy_matrix( b );
    free_matrix( b );

    print_matrix( c );
    free_matrix( c );
  }
}

void test_formats( matrix_t* a ) {
  matrix_t* b = copy_matrix( a );
  assert( b != NULL );
  print_matrix( b );

  int i, j, k, l, ret;
  for ( i = 0;i < 2;i++ ) { // from two types of base: 0 or 1
    for ( j = 0; j < matrix_format_t__MAX; j++ ) { // storage formats
      enum matrix_format_t old_format = b->format;
      printf( "%s -> %s (base %s -> %s)\n", enum2format[b->format].s, enum2format[j].s, enum2base[b->base].s, enum2base[i].s );
      ret = convert_matrix( b, ( enum matrix_format_t ) j, ( enum matrix_base_t ) i );
      if ( j == 0 ) { // INVALID
        assert( ret != 0 );
        assert( b->format == old_format );
        b->format = INVALID; // force it to be invalid so we can test the conversions
      }
      else {
        assert( ret == 0 );
        assert( b->format == ( enum matrix_format_t ) j );
        if (( b->format == DROW ) || ( b->format == DCOL ) )
          assert( b->base == FIRST_INDEX_ZERO );
        else
          assert( b->base == ( enum matrix_base_t ) i );
        assert( cmp_matrix( a, b ) == 0 );
      }

      for ( k = 0;k < 2;k++ ) { // to base
        for ( l = 1; l < matrix_format_t__MAX; l++ ) { // to all other storage formats
          matrix_t* c = copy_matrix( b );
          assert( c != NULL );

          printf( "  %s -> %s (base %s -> %s)\n", enum2format[j].s, enum2format[l].s, enum2base[i].s, enum2base[k].s );
          ret = convert_matrix( c, ( enum matrix_format_t ) l, ( enum matrix_base_t ) k );
          if (( j == 0 ) || ( l == 0 ) ) { // INVALID -> x OR x -> INVALID
            assert( ret != 0 );
            assert( c->format == ( enum matrix_format_t ) j );
          }
          else {
            print_matrix( b );
            print_matrix( c );
            assert( ret == 0 );
            assert( c->format == ( enum matrix_format_t ) l );
            if (( c->format == DROW ) || ( c->format == DCOL ) )
              assert( c->base == FIRST_INDEX_ZERO );
            else
              assert( c->base == ( enum matrix_base_t ) k );
            printf( "    cmp_matrix(b,c) = %d\n", cmp_matrix( b, c ) );
            printf( "    cmp_matrix(c,b) = %d\n", cmp_matrix( c, b ) );
            assert( cmp_matrix( b, c ) == 0 );
            assert( cmp_matrix( c, b ) == 0 );
          }

          free_matrix( c );
        }
      }
      // return to old format
      if ( b->format == INVALID )
        b->format = old_format;
    }
  }
  free_matrix( b );
}
// TODO check that if it starts as invalid, it stays invalid type
// TODO check that if it isn't invalid, it fails when asking to become invalid type

void test_basic() {
  // try mucking around with an INVALID matrix
  matrix_t* a = malloc_matrix();
  assert( a != NULL );
  *a = ( matrix_t ) {
    0
  };
  assert( validate_matrix( a ) == 0 );

  matrix_t* b = copy_matrix( a );
  assert( b != NULL );
  assert( validate_matrix( b ) == 0 );
  assert( cmp_matrix( a, b ) == 0 );

  free_matrix( b );
  free_matrix( a );

  // create a couple of test matrices
  // TODO generate random matrices?
  matrix_t* c = NULL;

  // try converting between all types
  build_test_matrix( &c, 0 );
  test_formats( c );


  build_test_matrix( &c, 1 );
  test_symmetry( c );

  build_test_matrix( &c, 0 );
  test_copy( c );


  // TODO do some cmp_matrix's that are supposed to fail in different ways

  // TODO test behaviour of an emtpy matrix (nz=0) in all formats
  free_matrix( c );
}

void build_test_matrix( matrix_t** m, int i ) {
  assert( m != NULL );
  free_matrix( *m );
  *m = malloc_matrix();
  assert( *m != NULL );
  matrix_t *const c = *m;
  *c = ( matrix_t ) {
    0
  };

  switch ( i ) {
    case 0:
      // unsymmetric matrix
      c->m = 8;
      c->n = 9;
      c->nz = 4;
      c->base = FIRST_INDEX_ZERO;
      c->format = SM_COO;
      // TODO try other symmetries/locations
      c->data_type = REAL_DOUBLE; // TODO try other data types
      c->dd = calloc( c->nz, sizeof( double ) );
      c->ii = calloc( c->nz, sizeof( unsigned int ) );
      c->jj = calloc( c->nz, sizeof( unsigned int ) );
      assert( c->dd != NULL );
      assert( c->ii != NULL );
      assert( c->jj != NULL );
      {
        int i;
        double* d = c->dd;
        for ( i = 0;i < c->nz; i++ ) {
          d[i] = ( double ) i + 10.0;
          c->ii[i] = i * 2;
          c->jj[i] = i + 4;
        }
      }
      break;

    case 1:
      c->m = 4;
      c->n = 4;
      c->nz = 2;
      c->base = FIRST_INDEX_ZERO;
      c->format = SM_COO;
      c->data_type = REAL_DOUBLE;
      c->dd = calloc( c->nz, sizeof( double ) );
      c->ii = calloc( c->nz, sizeof( unsigned int ) );
      c->jj = calloc( c->nz, sizeof( unsigned int ) );
      assert( c->dd != NULL );
      assert( c->ii != NULL );
      assert( c->jj != NULL );
      double* d = c->dd;
      d[0] = 1.0;
      d[1] = 2.0;
      c->ii[0] = 1;
      c->jj[0] = 1;
      c->ii[0] = 2;
      c->jj[0] = 3;
      c->sym = SM_SYMMETRIC;
      c->location = UPPER_TRIANGULAR;
      break;

    default:
      assert( 0 ); // bad i
  }
  printf( "  generated test matrix %d\n", i );
  assert( validate_matrix( c ) == 0 );
}

int main( int argc, char **argv ) {
  test_basic();
  return 0;
}
