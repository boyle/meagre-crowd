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
#include <stdio.h>
#include "umfpack.h"
int n = 5 ;
int Ap [ ] = {0, 2, 5, 9, 10, 12} ;
int Ai [ ] = { 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
double b [ ] = {8., 45., -3., 3., 19.} ;
double x [5] ;
double r [ ] = {1., 2., 3., 4., 5.};
double p = 1e-14;
int main( void ) {
  double *null = ( double * ) NULL ;
  int i ;
  void *Symbolic, *Numeric ;
  ( void ) umfpack_di_symbolic( n, n, Ap, Ai, Ax, &Symbolic, null, null ) ;
  ( void ) umfpack_di_numeric( Ap, Ai, Ax, Symbolic, &Numeric, null, null ) ;
  umfpack_di_free_symbolic( &Symbolic ) ;
  ( void ) umfpack_di_solve( UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null ) ;
  umfpack_di_free_numeric( &Numeric ) ;
  for ( i = 0 ; i < n ; i++ ) printf( "x [%d] = %g\n", i, x [i] ) ;
  for ( i = 0 ; i < n ; i++ ) if (( x[i] > r[i] + p ) || ( x[i] < r[i] - p ) ) {
      printf( "FAIL\n" );
      return 1;
    }
  printf( "PASS\n" );
  return ( 0 ) ;
}
