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
/* Example program using the C interface to the
* double precision version of MUMPS, dmumps_c.
* We solve the system A x = RHS with
* A = diag(1 2) and RHS = [1 4]ˆT
* Solution is [1 2]ˆT */
#include <stdio.h>
#include <assert.h>
#include "mpi.h"
#include "dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
int main( int argc, char ** argv ) {
  DMUMPS_STRUC_C id;
  int n = 2;
  int nz = 2;
  int irn[] = {1, 2};
  int jcn[] = {1, 2};
  double a[2];
  double rhs[2];
  int myid, ierr;

  ierr = MPI_Init( &argc, &argv );
  assert(ierr == 0);
  ierr = MPI_Comm_rank( MPI_COMM_WORLD, &myid );
  assert(ierr == 0);
  /* Define A and rhs */
  rhs[0] = 1.0;
  rhs[1] = 4.0;
  a[0] = 1.0;
  a[1] = 2.0;
  /* Initialize a MUMPS instance. Use MPI_COMM_WORLD. */
  id.job = JOB_INIT;
  id.par = 1;
  id.sym = 0;
  id.comm_fortran = USE_COMM_WORLD;
  dmumps_c( &id );
  /* Define the problem on the host */
  if ( myid == 0 ) {
    id.n = n;
    id.nz = nz;
    id.irn = irn;
    id.jcn = jcn;
    id.a = a;
    id.rhs = rhs;
  }
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
  /* No outputs */
  id.ICNTL( 1 ) = -1;
  id.ICNTL( 2 ) = -1;
  id.ICNTL( 3 ) = -1;
  id.ICNTL( 4 ) = 0;
  /* Call the MUMPS package. */
  id.job = 6;
  dmumps_c( &id );
  id.job = JOB_END;
  dmumps_c( &id ); /* Terminate instance */
  if ( myid == 0 ) {
    printf( "Solution is ( %.2f; %.2f ): %s\n", rhs[0], rhs[1], (( rhs[0] == a[0] ) && ( rhs[1] == a[1] ) ) ? "PASS" : "FAIL" );
  }
  ierr = MPI_Finalize();
  return 0;
}
