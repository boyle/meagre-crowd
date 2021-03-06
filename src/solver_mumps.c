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
#include "solver_mumps.h"
#include "solvers.h"
#include "matrix.h"
#include <dmumps_c.h>

#include <stdlib.h> // malloc, free
#include <string.h> // memcpy
#include <assert.h>


// mumps job controls
#define JOB_INIT -1
#define JOB_END -2

#define MUMPS_USE_COMM_WORLD -987654


// the functions called from the "solver" wrapper and defined in "solver_lookup.h"
// TODO probably need to see an example matrix so we can choose appropriate options, etc for solver
void solver_init_mumps( solver_state_t* s ) {
  // initialize MUMPS instance
  DMUMPS_STRUC_C* id = calloc( 1, sizeof( DMUMPS_STRUC_C ) ); // initialize to zero
  assert( id != NULL ); // calloc failure
  s->specific = id; // stored for future use

  id->job = JOB_INIT;
  id->par = 1; // host involved in factorization/solve
  id->sym = 0; // 0: general, 1: sym pos def, 2: sym (note: no hermitian support) // TODO support other matrix types
  // Note: if set to symmetric and matrix ISN'T, the redundant entries will be *summed*
  // TODO could convert the C communicator instead of using the fortran one (see MUMPS doc)
  id->comm_fortran = MUMPS_USE_COMM_WORLD;
#define INFOG(I) infog[(I)-1] // macro s.t. indices match documentation
  dmumps_c( id );
  assert( id->INFOG( 1 ) == 0 ); // check it worked
  // clears the rest of the unused values

  // set debug verbosity
#define ICNTL(I) icntl[(I)-1] // macro s.t. indices match documentation
  // No outputs
  if ( s->verbosity < 3 ) { // no debug
    id->ICNTL( 1 ) = -1;
    id->ICNTL( 2 ) = -1;
    id->ICNTL( 3 ) = -1;
    id->ICNTL( 4 ) = 0;
  }
  else { // debug
    id->ICNTL( 1 ) = 6; // err output stream
    id->ICNTL( 2 ) = 6; // warn/info output stream
    id->ICNTL( 3 ) = 6; // global output stream
    id->ICNTL( 4 ) = 4; // debug level 0:none, 1: err, 2: warn/stats 3:diagnostics, 4:parameters
  }
}

void solver_finalize_mumps( solver_state_t* s ) {
  assert( s != NULL );
  DMUMPS_STRUC_C* const id = s->specific;
  assert( id != NULL );

  id->job = JOB_END;
  dmumps_c( id ); // Terminate instance
  assert( id->INFOG( 1 ) == 0 ); // check it worked

  free( id->rhs );
  free( id );

  s->specific = NULL;
}

void solver_analyze_mumps( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  if ( s->mpi_rank == 0 )
    assert( A != NULL );
  else
    assert( A == NULL );

  DMUMPS_STRUC_C* const id = s->specific;

  // load A's pattern for analysis (COO format)
  if ( s->mpi_rank == 0 ) { // only rank=0 needs the initial data

    // TODO could do this smarter by giving MUMPS the symmetric matrix to solve,
    // instead we're just going straight for the unsymmetric solver
    assert( A->sym == SM_UNSYMMETRIC );
    assert( A->format == SM_COO );
    assert( A->base = FIRST_INDEX_ONE );

    assert( A->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO
    assert( A->m == A->n ); // square matrices only?


    // TODO really we should just copy this to be CORRECT/TYPESAFE (not worth being clever...)

    // mumps: irn=row indices, jcn=column idices, a=values, rhs=right-hand side, n = matrix order (on-a-side?) nz=non-zeros?
    id->n   = A->m; // A.m: rows, A.n: columns
    id->nz  = A->nz; // non-zeros
    id->irn = ( int* ) A->ii; // row    indices
    id->jcn = ( int* ) A->jj; // column indices
  }

  // Call the MUMPS package.
  // ICNTL(22) != 0: out-of-core
  //   requires OOC_TMPDIR, OOC_PREFIX -> tmp location/prefix
  // ICNTL(14): memory relaxation
  // ICNTL(5) = ICNTL(18) = 0; // centralized, assembled matrix load
  // set ICNTL(7)=1 for external ordering, requires id->PERM_IN
  //   0: AMD
  //   1: user provided, id->PERM_IN
  //   2: AMF
  //   3: SCOTCH
  //   4: PORD
  //   5: METIS
  //   6: QAMD (good when dense)
  //   7: auto
  //   for schur complement, only 0,1,5,7
  //   for elemental matrices, only 0,1,5,7
  //     both: 0,1
  //   INFOG(7) shows what was selected
  //   IGNORED if ICNTL(28)=2
  // ICNTL(8): scaling strategy (computed ICNTL(6) ICNTL(12) duing analysis)
  //   -2: computed during analysis
  //   -1: user provided COLSCA, ROWSCA
  //   0: none
  //   1: diagonal, during factorization
  //   2: row/column scaling during factorization
  //   3: column scaling during factorization
  //   4: row/column based on inf. norms
  //   5: another way (column)
  //   6: another way (row/col)
  //   7: iter. row/col
  //   8: similar to 8 but slower, more rigorous
  //   77: auto (analysis only), chooses ICNTL(8)
  //   INFOG(33) shows what was selected
  // ICNTL(28)=1 sequential analysis, 2: parallel analysis
  // set ICNTL(6)=5,6 for scaling (needs values in id->A)
  //   0: no column permutations calculated
  //   1: maximise # of diagonals
  //   2: maximise smallest diagonal entry
  //   3: same as 2, different performance
  //   4: maximize trace (sum of diagonal)
  //   5: maximize diagonal product, scale if ICNTL(8)=-2 or 77
  //   6: same as 5, different performance
  //   7: automatic scaling, choose appropriate method (default)
  //      INFOG(23) shows what was choosen
  // required, only on host:
  //   id->N, NZ, IRN, JCN (assembled matrix, ICNTL(5)=0 & ICNTL(18) != 0) or
  //     (where ICNTL(18)=1 or 2 -- distributed, IRN, JCN not required)???
  //     (where ICNTL(18)=3 -- distributed, N on host, and NZ_loc,
  //                           IRN_loc, JCN_loc on slaves)???
  //   id->N, NELT, ELTPTR, ELTVAR (elemental matrix, ICNTL(5)=1)
  //   ICNTL(19)=1,2,3, ICNTL(26) schur complement w/ reduced or condensed rhs
  //     requries id->SIZE_SCHUR, LISTVAR_SCHUR
  //     ICNTL(19)=1: centralized schur, 2: distributed lower schur (sym), or
  //               3: distributed complete schur
  //       ICNTL(19)=2,3 requires id->NPROW, NPCOL, MBLOCK, NBLOCK - processing conf
  //       ... sets SCHUR_MLOC, SCHUR_NLOC
  // id->WRITE_PROBLEM: store distributed in matrix market format
#define JOB_ANALYSE 1
  id->job = JOB_ANALYSE;
  dmumps_c( id );
  if ( id->INFOG( 1 ) != 0 ) fprintf( stderr, "warning: analysis failed\n" );
  assert( id->INFOG( 1 ) == 0 ); // check it worked

  // available info:
  // INFO(15)/INFOG(16/17): min/max/sum-over-all-cpus mem requried [in megabytes]
  // INFO(17): min mem for out-of-core (max, sum in INFOG(26,27))
  //   set ICNTL(23) for explicit max mem, per-proc [MB]
  //   set ICNTL(14) to limit mem increases
}

void solver_factorize_mumps( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  if ( s->mpi_rank == 0 )
    assert( A != NULL );
  else
    assert( A == NULL );
  DMUMPS_STRUC_C* const id = s->specific;

  // load A's *data* for factorization
  // Note: pattern must have remained the same
  if ( s->mpi_rank == 0 )
    id->a   = A->dd;

  // requires id->A if ICNTL(5)=0 (assembled matrix)
  // requires id->A_ELT if ICNTL(5)=1 (elemental matrix)
  // if (ICNTL(5)=0 && ICNTL(18)!=0) (assembled matrix, distributed load) ???i
  //   requires A_loc on slaves
  //   requires NZ_loc, IRN_loc, JCN_loc if ICNTL(18)=1 or 2 ???
  //       -- already passed in if ICNTL(18)=3
  // ICNTL(8)!=0 --> user supplied row/column scaling
  //   requires id->COLSCA, ROWSCA
  // ICNTL(19)=2,3 requires SCHUR_LLD, SCHUR
#define JOB_FACTORIZE 2
  id->job = JOB_FACTORIZE;
  dmumps_c( id );
  assert( id->INFOG( 1 ) == 0 ); // check it worked
}

void solver_evaluate_mumps( solver_state_t* s, matrix_t* b, matrix_t* x ) {
  assert( s != NULL );
  DMUMPS_STRUC_C* const id = s->specific;

  if ( s->mpi_rank == 0 ) {
    assert( x != NULL );
    assert( b != NULL );
    assert( b->data_type == REAL_DOUBLE );
    assert( b->m == id->n ); // rows of b match rows of A
    assert( (b->format == DCOL) || (b->format == SM_CSC) );

    id->lrhs = b->m; // rows of b
    id->nrhs = b->n; // columns of b
    id->rhs = malloc( id->n * id->nrhs * sizeof( double ) );
    assert( id->rhs != NULL ); // malloc failure
    if(b->format == SM_CSC) { // sparse RHS
      assert(b->base == FIRST_INDEX_ONE);
      id->ICNTL(20) = 1;
      id->nz_rhs = b->nz;             // non-zeros
      id->rhs_sparse = b->dd;         // data
      // TODO check the cast is safe: rows < max_int, nz < max_int, so we don't muck it up when we drop the signed-ness
      id->irhs_sparse = (int*) b->ii; // row indices
      id->irhs_ptr = (int*) b->jj;    // column ptrs
    }
    else { // dense RHS
      // need to copy 'b' in since it gets destroyed
      // TODO unless x == b && b != CSC
      memcpy( id->rhs, b->dd, id->n * id->nrhs * sizeof( double ) );
    }
  }
  else {
    assert( b == NULL );
    assert( x == NULL );
  }

  // solve Ax=b, AX=B OR A^t x=b, A^t X=B
  //   ICNTL(9)=1: A (default), 0: A^t
  // OR compute "null-space basis" if "null pivot row detection" was enabled
  //   ICNTL(24)=1 and INFOG(28)!=0
  // requires id->RHS
  // ICNTL(26)=0: solve internal problem, 1,2: reduced rhs on schur variables
  //   requires LREDRHS REDRHS
  // ICNTL(20)=0: centralized dense rhs, 1: sparse rhs
  //   1: NZ_RHS, NRHS, RHS_SPARSE, IRHS_SPARSE, IRHS_PTR
  // ICNTL(21)=0: centralized dense soln, 1: distributed soln
  //   1: SOL_LOC, LSOL_LOC, ISOL_LOC
  // ICNTL(10): iterative refinement
  // ICNTL(11): error analysis (expensive?)
  //   RINFOG(4): A's inf norm
  //   RINFOG(5): soln residual
  //   RINFOG(6): scaled soln residual
  //   RINFOG(7/8): backward error estimate
  //   RINFOG(9): soln err est
  //   RINFO(10/11): condition numbers
  // id->NRHS is number of rhs (optional), >1 disables itr refinement, err analysis
  // id->LRHS >= NRHS, =leading dimension of RHS (optional)
  //
#define JOB_SOLVE 3
  id->job = JOB_SOLVE;
  dmumps_c( id );
  assert( id->INFOG( 1 ) == 0 ); // check it worked

  // put the answer in a nice formatted bundle
  if ( s->mpi_rank == 0 ) {
    clear_matrix( x );
    x->m = id->n;
    x->n = id->nrhs;
    x->nz = x->m * x->n;
    x->format = DCOL;
    x->data_type = REAL_DOUBLE;
    // we can recycle this pointer:
    //  1. we allocated it just prior to the solve call, so it doesn't belong to 'b'
    //  2. the data was copied from 'b' then overwritten in the solve stage so no need to copy again
    x->dd = id->rhs;
    id->rhs = NULL; // transfer pointer ownership to x
  }
}
