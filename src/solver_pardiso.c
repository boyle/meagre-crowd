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
#include "solver_pardiso.h"
#include "solvers.h"
#include "matrix.h"
//#include <pardiso.h>

#include <stdlib.h> // malloc, free
#include <string.h> // memcpy
#include <assert.h>
#include <stdint.h> // int64_t
#include <unistd.h> // dup -> redirect stdout temporarily

// NOTE: pardiso doesn't include a header file for its library, use these prototypes instead
void pardisoinit( void   * PT, int    * MTYPE,   int * SOLVER, int * IPARM, double * DPARM, int *ERROR );
void pardiso( void   * PT, int    * MAXFCT,   int * MNUM, int * MTYPE,    int * PHASE, int * N,
              double * A, int    * IA,    int * JA, int * PERM,   int * NRHS, int * IPARM,
              int * MSGLVL, double * B, double * X, int * ERROR, double * DPARM );
void pardiso_chkmatrix( int *, int *, double *, int *, int *, int * );
void pardiso_chkvec( int *, int *, double *, int * );
void pardiso_printstats( int *, int *, double *, int *, int *, int *,
                         double *, int * );

enum pardiso_mtype { REAL_PATTERN_SYM = 1, REAL_SYM_POSDEF = 2, REAL_SYM_INDEF = -2,
                     COMPLEX_PATTERN_SYM = 3, COMPLEX_HERMITIAN_POSDEF = 4,
                     COMPLEX_HERMITIAN_INDEF = -4,
                     COMPLEX_SYM = 6,
                     REAL_UNSYM = 11,
                     COMPLEX_UNSYM = 13
                   };
enum pardiso_solver { SPARSE_DIRECT_SOLVER = 0, MULTIRECURSIVE_ITERATIVE_SOLVER = 1};

enum pardiso_solve_error { NO_ERROR = 0, INCONSISTENT_INPUT = -1,
                           MEM = -2, REORDERING = -3,
                           ZERO_PIVOT_NUMERICAL_FACT_OR_ITR_REFINEMENT = -4,
                           INTERNAL_ERR = -5, PREORDER_FAILED = -6, DIAG_MATRIX_PROB = -7,
                           INT32_OVERFLOW = -8,
                           // pardiso.lic problems
                           NO_LICENSE = -10, EXPIRED_LICENSE = -11, LICENSE_WRONG_HOST_USER = -12,
                           // Krylov subspace issues
                           MAX_ITER = -100, NO_CONVERGENCE = -101, // (w/in 25 iterations)
                           ITER_ERR = -102, ITER_BREAKDOWN = -103
                         };

#define IPARM(I) iparm[(I)-1]
#define DPARM(I) dparm[(I)-1]
typedef struct {
  int Arows;
  int Acols;
  int* Aii; // these are pointers to existing data (DON'T free)
  int* Ajj;
  double* Add;

  // Paradiso config
  int64_t *PT; // int *64 (Paradiso private ptr) -- 64b version expects 64b ints
  int MTYPE;
  int *iparm; // int * 64 - config
  double *dparm; // double * 64 - config
} solve_system_pardiso_t;

void solver_init_pardiso( solver_state_t* s ) {
  assert( s != NULL );
  solve_system_pardiso_t * const p = calloc( 1, sizeof( solve_system_pardiso_t ) );
  assert( p != NULL );
  s->specific = p;

  p->PT = malloc( 64 * sizeof( int64_t ) );
  assert( p->PT != NULL );
  p->MTYPE = REAL_UNSYM; // default to real/unsym
  int SOLVER = 0; // default to direct solver
  p->iparm = calloc( 64, sizeof( int ) );
  p->dparm = NULL; // calloc( 64, sizeof( double ) ); ONLY required for iterative solver!

  // Note: p->DPARM is used for the iterative CG solver.. not supported (yet? TODO)

  // Note: other values unused
  p->IPARM( 1 ) = 0; // use default values (2, 4..64), 3 MUST be set .. this will overwrite the following config with defaults (it mostly here for documentation)
  p->IPARM( 2 ) = 2; // fill-in reducing ordering (0: min-degree, 2: METIS)
  p->IPARM( 3 ) = 1; // number of processors: must match OMP_NUM_THREADS TODO  -- NOTE this is an *upper-limit* on the number of processors...
  p->IPARM( 4 ) = 0; // LU preconditioned CGS (10*L+K) where K=1:CGS,2:CG L=10^-L stopping threshold
  p->IPARM( 5 ) = 0; // user permutation PERM
  p->IPARM( 6 ) = 0; // overwrite b with x "write solution on x"
  // p->IPARM( 7 ) is the number of iterative refinement steps executed
  p->IPARM( 8 ) = 0; // iterative refinement steps requested (automatically does 2, if perturbed pivots are used) (<0 uses real*16 or complex*32 residual if supported (not gfortran))
  // p->IPARM( 9 ) is unused, must be 0
  // pivot perturbation (MTYPE=11,13,-2,-4,6=COMPLEX/REAL_UNSYM, REAL_SYM_INDEF, COMPLEX_HERMETIAN_INDEF, COMPLEX_SYM)
  // ( 10^-x * ||A||_inf where A is scaled and permuted) theshold
  if (( p->MTYPE == REAL_UNSYM ) || ( p->MTYPE == COMPLEX_UNSYM ) )
    p->IPARM( 10 ) = 13;
  else
    p->IPARM( 10 ) = 8;
  p->IPARM( 11 ) = 1; // scaling vectors (MTYPE=11,13, if(IPARM(13) then also -2,-4,6)
  // use IPARM(11) = IPARM(13) = 1 for highly indefinite matrices: requires numerical values of A in analysis
  p->IPARM( 12 ) = 0; // solve transpose matrix A^t x = b TODO support this
  // improved accuracy, 1: normal matchines, 2: advanced matchings (highly indefinite matrices)
  if (( p->MTYPE == REAL_UNSYM ) || ( p->MTYPE == COMPLEX_UNSYM ) )
    p->IPARM( 13 ) = 1;
  else
    p->IPARM( 13 ) = 0;
  // p->IPARM( 14 ) is number of perturbed pivots
  // p->IPARM( 15 ) peak analysis/factorization memory (kB)
  // p->IPARM( 16 ) permanent memory (allocated at the end of phase 1) (kB)
  // p->IPARM( 17 ) total double precision memory (kB)
  //   peak memory consumption is max( IPARM(15), IPARM(16) + IPARM(17) )
  p->IPARM( 18 ) = -1; // < 0 requests non-zero factor count, returned in p->IPARM(18)
  p->IPARM( 19 ) = 0; // <0 requests MFlops of factorization (10^6) (performance cost)
  // p->IPARM(20) CG diagnostics >0 = iterations, <0 -itr*10 - err, err=1 residual fluctuations, 2 slow convergence, 3 itr limit, 4 perturbed pivots resulted in itr refinement, 5 factorization too fast -> use IPARM(4)=0)
  p->IPARM( 21 ) = 1; // pivoting for sym indef matrices: 0=1x1, 1=1x1 and 2x2 Bunch-Kaufman (MTYPE=-2,-4,6)
  // p->IPARM(22) inertia: +ve eigenvalues
  // p->IPARM(23) inertia: -ve eigenvalues
  p->IPARM( 24 ) = 1; // parallel factorization 0:old method, 1:new method: two-level scheduling
  p->IPARM( 25 ) = 1; // parallel forward/backward solve
  p->IPARM( 26 ) = 0; // split forward/backward solve 0: LU (normal), 1: fwd (L or U^t), 2: backward (U or L^t) (see IPARM(12) solving A^t)
  // p->IPARM( 27 ) is unused, must be 0
  p->IPARM( 28 ) = 0; // TODO use parallel METIS (MPI)
  p->IPARM( 29 ) = 0; // 1:32b or 0:64b factorization (sym indef & real unsym MTYPE=-2,11 only)
  p->IPARM( 30 ) = 80; // column size of supernodes
  p->IPARM( 31 ) = 0; // partial solve rhs -- PERM(i)=1 compute this component / this rhs is non-zero (only want certain components of solution -- i.e. EIT!!) (sparse RHS only) (required at all phases)
  p->IPARM( 32 ) = 0; // choose solver 0: direct, 1: iterative
  p->IPARM( 33 ) = 0; // determinant of real sym indef in 1=compute --> result in DPARM(33), IPARM(33) = ln(abs(DPARM(33)))
  p->IPARM( 34 ) = 0; // bit identical solution (only for METIS (IPARM(2)=2), IPARM(24)=1, IPARM(25)=1) valid for IPARM(3) (upper processor limit) being constant
  // p->IPARM( 35 ..50 ) are unused, must be 0
  p->IPARM( 51 ) = 0; // TODO use MPI factorize/solve
  p->IPARM( 52 ) = 1; // number of compute nodes (MPI) if IPARM(51) = 1 TODO
  // p->IPARM( 53 .. 65 ) are unused, must be 0

  // launch pardiso
  // TODO do we need to know the matrix type when we start here?? might need to move to analyze stage...
  int error_code = 0;

  // we redirect stdout to /dev/null around this call to skip the licensing warnings
  // Very non-portable code...
  // check:
  //  printf("stdout is being redirected to /dev/null\n");
  // Method 1:
  //  int old_stdout = dup(STDOUT_FILENO);
  //  assert(old_stdout > 0); // or an error occurred
  //  FILE *fp1 = freopen("/dev/null","w", stdout);
  //  assert(fp1 != NULL);
  // Method 2:
  FILE *o = stdout;
  stdout = fopen( "/dev/null", "w" );
  assert( stdout != NULL );
  // Method 3:
  //  int o = dup(fileno(stdout);
  //  assert(o != -1);
  //  freopen = fopen("/dev/null","w");
  //  assert(freopen != NULL);
  // check:
  //  printf("can't see this. RIGHT?\n");

  pardisoinit( p->PT, &( p->MTYPE ), &SOLVER, p->iparm, p->dparm, &error_code );

  // the trick though, is to return stdout to normal afterwards
  // Note: this definitely is different on windows!
  // Method 1:
  //  FILE *fp2 = fdopen(old_stdout, "w");
  //  assert(fp2 != NULL);
  //  int ret = fclose(stdout);
  //  assert(ret == 0);
  //  stdout = fp2; // glibc
  //  OR *stdout = *fp2; // Solaris and MacOS X
  //  ret = close(old_stdout);
  //  assert(ret == 0);
  // Method 2:
  int ret = fclose( stdout );
  assert( ret == 0 );
  stdout = o;
  // Method 3:
  //  int ret = dup2(o, fileno(stdout));
  //  assert(ret != -1);
  //  ret = close(o);
  //  assert(ret == 0);
  // check:
  //  printf("stdout has been restored!\n");



  // check for errors from paradiso
  if ( error_code != NO_ERROR )
    fprintf( stderr, "error: pardiso initialization code %d\n", error_code ); // TODO decode
  assert( error_code == NO_ERROR );
}

void solver_analyze_pardiso( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  assert( A != NULL );
  solve_system_pardiso_t* const p = s->specific;
  assert( p != NULL );

  // TODO could do this smarter by giving MUMPS the symmetric matrix to solve,
  // instead we're just going straight for the unsymmetric solver
  // TODO this logic should really be moved to the wrapper solver() functions using the capabilities masks
  if ( A->sym == SM_SYMMETRIC )
    convert_matrix_symmetry( A, BOTH );
  // TODO for SYMMETRIC matrices: pardiso wants all diagonal entries to exist in A, even if they are zero
  // TODO pardiso can handle pattern symmetric matrices that have unsymmetric data

  // TODO speial handling if the matrix has a symmetric pattern but unsymmetric data .. copy, delete data, set to pattern and test

  // prepare the matrix
  int ierr = convert_matrix( A, SM_CSR, FIRST_INDEX_ONE );
  assert( ierr == 0 );
  assert(( A->sym == SM_UNSYMMETRIC ) || ( A->sym == SM_SYMMETRIC ) );
  assert( A->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO
  assert( A->m == A->n ); // TODO can only handle square matrices at present???

  // Compressed Row Format
  assert( A->ii[0] == 1 );
  assert( A->ii[A->n] == A->nz + 1 );

  int MAXFCT = 1;
  int MNUM = 1; // use which of MAXFACT stored factorizations in solution phase
  p->MTYPE = REAL_UNSYM; // TODO other types // TODO check it stays consistent in factorize stage
  int PHASE = 11; // 10*i + j, do phases i->j, 1: analysis, 2; factorize, 3: solve w/ itr refinement, 4: terminate and release mem for MNUM, <0: release all mem
  // Note: needs data in A if IPARM(11) = 1 and IPARM(13)=1 or 2
  int N = A->m; // rows of A
  // TODO might need to save this for later?
  // expects A->dd to be "double" either real or complex
  // expects "int" for A->ii, A->jj
  int* const PERM = NULL; // if(IPARM(5)=1, user fill-reducing permutation (ordering) (size N) B = P A P^t, PERM(i) = row i of A to column PERM(i) of B
  int NRHS = 0; // columns in rhs
  int MSGLVL = 0; // TODO set based on s->verbosity
  double* const B = NULL; // N x NRHS, replaced w/ X if IPARM(6) = 1 -- solution phase
  double* const X = NULL; // N x NRHS solution if IPRAM(6) = 0 -- solution phase
  int error_code = 0;
  pardiso( p->PT, &MAXFCT, &MNUM, &( p->MTYPE ), &PHASE, &N,
           A->dd, ( int* ) A->ii, ( int* ) A->jj, PERM, &NRHS, p->iparm,
           &MSGLVL, B, X, &error_code, p->dparm );
  if ( error_code != NO_ERROR )
    fprintf( stderr, "error: pardiso analyze code %d\n", error_code ); // TODO decode
  assert( error_code == NO_ERROR );
}

void solver_factorize_pardiso( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  assert( A != NULL );
  solve_system_pardiso_t* const p = s->specific;
  assert( p != NULL );

  // TODO symmetric matrix handling
  if ( A->sym == SM_SYMMETRIC )
    convert_matrix_symmetry( A, BOTH );

  // prepare the matrix
  int ierr = convert_matrix( A, SM_CSR, FIRST_INDEX_ONE );
  assert( ierr == 0 );
  assert(( A->sym == SM_UNSYMMETRIC ) || ( A->sym == SM_SYMMETRIC ) );
  assert( A->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO
  assert( A->m == A->n ); // TODO can only handle square matrices at present (UMFPACK?)


  // save A for the solve step.. needed for iterative refinement // TODO add A to solve stage to allow iterative refinement!, remove this (and from other solvers)
  p->Arows = A->m;
  p->Acols = A->n;
  p->Ajj = ( int* ) A->jj;
  p->Aii = ( int* ) A->ii;
  p->Add = A->dd;


  int MAXFCT = 1;
  int MNUM = 1; // use which of MAXFACT stored factorizations in solution phase
  if ( p->MTYPE == REAL_UNSYM ) {
    assert(( A->sym == SM_UNSYMMETRIC ) || (( A->sym == SM_SYMMETRIC ) && ( A->location == BOTH ) ) );
    assert( A->data_type == REAL_DOUBLE );
  }
  else {
    assert( 0 ); // TODO
  }
  int PHASE = 22; // 10*i + j, do phases i->j, 1: analysis, 2; factorize, 3: solve w/ itr refinement, 4: terminate and release mem for MNUM, <0: release all mem
  // Note: needs data in A if IPARM(11) = 1 and IPARM(13)=1 or 2
  int N = A->m; // rows of A
  // TODO might need to save this for later?
  // expects A->dd to be "double" either real or complex
  // expects "int" for A->ii, A->jj
  int* const PERM = NULL; // if(IPARM(5)=1, user fill-reducing permutation (ordering) (size N) B = P A P^t, PERM(i) = row i of A to column PERM(i) of B
  int NRHS = 0; // columns in rhs
  int MSGLVL = 0; // TODO set based on s->verbosity
  double* const B = NULL; // N x NRHS, replaced w/ X if IPARM(6) = 1 -- solution phase
  double* const X = NULL; // N x NRHS solution if IPRAM(6) = 0 -- solution phase
  int error_code = 0;
  pardiso( p->PT, &MAXFCT, &MNUM, &( p->MTYPE ), &PHASE, &N,
           A->dd, ( int* ) A->ii, ( int* ) A->jj, PERM, &NRHS, p->iparm,
           &MSGLVL, B, X, &error_code, p->dparm );
  if ( error_code != NO_ERROR )
    fprintf( stderr, "error: pardiso factorize code %d\n", error_code ); // TODO decode
  assert( error_code == NO_ERROR );
}

// TODO b can be sparse... ??
// TODO can b be a matrix (vs a vector)?
void solver_evaluate_pardiso( solver_state_t* s, matrix_t* b, matrix_t* x ) {
  assert( s != NULL );
  assert( b != NULL );
  assert( b != x ); // TODO allow this form
  solve_system_pardiso_t* const p = s->specific;
  assert( p != NULL );

  // and we have a valid 'x' and 'b'
  int ierr = convert_matrix( b, DCOL, FIRST_INDEX_ONE ); // TODO support for sparse rhs too
  assert( ierr == 0 );
  assert( b->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO

  // TODO if solver only handles single rhs, loops solver and collect answers...

  int MAXFCT = 1;
  int MNUM = 1; // use which of MAXFACT stored factorizations in solution phase
  int PHASE = 33; // 10*i + j, do phases i->j, 1: analysis, 2; factorize, 3: solve w/ itr refinement, 4: terminate and release mem for MNUM, <0: release all mem
  // Note: needs data in A if IPARM(11) = 1 and IPARM(13)=1 or 2
  int* const N = &( p->Arows ); // rows of A
  // TODO might need to save this for later?
  // expects A->dd to be "double" either real or complex
  // expects "int" for A->ii, A->jj
  int* const PERM = NULL; // if(IPARM(5)=1, user fill-reducing permutation (ordering) (size N) B = P A P^t, PERM(i) = row i of A to column PERM(i) of B
  // TODO add checks in master function to handle x == b (put answer in b) that b is big enough to hold x (resize if required)
  // TODO add checks in master function that rows of b == rows of A, rows of x == columns of A
  // TODO add checks in master function that x (or b if b==x) is big enough to hold answer (columns of x == columns of b, rows of x == rows of A)
  int NRHS = b->n; // columns in rhs // TODO handle b == x
  int MSGLVL = 0; // TODO set based on s->verbosity
  double* const B = b->dd; // N x NRHS, replaced w/ X if IPARM(6) = 1 -- solution phase (Note: must be big enough for answer if reusing it!)

  // allocate data space // TODO if required
  // TODO move this to master function
  // TODO handle sparse answers??
  clear_matrix( x );
  x->format = DCOL;
  x->sym = SM_UNSYMMETRIC;
  x->data_type = b->data_type;
  x->m = p->Arows;
  x->n = b->n;
  x->nz = x->m * x->n;
  x->dd = malloc(( x->nz ) * sizeof( double ) ); // TODO or handle complex!

  assert( p->IPARM( 6 ) == 0 ); // TODO don't handle returning data in "b" yet

  double* const X = x->dd; // N x NRHS solution if IPRAM(6) = 0 -- solution phase
  int error_code = 0;
  pardiso( p->PT, &MAXFCT, &MNUM, &( p->MTYPE ), &PHASE, N,
           p->Add, p->Aii, p->Ajj, PERM, &NRHS, p->iparm,
           &MSGLVL, B, X, &error_code, p->dparm );
  if ( error_code != NO_ERROR )
    fprintf( stderr, "error: pardiso solve code %d\n", error_code ); // TODO decode
  assert( error_code == NO_ERROR );
}

void solver_finalize_pardiso( solver_state_t* s ) {
  if ( s == NULL )
    return;

  solve_system_pardiso_t* const p = s->specific;

  // release memory
  if ( p != NULL ) {
//    int MAXFCT = 1;
//    int MNUM = 1; // use which of MAXFACT stored factorizations in solution phase
//    int PHASE = -1; // 10*i + j, do phases i->j, 1: analysis, 2; factorize, 3: solve w/ itr refinement, 4: terminate and release mem for MNUM, <0: release all mem
//    int MSGLVL = 0; // TODO set based on s->verbosity
    int error_code = NO_ERROR;
// TODO this seems to segfault... but if we don't, things aren't cleaned up nicely
//    pardiso( p->PT, &MAXFCT, &MNUM, NULL, &PHASE, NULL,
//             NULL, NULL, NULL, NULL, NULL, NULL,
//             &MSGLVL, NULL, NULL, &error_code, NULL );
    if ( error_code != NO_ERROR )
      fprintf( stderr, "error: pardiso finalize code %d\n", error_code ); // TODO decode
    assert( error_code == NO_ERROR );

    free( p->PT );
    free( p->iparm );
    free( p->dparm );
  }
  free( p );
}
