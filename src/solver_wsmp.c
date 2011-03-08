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
#include "solver_wsmp.h"
#include "solvers.h"
#include "matrix.h"
//#include <wsmp.h>

#include <stdlib.h> // malloc, free
#include <string.h> // memcpy
#include <assert.h>
#include <stdint.h> // int64_t
#include <unistd.h> // dup -> redirect stdout temporarily

#include <mpi.h>

// NOTE: wsmp doesn't include a header file for its library, use these prototypes instead
void wsetmaxthrds_(const int *);
void wsetmpicomm_(MPI_Fint *);
void wsmp_initialize();
void pwsmp_initialize();

void wsmp_clear();
void pwsmp_clear();

// TODO wgsmp, wssmp, pwssmp
void pwgsmp_( int* N, int* IA, int* JA, double* AVALS, double* B, int* LDB, int* NRHS, double* RMISC, int* IPARM, double* DPARM );
// N: col/rows per process (on this process)
// IA: N+1 row/col ptrs into JA, if N=0, JA[0]=0
// JA: NZ col/row indices
// RMISC only used if IPARM(25)=1

enum wsmp_solver_stages { WSMP_INIT = 0, WSMP_ANALYZE = 1, WSMP_FACTORIZE = 2, WSMP_EVALUATE = 3, WSMP_ITERATIVE_REFINEMENT = 4 }; // IPARM(2,3)
enum wsmp_solver_format { WSMP_CSR = 0, WSMP_CSC = 1 }; // IPARM(4)
enum wsmp_solver_errors { NO_ERROR = 0 };


#define IPARM(I) iparm[(I)-1]
#define DPARM(I) dparm[(I)-1]
typedef struct {
  int Arows;
  int Acols;
  int* Aii; // these are pointers to existing data (DON'T free)
  int* Ajj;
  double* Add;

  // Paradiso config
  int *iparm; // int * 64 - config
  double *dparm; // double * 64 - config
} solve_system_wsmp_t;

void solver_init_wsmp( solver_state_t* s ) {
  assert( s != NULL );
  solve_system_wsmp_t * const p = calloc( 1, sizeof( solve_system_wsmp_t ) );
  assert( p != NULL );
  s->specific = p;

  p->iparm = calloc( 64, sizeof( int ) );
  p->dparm = calloc( 64, sizeof( double ) );

  // TODO Note: as of version 10.1, IPARM(23, 24) have changed along with wsetmaxthrds and thread control, env variable WSMP_NUM_THREADS is new in version 9.8
  // TODO Note: requires a license file wsmp.lic in current directory or at WSMPLICPATH env variable
  // TODO can set work space to TMPDIR (env variable) defaults to /tmp
  // TODO Note: BLAS *must* run in serial threaded mode (bad performance) -- set OMP_NUM_THREADS = 1 MKL_NUM_THREADS, GOTO_NUM_THREADS
  // TODO might need to increase default stack/data size (hang/segfault) -- IPARM(64) = -102 -- can use command "limit" to increase
  // TODO can download memchck.c to see how much stack and data space is available
  // TODO note mpich or derivatibes may not be thread safe? (fine for symmetric solver, must use wsetmaxthrds = 1 to use only one thread per process) symm solver needs MPI_THREADS_FUNNELED and unsym solver needs MPI_THREAD_MULTIPLE support -- can force to MPI_THREAD_SINGLE
  // TODO set env var MALLOC_TRIM_THRESHOLD_ = -1 and MALLOC_MMAP_MAX_ = 0 except when using PWSMP: then MALLOC_TRIM_THRESHOLD=-1 can cause problems (Crashes) when there are more mpi processes on the same node -- set to MALLOC_TRIM_THRESHOLD=134217728 or wsetmaxthrds=1 and no other wsmp processes running on the same node
  // TODO suggests MKL, GOTO, ATLAS
  // TODO should use specific MPI for myrinet, infiniband, etc (vendor supplied)

  p->IPARM( 1 ) = 1; // 0: fill in default values (4..64), 2,3 MUST be set .. this will overwrite the following config with defaults (it mostly here for documentation)
  p->IPARM( 2 ) = 2; // start task, on output has next task (end task +1)
  p->IPARM( 3 ) = 1; // end task
  p->IPARM( 4 ) = 0; // matrix format (0:CSR, 1: CSC)
  p->IPARM( 5 ) = 1; // numbering style (base 0 or 1)
  p->IPARM( 6 ) = 3; // max iterative refinement steps
  p->IPARM( 7 ) = 3; // residual norm type (0, 1, 2, 3: double precision (same as remainder), 4, 5, 6, 7: residual in quadruple precision (double remainder) (0,4: do exactly IPARM(6) itr. refinement steps, 0: don't calculate DPARM(7) - backward error) (1,2,3,5,6,7: itr refinement until IPARM(6) - 2-norm, 1,5: 1-norm, 2,6: 2-norm, 3,7: inf.-norm), if NRHS >1 then max backward error is returned, if IPARM(10)=1 then err is wrt scaled system
  p->IPARM( 8 ) = 0; // max. matching use 0:auto, 1: permute, 2: don't permute - max. weight matching row permutations to maximize diagonal
  p->IPARM( 9 ) = 1; // scaling w/o matching  - simple equilibriation to set diagonal to 1.0
  p->IPARM( 10 ) = 1; // scaling w/ matching - 1: do scaling, if IPARM(9) is set too, IPARM(10) takes priority
  p->IPARM( 11 ) = 1; // thresh. pivoting opt. (1 && DPARM(11)=0.0 -- auto threshold)
  p->IPARM( 12 ) = 0; // pivot perturb. opt.
  // p->IPARM(13) = ? // output: # row/col exchanges
  // p->IPARM(14) = ? // output: # perturbations
  p->IPARM( 15 ) = 25; // # factorizations -- chooses how much effort to put into analysis, 0:lots, 1: only once, small small number: a few times
  p->IPARM( 16 ) = 1; // ordering option 1 (-1: no ordering, -2: reverse Cuthill-KcKee, >0: graph partitioning, 0: use all options, speed=3) (1,2,3: use IPARM(17--20) (1: slowest/best, 2: intermediate, 3: fastest/worst) (only a couple of solns? use 2 or 3...)
  p->IPARM( 17 ) = 0; // ordering option 2 (max subgraph nodes for using min local fill w/o further partitioning) (0: auto, >50--200 might be better in some cases?)
  p->IPARM( 18 ) = 0; // ordering option 3 (1: force ordering to compute min local fill as well as graph bisection, choose best) -- ignored for MPI (bad performance)
  p->IPARM( 19 ) = 0; // ordering option 4 random number seed for choosing ordering
  p->IPARM( 20 ) = 0; // ordering option 5 known matrix characteristics 0:auto, 1: a few dense rows/cols, 2: FEM like graph (many degrees of freedom, repeated structure, uses a compressed graph -- much faster!)
  p->IPARM( 21 ) = 1; // block triangular form (0: suppress block triangulation, sometimes helps with soln accuracy)
  // p->IPARM(22) = ? // # blocks in block triangular form
  // p->IPARM(23) = ? // factorization nnz_l + nnz_u (in thousands) -- WSMP uses relaxed supernodes to maximize BLAS-lvl-3 operations -- artificial zeros into matrices to improve BLAS efficiency
  // p->IPARM(24) = ? // symbolic nnz_l + nnz_u (in thousands)
  p->IPARM( 25 ) = 0; // RMISC use (1: return component-wise backward error (iterative refinement stage) in RMISC where RMISC is "double *RMISC[N];")
  // p->IPARM(26) = ? // #itr ref steps
  p->IPARM( 27 ) = 0; // # fact. before re-analyze (0: auto, +N: re-analyze at-least every N factorizations)
  p->IPARM( 28 ) = 0; // rook pivoting (note S, T) == ??? (0: default pivoting, choose pivot so |diagonal elem| >= pivot thres * max pivot column, 1: check row and column (rook pivoting, much slower, might need IPARM(21)=0) -- increased accuracy) -- not available in MPI version
  p->IPARM( 29 ) = 0; // garbage collection (1: return with properly size factorization, garbage collect) (might improve solve speed for multiple solves)
  p->IPARM( 30 ) = 0; // solve option (x overwrites b) (0: x = A^-1 b, 1: x = L^-1 b, 2: x = U^-1 b, 4: x = A^-T b, 5: x = U^-T b, 6: x = L^-T b) calling solve(1) + solve(2) = solve(0) w/o itr refinement, no backward error, requires IPARM(21)=0 (slower) -- also (4,5,6) == IPARM(4)=0 + (0,1,2)
  p->IPARM( 31 ) = 1; // # solves per factor (!NRHS, its calls to the evaluate stage)
  p->IPARM( 32 ) = 0; // block size (note P) == ???, MPI only, block size for dense matrix computations for 2d decomposition of frontal and update matrices (0: auto, N: use <= 2^N
  // p->IPARM(33) = ? // output: # CPUs used in SMP mode (threads), local to each MPI process
  p->IPARM( 34 ) = 10; // DAG manip. option (note T,P) == ??? , to improve load-imbalance vs. fill-in (0: none, 0--log2(NCPUs): communication/load-imbalance reduction effort)
  // p->IPARM(35 -- 63) // reserved, must = 0
  // p->IPARM(64) = ? // output: error code (0: success, lowest 3 digits: err code, most-siginificant digits: MPI process) all processes get the same code, <0: invalid input, |n| is IPARM(n) err location, < -99 non-numerical runtime error, -102 = mem allocation err (see note 3.2,3.4,3.5), -200 = MPI -- problem too small for number of MPI processes, try DPARM(25), also if MPI isn't initialized, -300: previous stages not completed, -700: internal error (please report to wsmp@watson.ibm.com) can be due to indicies > 2^31, some platforms have 8-byte integer version available, -900: license expired, invalid, missing, +1--N: computational err (first pivot to = 0)

  // other fields unused/reserved  = 0
  // p->DPARM(1--3) unused
  // p->DPARM(4) = ? // output: largest pivot
  // p->DPARM(5) = ? // output: smallest pivot
  p->DPARM( 6 ) = 2e-15; // backward substitution error limit
  // p->DPARM(7) = ? // output: backward error
  p->DPARM( 10 ) = 1e-18; // backward substitution error limit
  p->DPARM( 11 ) = 0.01; // pivot threshold
  p->DPARM( 12 ) = 2e-8; // small pivot threshold
  // p->DPARM(13) = ? // output: number of super nodes
  // p->DPARM(14) = ? // output: # data-DAG edges
  // p->DPARM(15--20) unused
  // p->DPARM(21) = ? // output: structural symmetry
  p->DPARM( 22 ) = 2e-8; // small pivot replacement
  // p->DPARM(23) = ? // output: actual factorization operations
  // p->DPARM(24) = ? // output: symbolic factorization operations
  p->DPARM( 25 ) = 5e6; // 1.0; // TODO FIXME default is 5e6; // minimum parallel task size // TODO FIXME -- forcing use of parallel solve even when its not optimal! (avoid error_code = -200)
  p->DPARM( 26 ) = 1.0; // supernode amalgamation
  p->DPARM( 27 ) = 1.0; // re-analyze condition
  // p->DPARM(28 -- 32) unused
  // p->DPARM(33) = ? // output: load imbalance
  // p->DPARM(34) unused
  // p->DPARM(35 -- 63) reserved = 0
  // p->DPARM(64) unused


  // must call wsetmaxthrds_() before initialize
//  MPI_Fint inpcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
//  wsetmpicomm_(&inpcomm);
  const int threads = 1;
  wsetmaxthrds_(&threads); // TODO FIXME ... if using MPI decide if MPI is thread safe, if not, must limit wsmp to one thread in unsymmetric solver (safe for symmetric solver)
  // TODO p?wsmp_intialize()

  // make wsmp quiet
  FILE *o = stderr;
  if(s->verbosity < 4) {
    stderr = fopen( "/dev/null", "w" );
    assert( stderr != NULL );
  }

  pwsmp_initialize(); // MPI

  // undo quieting
  if(s->verbosity < 4) {
    int ret = fclose( stderr );
    assert( ret == 0 );
    stderr = o;
  }
}

// TODO split analyze stage into ordering and symbolic factorization stages?
void solver_analyze_wsmp( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  assert( A != NULL );
  solve_system_wsmp_t* const p = s->specific;
  assert( p != NULL );

  // TODO could do this smarter by giving MUMPS the symmetric matrix to solve,
  // instead we're just going straight for the unsymmetric solver
  // TODO this logic should really be moved to the wrapper solver() functions using the capabilities masks
  if ( A->sym == SM_SYMMETRIC )
    convert_matrix_symmetry( A, BOTH );
  // TODO for SYMMETRIC matrices: wsmp wants all diagonal entries to exist in A, even if they are zero

  // prepare the matrix
  int ierr = convert_matrix( A, SM_CSR, FIRST_INDEX_ONE );
  assert( ierr == 0 );
  assert(( A->sym == SM_UNSYMMETRIC ) || ( A->sym == SM_SYMMETRIC ) );
  assert( A->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO
  assert( A->m == A->n ); // TODO can only handle square matrices at present???

  // Compressed Row Format
  assert( A->ii[0] == 1 );
  assert( A->ii[A->n] == A->nz + 1 );

  int N = A->m; // rows of A
  p->IPARM(2) = WSMP_ANALYZE;
  p->IPARM(3) = WSMP_ANALYZE;
  int TEST = 1;
  pwgsmp_( &N, ( int* ) A->ii, ( int* ) A->jj, A->dd, NULL, &TEST, &TEST, NULL, p->iparm, p->dparm);
  const int error_code = p->IPARM(64);
  if ( error_code != NO_ERROR )
    fprintf( stderr, "error: wsmp analyze code %d\n", error_code ); // TODO decode
  assert( error_code == NO_ERROR );
}

void solver_factorize_wsmp( solver_state_t* s, matrix_t* A ) {
  assert( s != NULL );
  assert( A != NULL );
  solve_system_wsmp_t* const p = s->specific;
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

  int N = A->m; // rows of A
  p->IPARM(2) = WSMP_FACTORIZE;
  p->IPARM(3) = WSMP_FACTORIZE;
  int TEST = 1;
  pwgsmp_( &N, ( int* ) A->ii, ( int* ) A->jj, A->dd, NULL, &TEST, &TEST, NULL, p->iparm, p->dparm);
  const int error_code = p->IPARM(64);
  if ( error_code != NO_ERROR )
    fprintf( stderr, "error: wsmp factorize code %d\n", error_code ); // TODO decode
  assert( error_code == NO_ERROR );
}

// TODO b can be sparse... ??
// TODO can b be a matrix (vs a vector)?
void solver_evaluate_wsmp( solver_state_t* s, matrix_t* b, matrix_t* x ) {
  assert( s != NULL );
  assert( b != NULL );
  assert( b != x ); // TODO allow this form
  solve_system_wsmp_t* const p = s->specific;
  assert( p != NULL );

  // and we have a valid 'x' and 'b'
  int ierr = convert_matrix( b, DCOL, FIRST_INDEX_ONE ); // TODO support for sparse rhs too
  assert( ierr == 0 );
  assert( b->data_type == REAL_DOUBLE ); // don't handle complex... yet TODO

  // TODO if solver only handles single rhs, loops solver and collect answers...

  // allocate data space // TODO if required
  // TODO move this to master function
  if(b == x) {
    // go ahead, we're expecting b to be destroyed
    // need to put data from b into x, clear x's ptr
    clear_matrix( x );
    x->format = DCOL;
    x->sym = SM_UNSYMMETRIC;
    x->data_type = b->data_type;
    x->m = p->Arows;
    x->n = b->n;
    x->nz = x->m * x->n;
    if(b->nz > x->nz) {
      x->dd = malloc(b->nz * sizeof( double ));
    }
    else {
      x->dd = malloc(x->nz * sizeof( double ));
    }
    x->dd = b->dd;
    b->dd = NULL;
    clear_matrix(b);
  }
  else {
    // need to copy b since it's destroyed in the process
    clear_matrix(x);
    matrix_t* t = copy_matrix(b);
    // push copied t contents into x
    x->format = t->format;
    x->sym = t->sym;
    x->data_type = b->data_type;
    x->m = p->Arows;
    x->n = b->n;
    x->nz = x->m * x->n;
    x->dd = t->dd;
    free(t);
    if(b->nz > x->nz) {
      double *t = realloc(x->dd, b->nz * sizeof( double ));
      assert(t != NULL); // realloc failure -- TODO proper error code
      x->dd = t;
    }
    assert(x->dd != NULL);
  }

  int N = p->Arows; // rows of A
  double* B = x->dd; // N x NRHS
  int LDB = b->m;  // rows of B, must be >= N
  int NRHS = b->n; // columns of B
  p->IPARM(2) = WSMP_EVALUATE;
//  p->IPARM(3) = WSMP_ITERATIVE_REFINEMENT; // TODO split into seperate stage
  p->IPARM(3) = WSMP_EVALUATE; // TODO split into seperate stage
  double TEST = 1;
  assert(&N != NULL);
  assert(p->Aii != NULL);
  assert(p->Ajj != NULL);
  assert(p->Add != NULL);
  assert(B != NULL);
  assert(&LDB != NULL);
  assert(&NRHS != NULL);
  assert(&TEST != NULL);
  assert(p->iparm != NULL);
  assert(p->dparm != NULL);
  pwgsmp_( &N, ( int* ) p->Aii, ( int* ) p->Ajj, p->Add, B, &LDB, &NRHS, &TEST, p->iparm, p->dparm);
  const int error_code = p->IPARM(64);
  if ( error_code != NO_ERROR )
    fprintf( stderr, "error: wsmp evaluate code %d\n", error_code ); // TODO decode
  assert( error_code == NO_ERROR );
}

void solver_finalize_wsmp( solver_state_t* s ) {
  if ( s == NULL )
    return;

  solve_system_wsmp_t* const p = s->specific;

  // release memory
  if ( p != NULL ) {
    int error_code = NO_ERROR;

    // for symmetric matrices its wssmp/pwssmp
    // wsmp_clear(); // for SMP or
    pwsmp_clear(); // MPI
    if ( error_code != NO_ERROR )
      fprintf( stderr, "error: wsmp finalize code %d\n", error_code ); // TODO decode
    assert( error_code == NO_ERROR );

    free( p->iparm );
    free( p->dparm );
  }
  free( p );
  s->specific = NULL;
}
