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
#ifndef _SOLVER_LOOKUP_H_
#define _SOLVER_LOOKUP_H_

#include "solver_mumps.h"
#include "solver_umfpack.h"
#include "solver_cholmod.h"
#include "solver_taucs.h"
#include "solver_pardiso.h"
#include "solver_wsmp.h"


// the static solver lookup
struct  solver_properties_t {
  char* shortname; // lowercase, single word
  char* name; // full name
  char* author; // author
  char* organization; // author
  char* version; // version info
  // TODO find version from library rather than hard-coding it?
  char* license; // version info
  char* url;
  void ( *init )( solver_state_t* s );
  void ( *analyze )( solver_state_t*, matrix_t* );
  void ( *factorize )( solver_state_t*, matrix_t* );
  void ( *evaluate )( solver_state_t*, matrix_t*, matrix_t* );
  void ( *finalize )( solver_state_t* s );
  unsigned int capabilities; // flags (needs 32-bits)
  unsigned int multicore; // flags (needs 32-bits)
  char* references;
};


#define SOLVES_BASE_ZERO (1<<0)
#define SOLVES_BASE_ONE  (1<<1)

#define SOLVES_FORMAT_DROW (1<<2)
#define SOLVES_FORMAT_DCOL (1<<3)
#define SOLVES_FORMAT_COO  (1<<4)
#define SOLVES_FORMAT_CSC  (1<<5)
#define SOLVES_FORMAT_CSR  (1<<6)

#define SOLVES_SQUARE_ONLY                          (1<<7)
#define SOLVER_SYM_REQUIRES_DIAGONAL                (1<<8)
// Paradiso requires the entries on the diagonal to be allocated, even if they're zero! (WSMP too)
// This is for symmetric matrices that are not pos-def. (Pos-def matrices have entries on all diagonals.)

#define SOLVES_UNSYMMETRIC                          (1<<9)
#define SOLVES_SYMMETRIC_POSITIVE_DEFINITE_ONLY     (1<<10)
#define SOLVES_SYMMETRIC_UPPER_TRIANGULAR           (1<<11)
#define SOLVES_SYMMETRIC_LOWER_TRIANGULAR           (1<<12)
#define SOLVES_SYMMETRIC_BOTH                       (1<<13)
#define SOLVES_SYMMETRIC                            (SOLVES_SYMMETRIC_UPPER_TRIANGULAR | SOLVES_SYMMETRIC_LOWER_TRIANGULAR | SOLVES_SYMMETRIC_BOTH)
#define SOLVES_SKEW_SYMMETRIC_UPPER_TRIANGULAR      (1<<14)
#define SOLVES_SKEW_SYMMETRIC_LOWER_TRIANGULAR      (1<<15)
#define SOLVES_SKEW_SYMMETRIC                       (SOLVES_SKEW_SYMMETRIC_UPPER_TRIANGULAR | SOLVES_SKEW_SYMMETRIC_LOWER_TRIANGULAR)
#define SOLVES_HERMITIAN_SYMMETRIC_UPPER_TRIANGULAR (1<<16)
#define SOLVES_HERMITIAN_SYMMETRIC_LOWER_TRIANGULAR (1<<17)
#define SOLVES_HERMITIAN_SYMMETRIC                  (SOLVES_HERMITIAN_SYMMETRIC_UPPER_TRIANGULAR | SOLVES_HERMITIAN_SYMMETRIC_LOWER_TRIANGULAR)

#define SOLVES_DATA_TYPE_REAL_DOUBLE    (1<<18)
#define SOLVES_DATA_TYPE_REAL_SINGLE    (1<<19)
#define SOLVES_DATA_TYPE_COMPLEX_DOUBLE (1<<20)
#define SOLVES_DATA_TYPE_COMPLEX_SINGLE (1<<21)

#define SOLVES_RHS_BASE_ZERO (1<<22)
#define SOLVES_RHS_BASE_ONE  (1<<23)

#define SOLVES_RHS_DROW (1<<24)
#define SOLVES_RHS_DCOL (1<<25)
#define SOLVES_RHS_COO  (1<<26)
#define SOLVES_RHS_CSC  (1<<27)
#define SOLVES_RHS_CSR  (1<<28)

#define SOLVES_RHS_VECTOR_ONLY (1<<29)


// how multicore is each solver?
#define SOLVER_SINGLE_THREADED_ONLY 0

#define SOLVER_REQUIRES_MPI (1<<0)
#define SOLVER_CAN_USE_MPI  (1<<1)

#define SOLVER_REQUIRES_OMP (1<<2)
#define SOLVER_CAN_USE_OMP  (1<<3)


// not the actual longest string but something to make it safe to compare strings from the command line
#define SOLVER_SHORTNAME_MAX_LEN 15

static const struct solver_properties_t solver_lookup[] = {
  { "umfpack", "UMFPACK", "Tim Davis et al", "University of Florida", "5.5.0", "GPL",
    "http://www.cise.ufl.edu/research/sparse/umfpack",
    // TODO 5.5.1 is available
    &solver_init_umfpack,
    &solver_analyze_umfpack,
    &solver_factorize_umfpack,
    &solver_evaluate_umfpack,
    &solver_finalize_umfpack,
    SOLVES_FORMAT_CSC | SOLVES_BASE_ZERO | SOLVES_UNSYMMETRIC |
    SOLVES_DATA_TYPE_REAL_DOUBLE |
    SOLVES_SQUARE_ONLY | // TODO is umfpack really restricted to square matrices?
    SOLVES_RHS_DCOL | SOLVES_RHS_VECTOR_ONLY,
    SOLVER_SINGLE_THREADED_ONLY,
    // Note: Adjacent constant strings will be concatentated
    "    * A column pre-ordering strategy for the unsymmetric-pattern multifrontal method,\n"
    "      T. A. Davis, ACM Transactions on Mathematical Software,\n"
    "      vol 30, no. 2, June 2004, pp. 165-195.\n"
    "    * Algorithm 832: UMFPACK, an unsymmetric-pattern multifrontal method,\n"
    "      T. A. Davis, ACM Transactions on Mathematical Software,\n"
    "      vol 30, no. 2, June 2004, pp. 196-199.\n"
    "    * A combined unifrontal/multifrontal method for unsymmetric sparse matrices,\n"
    "      T. A. Davis and I. S. Duff, ACM Transactions on Mathematical Software,\n"
    "      vol. 25, no. 1, pp. 1-19, March 1999.\n"
    "    * An unsymmetric-pattern multifrontal method for sparse LU factorization,\n"
    "      T. A. Davis and I. S. Duff, SIAM Journal on Matrix Analysis and Applications,\n"
    "      vol 18, no. 1, pp. 140-158, Jan. 1997.\n" },

  { "mumps", "MUMPS", "Patrick Amestoy et al", "Université de Toulouse, et. al", "4.9.2", "public domain",
    "http://graal.ens-lyon.fr/MUMPS",
    &solver_init_mumps,
    &solver_analyze_mumps,
    &solver_factorize_mumps,
    &solver_evaluate_mumps,
    &solver_finalize_mumps,
    SOLVES_FORMAT_COO | SOLVES_BASE_ONE | SOLVES_UNSYMMETRIC |
    SOLVES_DATA_TYPE_REAL_DOUBLE |
    SOLVES_SQUARE_ONLY |
    // TODO is mumps really restricted to square matrices?
    // TODO mumps can handle sparse rhs and can solve multiple right hand sides
    SOLVES_RHS_DCOL | SOLVES_RHS_VECTOR_ONLY,
    SOLVER_REQUIRES_MPI,
    "    [1] P. R. Amestoy, I. S. Duff, J. Koster and J.-Y. L'Excellent,\n"
    "        A fully asynchronous multifrontal solver using distributed dynamic scheduling,\n"
    "        SIAM Journal of Matrix Analysis and Applications, Vol 23, No 1, pp 15-41 (2001).\n"
    "    [2] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and S. Pralet,\n"
    "        Hybrid scheduling for the parallel solution of linear systems.\n"
    "        Parallel Computing Vol 32 (2), pp 136-156 (2006).\n" },

  { "cholmod", "CHOLMOD", "Tim Davis, William Hager", "University of Florida", "1.7.1", "LGPL",
    // version 1.7.3 is available...
    "http://www.cise.ufl.edu/research/sparse/cholmod",
    &solver_init_cholmod,
    &solver_analyze_cholmod,
    &solver_factorize_cholmod,
    &solver_evaluate_cholmod,
    &solver_finalize_cholmod,
    SOLVES_FORMAT_CSC | SOLVES_BASE_ZERO |
    // TODO can handle COO matrices and dense matrices (DROW?)
    SOLVES_SQUARE_ONLY |
    SOLVES_SYMMETRIC_POSITIVE_DEFINITE_ONLY |
    SOLVES_SYMMETRIC |
    // upper or lower triangular symmetric or BOTH or unsymmetric but must still be SPD
    // TODO keep track of when a matrices' entries have been sorted!
    SOLVES_DATA_TYPE_REAL_DOUBLE |
    // TODO is cholmod really restricted to square matrices?
    SOLVES_RHS_DCOL | SOLVES_RHS_CSC,
    SOLVER_SINGLE_THREADED_ONLY,
    "    * Dynamic supernodes in sparse Cholesky update/downdate and triangular\n"
    "      solves, T. A. Davis and W. W. Hager, ACM Trans. Math. Software,\n"
    "      Vol 35, No. 4, 2009. (as CISE Tech Report)\n"
    "    * Algorithm 887: CHOLMOD, supernodal sparse Cholesky factorization and\n"
    "      update/downdate, Y. Chen, T. A. Davis, W. W. Hager, and \n"
    "      S. Rajamanickam, ACM Trans. Math. Software, Vol 35, No. 3, 2009.\n"
    "      (as CISE Tech Report)\n"
    "    * Row modifications of a sparse Cholesky factorization, T. A. Davis\n"
    "      and W. W. Hager, SIAM Journal on Matrix Analysis and Applications,\n"
    "      vol 26, no 3, pp. 621-639, 2005.\n"
    "    * Multiple-rank modifications of a sparse Cholesky factorization,\n"
    "      T. A. Davis and W. W. Hager, SIAM Journal on Matrix Analysis and\n"
    "      Applications, vol. 22, no. 4, pp. 997-1013, 2001.\n"
    "    * Modifying a sparse Cholesky factorization, T. A. Davis and W. W.\n"
    "      Hager, SIAM Journal on Matrix Analysis and Applications, vol. 20,\n"
    "      no. 3, pp. 606-627, 1999.\n" },


  // "Matlab 7 will use TAUCS' in-core sparse Cholesky factorization within the backslash linear solver." (TODO true?)
  // effectively single threaded
  // mostly symmetric solvers: out-of-core and in-core, multithreaded w/ CILK (SMP -- taucs.cilk.nproc=N) + 1 out-of-core unsym solver (slow?)...
  // entry added 2011-03-01, software last updated Sept, 2003 (v2.2)
  { "taucs", "TAUCS", "Sivan Toledo", "Tel-Aviv University", "2.2", "LGPL",
    "http://www.tau.ac.il/~stoledo/taucs/",
    &solver_init_taucs,
    &solver_analyze_taucs,
    &solver_factorize_taucs,
    &solver_evaluate_taucs,
    &solver_finalize_taucs,
    SOLVES_FORMAT_CSC | SOLVES_BASE_ZERO |
    SOLVES_SYMMETRIC_UPPER_TRIANGULAR | // SOLVES_UNSYMMETRIC | // TODO or lower triangular?
    SOLVES_DATA_TYPE_REAL_DOUBLE | // TODO and REAL_COMPLEX, double and single precision, and handles hermitian
    SOLVES_RHS_DCOL,
    0, // uses CILK? not MPI or openMP
    "        TODO citations\n" }, // TODO


  // can calculate bit-identical solution on multicore vs. clusters of multicores
  // MPI based solver is for symmetric indefinite!
  // supports: unsymmetric, structurally symmetric, real/complex, hermitian
  // LU w/ complete pivoting
  // combines iterative and direct solvers... for very large 3D problems
  // integrated METIS -- what version?
  // email: pardiso-informatik@unibas.ch about where this solver is used"
  // license is specific to host & username
  // uses left and right-looking L3 BLAS w/ supernodes (openMP and MPI)
  //  symmetric: orderings (min-degree or nested-disection METIS), then
  //     parallel left-right Choleksy (PAP^t=LL^t or indefinite PAP^t = LDL^t)
  //     1x1 and 2x2 Bunch-Kaufman pivoting for sym-indef, then fwd-backward subst & itr refinement
  //     perturbs coefficient matrix when pivots can't be found, corrected by iterative refinement steps
  //      -- as accurate as complete sparse pivoting techniques
  //     can improve pivoting accuracy by identifying large entries and permuting closer to diagonal (fewer pivot perturbations required)
  //       -- maximum weighted matchings
  //     also computes inertia for real sym indef matrices
  //   structurally sym: PAP^t = QLU^t do symmetric fill-in reducing ordering, then apply to parallel unsym factorization
  //     partial pivoting at supernodes and itr refinement
  //   unsym: permute and scale (large entries near diagonal to improve numerical reliability), fill-in reducing permutation (P_mps A +(P_mps A)^t)
  //     parallel factorization QLUR = A' = P P_mps D_r A D_c P (super-node pivoting Q R) P=P P_mps to keep "sufficiently large cycles of P_mps in one diagonal block
  //     can't factorize anymore? pivot perturbation
  //     pivot thres \alpha = e ||A'||_inf (e=machine epsilon), tiny pivots set to sgn(l)*e*||A'|_inf (numerical stability vs. small pivots)
  //      -- so needs iterative refinement to get a good answer since the factorization has errors
  //   for real sym indef matrices: MPI based solver
  //   for sym indef: can use *iterative* solver (weighted matchings, algebraic multilevel incomplete factorization) (Krylov subspace techniques) (multilevel incomplete factorization preconditioners)
  //   unsym: combine direct and iterative for unsym (same sparsity pattern, slowly changing system)
  //     solves first factorization to LU, then uses these as preconditioned krylov subspace iterations, switch back if not converging
  //     IPARM(4), IPARM(20)
  { "pardiso", "Pardiso", "Olaf Schenk, Klaus Gärtner", "University Basel", "4.1.0", "academic/commercial",
    "http://www.pardiso-project.org",
    &solver_init_pardiso,
    &solver_analyze_pardiso,
    &solver_factorize_pardiso,
    &solver_evaluate_pardiso,
    &solver_finalize_pardiso,
    SOLVES_FORMAT_CSR | SOLVES_BASE_ONE | SOLVER_SYM_REQUIRES_DIAGONAL | // TODO diagonal only for sym matrices?
    SOLVES_SYMMETRIC_UPPER_TRIANGULAR | SOLVES_UNSYMMETRIC | // TODO or lower triangular?
    SOLVES_DATA_TYPE_REAL_DOUBLE | // TODO and REAL_COMPLEX
    SOLVES_RHS_DCOL | SOLVES_RHS_CSR,
    SOLVER_REQUIRES_OMP | SOLVER_CAN_USE_MPI,
    "    [1] O. Schenk and K. Gärtner, Solving Unsymmetric Sparse Systems of Linear\n"
    "        Equations with PARDISO, Journal of Future Generation Computer Systems,\n"
    "        20(3):475--487, 2004.\n"
    "    [2] O. Schenk and K. Gärtner, On fast factorization pivoting methods for\n"
    "        symmetric indefinite systems, Elec. Trans. Numer. Anal., 23:158--179, 2006.\n"
    "    [3] G.Karypis and V.Kumar, A fast and high quality multilevel scheme for\n"
    "        partitioning irregular graphs, SIAM Journal on Scientific Computing, 1998\n"
    "        (20) 1, 359-392\n"
    "                            Version 4.1.0 related:\n"
    "    [4] O. Schenk, M. Bollhoefer, and R. Roemer, On large-scale diagonalization\n"
    "        techniques for the Anderson model of localization. SIAM Review 50 (2008),\n"
    "        pp. 91-112.\n"
    "    [5] O. Schenk, A. Waechter, and M. Hagemann, Matching-based Preprocessing\n"
    "        Algorithms to the Solution of Saddle-Point Problems in Large-Scale\n"
    "        Nonconvex Interior-Point Optimization. Journal of Computational\n"
    "        Optimization and Applications, pp. 321-341, Volume 36, Numbers 2-3 /\n"
    "        April, 2007.\n" },

  { "wsmp", "WSMP", "Anshul Gupta", "IBM/Univ. Minnesota", "11.01.19", "commercial?",
    "http://www-users.cs.umn.edu/~agupta/wsmp.html",
    &solver_init_wsmp,
    &solver_analyze_wsmp,
    &solver_factorize_wsmp,
    &solver_evaluate_wsmp,
    &solver_finalize_wsmp,
    SOLVES_FORMAT_CSR | SOLVES_BASE_ZERO | SOLVER_SYM_REQUIRES_DIAGONAL | // TODO CSC lower triangular, base 1
    SOLVES_SYMMETRIC_UPPER_TRIANGULAR | SOLVES_UNSYMMETRIC | // TODO solvers
    SOLVES_DATA_TYPE_REAL_DOUBLE | // TODO and REAL_COMPLEX -- 8Byte floating pt, 4 byte ints
    SOLVES_RHS_DCOL,
    SOLVER_CAN_USE_OMP | SOLVER_REQUIRES_MPI, // TODO: SOLVER_CAN_USE_MPI (select non-MPI solver...) // TODO non-MPI/OMP solvers (switch, based on number of threads/nodes)
    "    [1] A. Gupta, G. Karypis, V. Kumar, A Highly Scalable Parallel Algorithm for\n"
    "        Sparse Matrix Factorization, IEEE Transactions on Parallel and\n"
    "        Distributed Systems, 8(5):502–520, May 1997.\n" },

  // Note: MUST have null entry at the end of the list!
  {0}
};

#endif
