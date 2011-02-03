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

#define SOLVES_UNSYMMETRIC                          (1<<8)
#define SOLVES_SYMMETRIC_POSITIVE_DEFINITE_ONLY     (1<<9)
#define SOLVES_SYMMETRIC_UPPER_TRIANGULAR           (1<<10)
#define SOLVES_SYMMETRIC_LOWER_TRIANGULAR           (1<<11)
#define SOLVES_SYMMETRIC_BOTH                       (1<<12)
#define SOLVES_SYMMETRIC                            (SOLVES_SYMMETRIC_UPPER_TRIANGULAR | SOLVES_SYMMETRIC_LOWER_TRIANGULAR | SOLVES_SYMMETRIC_BOTH)
#define SOLVES_SKEW_SYMMETRIC_UPPER_TRIANGULAR      (1<<13)
#define SOLVES_SKEW_SYMMETRIC_LOWER_TRIANGULAR      (1<<14)
#define SOLVES_SKEW_SYMMETRIC                       (SOLVES_SKEW_SYMMETRIC_UPPER_TRIANGULAR | SOLVES_SKEW_SYMMETRIC_LOWER_TRIANGULAR)
#define SOLVES_HERMITIAN_SYMMETRIC_UPPER_TRIANGULAR (1<<15)
#define SOLVES_HERMITIAN_SYMMETRIC_LOWER_TRIANGULAR (1<<16)
#define SOLVES_HERMITIAN_SYMMETRIC                  (SOLVES_HERMITIAN_SYMMETRIC_UPPER_TRIANGULAR | SOLVES_HERMITIAN_SYMMETRIC_LOWER_TRIANGULAR)

#define SOLVES_DATA_TYPE_REAL_DOUBLE    (1<<17)
#define SOLVES_DATA_TYPE_REAL_SINGLE    (1<<18)
#define SOLVES_DATA_TYPE_COMPLEX_DOUBLE (1<<19)
#define SOLVES_DATA_TYPE_COMPLEX_SINGLE (1<<20)

#define SOLVES_RHS_BASE_ZERO (1<<21)
#define SOLVES_RHS_BASE_ONE  (1<<22)

#define SOLVES_RHS_DROW (1<<23)
#define SOLVES_RHS_DCOL (1<<24)
#define SOLVES_RHS_COO  (1<<25)
#define SOLVES_RHS_CSC  (1<<26)
#define SOLVES_RHS_CSR  (1<<27)

#define SOLVES_RHS_VECTOR_ONLY (1<<28)


// how multicore is each solver?
#define SOLVER_SINGLE_THREADED_ONLY 0

#define SOLVER_REQUIRES_MPI (1<<0)
#define SOLVER_CAN_USE_MPI  (1<<1)

#define SOLVER_REQUIRES_OMP (1<<2)
#define SOLVER_CAN_USE_OMP  (1<<3)


// Note: remember to upate this value when adding more solvers!
#define SOLVER_INDEX_COUNT 3
#define SOLVER_SHORTNAME_MAX_LEN 7

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
    "      no. 3, pp. 606-627, 1999.\n" }

};

// Note: did the SOLVER_INDEX_COUNT get updated when adding more solvers??

#endif
