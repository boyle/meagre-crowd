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
  void (*analyze)(solver_state_t*, matrix_t*);
  void (*factorize)(solver_state_t*, matrix_t*);
  void (*evaluate)(solver_state_t*, matrix_t*, matrix_t*);
  unsigned int capabilities; // flags (needs 32-bits)
  char* references;
};


#define SOLVES_BASE_ZERO (1<<0)
#define SOLVES_BASE_ONE  (1<<1)

#define SOLVES_FORMAT_DROW (1<<2)
#define SOLVES_FORMAT_DCOL (1<<3)
#define SOLVES_FORMAT_COO  (1<<4)
#define SOLVES_FORMAT_CSC  (1<<5)
#define SOLVES_FORMAT_CSR  (1<<6)

#define SOLVES_SQUARE                               (1<<7)

#define SOLVES_UNSYMMETRIC                          (1<<8)
#define SOLVES_SYMMETRIC_POSITIVE_DEFINITE_ONLY     (1<<9)
#define SOLVES_SYMMETRIC_UPPER_TRIANGULAR           (1<<10)
#define SOLVES_SYMMETRIC_LOWER_TRIANGULAR           (1<<11)
#define SOLVES_SYMMETRIC                            (SOLVES_SYMMETRIC_UPPER_TRIANGULAR | SOLVES_SYMMETRIC_LOWER_TRIANGULAR) 
#define SOLVES_SKEW_SYMMETRIC_UPPER_TRIANGULAR      (1<<12)
#define SOLVES_SKEW_SYMMETRIC_LOWER_TRIANGULAR      (1<<13)
#define SOLVES_SKEW_SYMMETRIC                       (SOLVES_SKEW_SYMMETRIC_UPPER_TRIANGULAR | SOLVES_SKEW_SYMMETRIC_LOWER_TRIANGULAR)
#define SOLVES_HERMITIAN_SYMMETRIC_UPPER_TRIANGULAR (1<<14)
#define SOLVES_HERMITIAN_SYMMETRIC_LOWER_TRIANGULAR (1<<15)
#define SOLVES_HERMITIAN_SYMMETRIC                  (SOLVES_HERMITIAN_SYMMETRIC_UPPER_TRIANGULAR | SOLVES_HERMITIAN_SYMMETRIC_LOWER_TRIANGULAR)

#define SOLVES_DATA_TYPE_REAL_DOUBLE    (1<<16)
#define SOLVES_DATA_TYPE_REAL_SINGLE    (1<<17)
#define SOLVES_DATA_TYPE_COMPLEX_DOUBLE (1<<18)
#define SOLVES_DATA_TYPE_COMPLEX_SINGLE (1<<19)

#define SOLVES_RHS_BASE_ZERO (1<<20)
#define SOLVES_RHS_BASE_ONE  (1<<21)

#define SOLVES_RHS_DROW (1<<22)
#define SOLVES_RHS_DCOL (1<<23)
#define SOLVES_RHS_COO  (1<<24)
#define SOLVES_RHS_CSC  (1<<25)
#define SOLVES_RHS_CSR  (1<<26)

#define SOLVES_RHS_VECTOR_ONLY (1<<27)


// Note: remember to upate this value when adding more solvers!
#define SOLVER_INDEX_COUNT 2
#define SOLVER_SHORTNAME_MAX_LEN 7

static const struct solver_properties_t solver_lookup[] = {
  { "umfpack", "UMFPACK", "Tim Davis et al", "University of Florida", "5.5.0", "GPL",
    "http://www.cise.ufl.edu/research/sparse/umfpack",
    // TODO 5.5.1 is available
    NULL, NULL, NULL,
    SOLVES_FORMAT_CSC | SOLVES_BASE_ZERO | SOLVES_UNSYMMETRIC |
    SOLVES_DATA_TYPE_REAL_DOUBLE |
    SOLVES_SQUARE | // TODO is umfpack really restricted to square matrices?
    SOLVES_RHS_DCOL | SOLVES_RHS_VECTOR_ONLY,
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
    NULL, NULL, NULL,
    SOLVES_FORMAT_COO | SOLVES_BASE_ONE | SOLVES_UNSYMMETRIC |
    SOLVES_DATA_TYPE_REAL_DOUBLE |
    SOLVES_RHS_DCOL | SOLVES_RHS_VECTOR_ONLY,
    "    [1] P. R. Amestoy, I. S. Duff, J. Koster and J.-Y. L'Excellent,\n"
    "        A fully asynchronous multifrontal solver using distributed dynamic scheduling,\n"
    "        SIAM Journal of Matrix Analysis and Applications, Vol 23, No 1, pp 15-41 (2001).\n"
    "    [2] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and S. Pralet,\n"
    "        Hybrid scheduling for the parallel solution of linear systems.\n"
    "        Parallel Computing Vol 32 (2), pp 136-156 (2006).\n" }
};

// Note: did the SOLVER_INDEX_COUNT get updated when adding more solvers??

#endif
