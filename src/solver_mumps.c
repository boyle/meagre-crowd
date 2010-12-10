#include "solver_mumps.h"


#include <stdlib.h> // malloc, free
#include <string.h> // memcpy
#include <assert.h>

#include <bebop/smc/coo_matrix.h>
#include <bebop/util/enumerations.h>


// mumps job controls
#define JOB_INIT -1
#define JOB_END -2

#define MUMPS_USE_COMM_WORLD -987654



// Initialize a MUMPS instance. Use MPI_COMM_WORLD.
// Note: only using A to determine matrix type
DMUMPS_STRUC_C* solver_init_dmumps(struct parse_args* args, perftimer_t* timer, struct sparse_matrix_t* A) {
  // initialize MUMPS instance
  DMUMPS_STRUC_C* id = calloc(1,sizeof(DMUMPS_STRUC_C)); // initialize to zero
  assert(id != NULL); // calloc failure
  id->job=JOB_INIT;
  id->par=1; // host involved in factorization/solve
  id->sym=0; // 0: general, 1: sym pos def, 2: sym (note: no hermitian support) // TODO support other matrix types
  // Note: if set to symmetric and matrix ISN'T, the redundant entries will be *summed*
  // TODO could convert the C communicator instead of using the fortran one (see MUMPS doc)
  id->comm_fortran=MUMPS_USE_COMM_WORLD;
  #define INFOG(I) infog[(I)-1] // macro s.t. indices match documentation
  dmumps_c(id);  assert(id->INFOG(1) == 0); // check it worked
  // clears the rest of the unused values

  // set debug verbosity
  #define ICNTL(I) icntl[(I)-1] // macro s.t. indices match documentation
  // No outputs
  if(args->verbosity < 3) { // no debug
    id->ICNTL(1)=-1; id->ICNTL(2)=-1; id->ICNTL(3)=-1; id->ICNTL(4)=0;
  }
  else { // debug
    id->ICNTL(1)=6; // err output stream
    id->ICNTL(2)=6; // warn/info output stream
    id->ICNTL(3)=6; // global output stream
    id->ICNTL(4)=4; // debug level 0:none, 1: err, 2: warn/stats 3:diagnostics, 4:parameters
  }
  return id;
}

// if NULL, do nothing to A or b
// TODO b can be sparse...
void solver_data_prep_dmumps(DMUMPS_STRUC_C* id, struct sparse_matrix_t* A, double* b) {
  if(A != NULL) {
    // prepare the matrix
    int ierr = sparse_matrix_convert(A, COO); assert(ierr == 0);
    struct coo_matrix_t* Acoo = A->repr;
    coo_c_to_fortran(Acoo); assert(Acoo != NULL);
    assert(Acoo->index_base == ONE); // index zero is the first entry
    assert(Acoo->symmetry_type == UNSYMMETRIC);
    assert(Acoo->value_type == REAL); // don't handle complex... yet TODO

    // TODO do soemthing with A.ownership, so we can tell bebop to clean itself up, but not have to copy the elements
    // TODO really we should just copy this to be CORRECT/TYPESAFE (not worth being clever...)

    // mumps: irn=row indices, jcn=column idices, a=values, rhs=right-hand side, n = matrix order (on-a-side?) nz=non-zeros?
    id->n   = Acoo->m; // A.m: rows, A.n: columns
    id->nz  = Acoo->nnz; // non-zeros
    id->irn = Acoo->II; // row    indices
    id->jcn = Acoo->JJ; // column indices
    id->a   = Acoo->val;
  }

  if(b != NULL) {
    free(id->rhs); // no-op if NULL
    id->rhs = malloc(id->n * sizeof(double));
    assert(id->rhs != NULL); // malloc failure
    memcpy(id->rhs, b, id->n);
  }
}

void solver_solve_dmumps(DMUMPS_STRUC_C* id, struct parse_args* args, perftimer_t* timer) {
  // TODO rm?  perftimer_inc(timer,"solver",-1);
  // TODO rm?  perftimer_adjust_depth(timer,+1);

  // Call the MUMPS package.
  perftimer_inc(timer,"analyze",-1);
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
  id->job=JOB_ANALYSE;
  dmumps_c(id);
  if(id->INFOG(1) != 0) fprintf(stderr, "warning: analysis failed\n");
  assert(id->INFOG(1) == 0); // check it worked

  // available info:
  // INFO(15)/INFOG(16/17): min/max/sum-over-all-cpus mem requried [in megabytes]
  // INFO(17): min mem for out-of-core (max, sum in INFOG(26,27))
  //   set ICNTL(23) for explicit max mem, per-proc [MB]
  //   set ICNTL(14) to limit mem increases

  perftimer_inc(timer,"factorize",-1);
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
  id->job=JOB_FACTORIZE;
  dmumps_c(id);
  assert(id->INFOG(1) == 0); // check it worked

  perftimer_inc(timer,"solve",-1);
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
  id->job=JOB_SOLVE;
  dmumps_c(id);
  assert(id->INFOG(1) == 0); // check it worked
  perftimer_inc(timer,"done",-1);
}


void solver_finalize_dmumps(DMUMPS_STRUC_C* id) {
  id->job=JOB_END;
  dmumps_c(id); // Terminate instance
  assert(id->INFOG(1) == 0); // check it worked

  free(id->rhs);
  free(id);
}
