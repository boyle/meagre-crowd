/* Example program using the C interface to the
* double precision version of MUMPS, dmumps_c.
* We solve the system A x = RHS with
* A = diag(1 2) and RHS = [1 4]ˆT
* Solution is [1 2]ˆT */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <argp.h>
#include <libgen.h>
#include <string.h>

#include <mpi.h>
#include <dmumps_c.h>

#include "config.h"

#include <bebop/util/init.h>
#include <bebop/util/enumerations.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h>
#include <bebop/smc/coo_matrix.h>

#include "perftimer.h"

// mumps job controls
#define JOB_INIT -1
#define JOB_END -2

// global argp variables
// TODO somehow hide these?
const char* argp_program_version = PACKAGE_STRING;
const char* argp_program_bug_address = PACKAGE_BUGREPORT;
// command line option
struct parse_args {
  char* input;
  char* output;
  unsigned int timing_enabled;
};

error_t parse_opt(int key, char *arg, struct argp_state *state);
error_t parse_opt(int key, char *arg, struct argp_state *state) {
  struct parse_args *args = state->input;
  switch (key) {
    case 'p': args->timing_enabled++; break;
    case 'h': case'?': argp_state_help(state,state->out_stream,ARGP_HELP_STD_HELP); break; // help
    case -1: argp_state_help(state,state->out_stream,ARGP_HELP_SHORT_USAGE | ARGP_HELP_EXIT_OK); break; // usage
    case 'V': printf("%s\n",PACKAGE_STRING); exit(0); break; // version

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}



// test result of matrix computations
int results_match(const double * const expected, const double const * result, const int n, const double precision);
int results_match(const double * const expected, const double const * result, const int n, const double precision) {
  int i;
  for(i=0;i<n;i++) {
    if( (result[i] < (expected[i] - precision)) ||
        (result[i] > (expected[i] + precision)) ) {
      return 0;
    }
  }
  return 1;
}


#define USE_COMM_WORLD -987654
int main(int argc, char ** argv) {
  int myid, ierr;
  int retval = 0;

  perftimer_t* timer = perftimer_malloc();
  perftimer_inc(timer,"initialization",-1);
  perftimer_adjust_depth(timer,+1);

  perftimer_inc(timer,"command-line parsing",-1);
  struct parse_args args = {0};

  // parse command line
  {
    static char doc[] = "\n\
Solves Ax=b for x, where A is sparse.\n\
Given an input matrix A and right-hand side b, find x.\n\
\n\
This is intended as a generic performance measurement platform\n\
providing timing and other metrics for distributed matrix solvers.\n\
This includes legacy single threaded solvers for benchmarking,\n\
Shared Memory solvers (SMP, generally OpenMP based), and\n\
heterogeneous solvers (generally using SMP).\n\
\n\
Matrices are available through the Harwell-Boeing sparse matrix\n\
collection and the University of Florida sparse matrix collection.\n\
\n\
Limitations: currently only 'Matrix Market' format is supported (*.mm).\n\
  (The Rutherford/Harwell-Boeing format loader is broken. *.hb, *.rb)\n\
\n\
Options:"
;
    // TODO describe fields..
    // "long", 'l', "value", flags, "desc", groupid
    static const struct argp_option opt[] = {
      {"help", 'h',0,0,"Give this help list"},
      {0,      '?',0,OPTION_ALIAS},
      {"usage",-1, 0,0,"Show usage information",-1},
      {"version",'V', 0,0,"Show version information",-1},
      {"input", 'i',0,0,"Input matrix file",10},
      // TODO these should just be all the other command line components (no arg required ala gcc)
      {"perf",'p',0,0,"Show performance/timing information",11},
      // TODO add note to man page: -p, -pp, -ppp, etc for more detail
      // TODO option to change output format "--perf-csv"
      {"output",'o',0,0,"Output file (default: stdout)",20},
      { 0 } // null termintated list
    };
    // argp_option*, argp_parser, extra-usage line options, pre-help, // optional: argp_child, *help_filter, argp_domain
    const struct argp p = { opt, parse_opt, 0, doc};
    if(argc == 1) { // there's no arguments, exit!
      char *prog = strndup(argv[0],100);
      argp_help(&p, stderr, ARGP_HELP_SHORT_USAGE, basename(prog));
      free(prog);
      perftimer_free(timer);
      return EXIT_FAILURE; // so, exit for reals
    }
    // error_t argp_parse (argp*, argc, **argv, unsigned flags, int *arg_index, void *input)
    argp_parse(&p, argc, argv, ARGP_NO_HELP, 0, &args);
  }

  // initialize MPI
  perftimer_inc(timer,"MPI init",-1);
  ierr = MPI_Init(&argc, &argv);               assert(ierr == 0);
  perftimer_inc(timer,"MPI comms",-1);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid); assert(ierr == 0);

  // Initialize a MUMPS instance. Use MPI_COMM_WORLD.
  perftimer_inc(timer,"MUMPS",-1);
  DMUMPS_STRUC_C id;
  id.job=JOB_INIT;
  id.par=1; // host involved in factorization/solve
  id.sym=0; // 0: general, 1: sym pos def, 2: sym (note: no hermitian support)
  // Note: if set to symmetric and matrix ISN'T, the redundant entries will be *summed*
  // TODO could convert the C communicator instead of using the fortran one (see MUMPS doc)
  id.comm_fortran=USE_COMM_WORLD;
  #define INFOG(I) infog[(I)-1] // macro s.t. indices match documentation
  dmumps_c(&id);  assert(id.INFOG(1) == 0); // check it worked
  // clears the rest of the unused values

  // Define the problem on the host
  if (myid == 0) {
    perftimer_inc(timer,"sparse file handling",-1);
    bebop_default_initialize (argc, argv, &ierr); assert (ierr == 0);

    perftimer_adjust_depth(timer,-1);
    perftimer_inc(timer,"input",-1);
    perftimer_adjust_depth(timer,+1);

    struct sparse_matrix_t* A;
    perftimer_inc(timer,"load",-1);
    A = load_sparse_matrix (MATRIX_MARKET, "test.mm"); assert (A != NULL);
    ierr = sparse_matrix_convert(A, COO); assert(ierr == 0);

    // check somethings are as expected
    perftimer_inc(timer,"sanity check",-1);
    struct coo_matrix_t* Acoo = A->repr;
    coo_c_to_fortran(Acoo); assert(Acoo != NULL);
    assert(Acoo->index_base == ONE); // index zero is the first entry
    assert(Acoo->symmetry_type == UNSYMMETRIC);
    assert(Acoo->value_type == REAL); // don't handle complex... yet TODO

    // TODO do soemthing with A.ownership, so we can tell bebop to clean itself up, but not have to copy the elements

    perftimer_inc(timer,"reformat",-1);
    // mumps: irn=row indices, jcn=column idices, a=values, rhs=right-hand side, n = matrix order (on-a-side?) nz=non-zeros?
    id.n   = Acoo->m; // A.m: rows, A.n: columns
    id.nz  = Acoo->nnz; // non-zeros
    id.irn = Acoo->II; // row    indices
    id.jcn = Acoo->JJ; // column indices
    id.a   = Acoo->val;
    // show what we got
    printf("Input A is\n");
    int i;
    for(i=0;i<id.nz;i++) {
      printf("  (%i,%i)=%.2f\n",id.irn[i],id.jcn[i],id.a[i]);
    }
    perftimer_inc(timer,"rhs",-1);
    // allocate an all-zeros right-hand side of A.m rows
    id.rhs = malloc(id.n * sizeof(double)); assert(id.rhs != NULL);
    printf("Right-hand side is\n");
    for(i=0;i<id.n;i++) {
      id.rhs[i] = i;
      printf("  %.2f\n", id.rhs[i]);
    }
    perftimer_adjust_depth(timer,-1);
    //destroy_sparse_matrix (A); // TODO can't release it unless we're copying it...
  }

  perftimer_inc(timer,"solve",-1);
  perftimer_adjust_depth(timer,+1);
  perftimer_inc(timer,"MUMPS config",-1);
  #define ICNTL(I) icntl[(I)-1] // macro s.t. indices match documentation
  // No outputs
  if(1) { // no debug
    id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
  }
  else { // debug
    id.ICNTL(1)=6; // err output stream
    id.ICNTL(2)=6; // warn/info output stream
    id.ICNTL(3)=6; // global output stream
    id.ICNTL(4)=4; // debug level 0:none, 1: err, 2: warn/stats 3:diagnostics, 4:parameters
  }
  // Call the MUMPS package.
  perftimer_inc(timer,"MUMPS analyze",-1);
  // ICNTL(22) != 0: out-of-core
  //   requires OOC_TMPDIR, OOC_PREFIX -> tmp location/prefix
  // ICNTL(14): memory relaxation
  // ICNTL(5) = ICNTL(18) = 0; // centralized, assembled matrix load
  // set ICNTL(7)=1 for external ordering, requires id.PERM_IN
  //   0: AMD
  //   1: user provided, id.PERM_IN
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
  // set ICNTL(6)=5,6 for scaling (needs values in id.A)
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
  //   id.N, NZ, IRN, JCN (assembled matrix, ICNTL(5)=0 & ICNTL(18) != 0) or
  //     (where ICNTL(18)=1 or 2 -- distributed, IRN, JCN not required)???
  //     (where ICNTL(18)=3 -- distributed, N on host, and NZ_loc,
  //                           IRN_loc, JCN_loc on slaves)???
  //   id.N, NELT, ELTPTR, ELTVAR (elemental matrix, ICNTL(5)=1)
  //   ICNTL(19)=1,2,3, ICNTL(26) schur complement w/ reduced or condensed rhs
  //     requries id.SIZE_SCHUR, LISTVAR_SCHUR
  //     ICNTL(19)=1: centralized schur, 2: distributed lower schur (sym), or
  //               3: distributed complete schur
  //       ICNTL(19)=2,3 requires id.NPROW, NPCOL, MBLOCK, NBLOCK - processing conf
  //       ... sets SCHUR_MLOC, SCHUR_NLOC
  // id.WRITE_PROBLEM: store distributed in matrix market format
  #define JOB_ANALYSE 1
  id.job=JOB_ANALYSE;
  dmumps_c(&id);
  assert(id.INFOG(1) == 0); // check it worked

  // available info:
  // INFO(15)/INFOG(16/17): min/max/sum-over-all-cpus mem requried [in megabytes]
  // INFO(17): min mem for out-of-core (max, sum in INFOG(26,27))
  //   set ICNTL(23) for explicit max mem, per-proc [MB]
  //   set ICNTL(14) to limit mem increases

  perftimer_inc(timer,"MUMPS factorize",-1);
  // requires id.A if ICNTL(5)=0 (assembled matrix)
  // requires id.A_ELT if ICNTL(5)=1 (elemental matrix)
  // if (ICNTL(5)=0 && ICNTL(18)!=0) (assembled matrix, distributed load) ???i
  //   requires A_loc on slaves
  //   requires NZ_loc, IRN_loc, JCN_loc if ICNTL(18)=1 or 2 ???
  //       -- already passed in if ICNTL(18)=3
  // ICNTL(8)!=0 --> user supplied row/column scaling
  //   requires id.COLSCA, ROWSCA
  // ICNTL(19)=2,3 requires SCHUR_LLD, SCHUR
  #define JOB_FACTORIZE 2
  id.job=JOB_FACTORIZE;
  dmumps_c(&id);
  assert(id.INFOG(1) == 0); // check it worked

  perftimer_inc(timer,"MUMPS solve",-1);
  // solve Ax=b, AX=B OR A^t x=b, A^t X=B
  //   ICNTL(9)=1: A (default), 0: A^t
  // OR compute "null-space basis" if "null pivot row detection" was enabled
  //   ICNTL(24)=1 and INFOG(28)!=0
  // requires id.RHS
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
  // id.NRHS is number of rhs (optional), >1 disables itr refinement, err analysis
  // id.LRHS >= NRHS, =leading dimension of RHS (optional)
  //
  #define JOB_SOLVE 3
  id.job=JOB_SOLVE;
  dmumps_c(&id);
  assert(id.INFOG(1) == 0); // check it worked

  perftimer_inc(timer,"MUMPS clean up",-1);
  id.job=JOB_END; dmumps_c(&id); // Terminate instance
  assert(id.INFOG(1) == 0); // check it worked


  perftimer_adjust_depth(timer,-1);
  perftimer_inc(timer,"output",-1);
  if (myid == 0) {
    printf("Solution is\n");
    int i;
    for(i=0;i<id.n;i++) {
      printf("  %.2f\n", id.rhs[i]);
      // Octave says the answer should be
      // ans =
      //	   1.000000
      //	   0.434211
      //	   0.250000
      //	  -0.092105
      //	   0.250000
      //
    }
    // test result
    // TODO make this a command-line option with file to compare
    perftimer_inc(timer,"test",-1);
    const double const e[5] = {1.0, 0.434211, 0.25, -0.092105, 0.25};
    if( results_match(e,id.rhs,id.n,0.001) ) {
      retval = 0; printf("PASS\n");
    }
    else {
      retval = 1; printf("FAIL\n");
    }

    // clean up
    free(id.rhs);
  }
  perftimer_inc(timer,"clean up",-1);
  perftimer_adjust_depth(timer,+1);

  perftimer_inc(timer,"MPI",-1);
  ierr = MPI_Finalize(); assert(ierr == 0);

  // show timing info, if requested, to depth N
  perftimer_adjust_depth(timer,-1);
  perftimer_inc(timer,"finished",-1);
  if(args.timing_enabled != 0)
    perftimer_printf(timer,args.timing_enabled-1);
  perftimer_free(timer);
  return retval;
}
