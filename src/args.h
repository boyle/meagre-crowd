#ifndef _ARGS_H_
#define _ARGS_H_

#include "config.h"
#include <argp.h>

enum solver_types_t {
  // MPI based solvers
  MUMPS = 0,
  // multithreaded solvers
  // single threaded solvers
  UMFPACK = 10 };

// command line options
struct parse_args {
  char* input;
  char* output;
  unsigned int timing_enabled;
  unsigned int verbosity;
  unsigned int rep;
  int mpi_rank; // id
  enum solver_types_t solver;
};

int parse_args(int argc, char ** argv, struct parse_args* args);

#endif