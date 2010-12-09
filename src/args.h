#ifndef _ARGS_H_
#define _ARGS_H_

#include "config.h"
#include <argp.h>



// command line options
struct parse_args {
  char* input;
  char* output;
  unsigned int timing_enabled;
  unsigned int verbosity;
  unsigned int rep;
  int mpi_rank; // id
};

int parse_args(int argc, char ** argv, struct parse_args* args);

#endif