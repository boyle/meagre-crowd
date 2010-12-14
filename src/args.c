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
#include "config.h"
#include "args.h"
#include <stdlib.h>
#include <string.h>

// global argp variables
// TODO somehow hide these?
const char* argp_program_version = PACKAGE_STRING;
const char* argp_program_bug_address = PACKAGE_BUGREPORT;

error_t parse_opt(int key, char *arg, struct argp_state *state);
error_t parse_opt(int key, char *arg, struct argp_state *state) {
  struct parse_args *args = state->input;
  switch (key) {
    // help, etc
    case 'h': case'?': argp_state_help(state,state->out_stream,ARGP_HELP_STD_HELP); break; // help
    case -1: argp_state_help(state,state->out_stream,ARGP_HELP_SHORT_USAGE | ARGP_HELP_EXIT_OK); break; // usage
    case 'V': printf("%s\n",PACKAGE_STRING); exit(EXIT_SUCCESS); break; // version
    // output controls
    case 't': args->timing_enabled++; break;
    case 'v': args->verbosity++; break;

    // solvers
    case 's': {
        if(strncmp(arg,"mumps",6)==0)
          args->solver = MUMPS;
        else if(strncmp(arg,"umfpack",8)==0)
          args->solver = UMFPACK;
        else {
          fprintf(stderr,"invalid solver (-s)");
          exit(EXIT_FAILURE);
        }
      }
      break;

    case 'r': {
        int i = atoi(arg);
        i = (i < 0)?0:i; // > 0
        args->rep = i;
      }
      break;
    // file I/O
    case 'i':
      args->input = arg;
      FILE* f = fopen(args->input, "r");
      if(f == NULL) {
        perror("input error");
        exit(EXIT_FAILURE);
      }
      fclose(f);
      break;

    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}


int parse_args(int argc, char ** argv, struct parse_args* args) {

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
      {"input", 'i',"FILE",0,"Input matrix FILE",10},
      // TODO these should just be all the other command line components (no arg required ala gcc)
      {"verbose",'v',0,0,"Increase verbosity",11},
      // TODO add note to man page: -v, -vv, -vvv, etc for more detail
      // none: no output, -v: matrix info & any available stats (i.e. cond. number),
      // -vv: more detail(?), -vvv: max debug
      {"timing",'t',0,0,"Show/increase timing information",12},
      {"repeat",'r',"N",0,"Repeat calculations N times",15},
      {"solver",'s',"<mumps|umfpack>",0,"Select solver",15},
      // TODO add note to man page: -t, -tt, -ttt for more detail
      // none: no output, -t: single-line (csv), -tt: chart, -ttt: greater detail
      //{"output",'o',0,0,"Output file (default: stdout)",20},
      { 0 } // null termintated list
    };
    // argp_option*, argp_parser, extra-usage line options, pre-help, // optional: argp_child, *help_filter, argp_domain
    const struct argp p = { opt, parse_opt, 0, doc};
    if(argc == 1) { // there's no arguments, exit!
      char *prog = strndup(argv[0],100);
      argp_help(&p, stderr, ARGP_HELP_SHORT_USAGE, basename(prog));
      free(prog);
      return EXIT_FAILURE; // so, exit for reals
    }
    // error_t argp_parse (argp*, argc, **argv, unsigned flags, int *arg_index, void *input)
    argp_parse(&p, argc, argv, ARGP_NO_HELP, 0, args);
  }

  return EXIT_SUCCESS;
}
