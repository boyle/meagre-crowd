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
#ifndef _ARGS_H_
#define _ARGS_H_

#include "config.h"
#include <argp.h>
#include "solvers.h"

// command line options
struct parse_args {
  char* input;
  char* output;
  char* rhs;
  char* expected;
  double expected_precision;
  unsigned int timing_enabled;
  unsigned int verbosity;
  unsigned int rep;
  int mpi_rank; // id
  int solver;
};

int parse_args( int argc, char ** argv, struct parse_args* args );

#endif
