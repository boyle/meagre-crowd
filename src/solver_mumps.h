#ifndef _SOLVER_MUMPS_H_
#define _SOLVER_MUMPS_H_

#include "config.h"
#include <dmumps_c.h>
#include <bebop/smc/sparse_matrix.h>
#include "args.h"
#include "perftimer.h"


DMUMPS_STRUC_C* solver_init_dmumps(struct parse_args* args, perftimer_t* timer, struct sparse_matrix_t* A);

void solver_data_prep_dmumps(DMUMPS_STRUC_C* id, struct sparse_matrix_t* A, double* b);

void solver_solve_dmumps(DMUMPS_STRUC_C* id, struct parse_args* args, perftimer_t* timer);

void solver_finalize_dmumps(DMUMPS_STRUC_C* id);

#endif