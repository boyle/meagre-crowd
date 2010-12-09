#ifndef _FILE_H_
#define _FILE_H_

#include "config.h"
#include <bebop/smc/sparse_matrix.h>

// load a matrix from file "n" into matrix A
// returns 0: success, 1: failure
int load_matrix(char* n, struct sparse_matrix_t** A);

// returns number of rows in matrix A
unsigned int matrix_rows(struct sparse_matrix_t* A);

#endif
