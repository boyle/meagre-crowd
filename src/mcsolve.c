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
/*#include "config.h"*/

#include "mex.h"

#include <assert.h>
#include <string.h> // memcpy()
#include <stdlib.h> // malloc, free, etc.
#include "solvers.h"
#include "matrix.h"

#define MATLAB_BASE_INDEX   FIRST_INDEX_ZERO
// TODO check MatLab base indexing

/* MatLab gateway routine
 * (Always has the same interface signature. Function name comes from
 * compiled .mex file name.)
 * Matlab function is:
 *   [X] = mcsolve(A,B) --> X = A \ B --> X = A^-1 B   (from AX=B)
 *   mcsolver() will show list of solvers  */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if(nrhs == 0) {
    printf_solvers(1); // TODO get a string and print with mexPrintf - so it works with Windows
    return;
  }

  if(nrhs < 2 || nrhs > 3) {
    mexErrMsgTxt("Two or three input arguments required.");
  }
  else if(nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.");
  }
  else if(!mxIsDouble(prhs[0]) || mxIsSparse(prhs[0]) || mxIsComplex(prhs[0])) {
    mexErrMsgTxt("Expected A to be a real dense matrix.");
  }
  else if(!mxIsDouble(prhs[1]) || mxIsSparse(prhs[1]) || mxIsComplex(prhs[1])) {
    mexErrMsgTxt("Expected B to be a real dense matrix.");
  }
  else if(nrhs == 3 && !mxIsChar(prhs[2])) {
    mexErrMsgTxt("Expected solver string for 3rd parameter, run mcsolve() for a list.");
  }

  assert(sizeof(mwSize) == sizeof(int));

  // make copies of A and b just to be safe
  // TODO look through code and remove this after checking it will always be safe (i.e. copied before changed)
  matrix_t A = {0};
  matrix_t b = {0};
  A.base = MATLAB_BASE_INDEX;
  A.format = DCOL;
  A.sym = SM_UNSYMMETRIC;
  A.data_type = SM_REAL;

  A.m = mxGetM(prhs[0]);
  A.n = mxGetN(prhs[0]);
  A.nz = A.m * A.n;
  A.dd = malloc(A.nz * sizeof(double));
  if(A.dd == NULL) {
    mexErrMsgTxt("Out of memory.");
  }
  memcpy(A.dd, mxGetPr(prhs[0]), A.nz * sizeof(double));
  
  b = A;
  b.m = mxGetM(prhs[1]);
  b.n = mxGetN(prhs[1]);
  b.nz = b.m * b.n;
  b.dd = malloc(b.nz * sizeof(double));
  if(b.dd == NULL) {
    mexErrMsgTxt("Out of memory.");
  }
  memcpy(b.dd, mxGetPr(prhs[1]), b.nz * sizeof(double));

  int solver = 0; // default, TODO select 3rd parameter (string)
  if(nrhs == 3) {
    char * s = mxArrayToString(prhs[2]);
    solver = lookup_solver_by_shortname(s);
    mxFree(s);
  }
  mexPrintf("A is %dx%d, B is %dx%d with %s\n", A.m, A.n, b.m, b.n, solver2str(solver));

solver = 0; // TODO rm (only UMFPACK works until MPI loader is working)
mexPrintf("TODO Only solving with UMFPACK so far...\n");

  // TODO initialize MPI if required - can we keep them alive between calls?
  // TODO launch forks for slaves?
  
  // We'd like to create x's data using mxCreateArray ahead of time
  // so that we don't have to do a memcpy after the solve but we might
  // lose the pointer through conversions and we use malloc/free inside
  matrix_t x = {0}; 

  // run solver
  const int verbosity = 0; // silent
  const int mpi_rank = 0; // TODO
  mc_solver(solver, verbosity, mpi_rank, &A, &b, &x); 

  // TODO close down mpi if required
  // TODO end forked slaves? or prep for next call

  // store x into a format that MatLab likes
  convert_matrix(&x, DCOL, MATLAB_BASE_INDEX);
  plhs[0] = mxCreateDoubleMatrix(x.m, x.n, mxREAL);
  double * const d = mxGetPr(plhs[0]);
  memcpy(d, x.dd, x.m * x.n * sizeof(double));

  clear_matrix(&A);
  clear_matrix(&b);
  clear_matrix(&x);
}
