#include "config.h"
#include "file.h"
#include <stdio.h> // fprintf
#include <string.h> // strnlen, strcmp
#include <assert.h>

#include <bebop/util/init.h>
#include <bebop/util/enumerations.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h> // load_sparse_matrix
#include <bebop/smc/coo_matrix.h> // convert to coo

int _identify_format_from_extension(char* n, enum sparse_matrix_file_format_t* ext);

// load a matrix from file "n" into matrix A
// returns 0: success, 1: failure
int load_matrix(char* n, struct sparse_matrix_t** A) {
  if(n == NULL) {
    fprintf(stderr,"input error: No input specified (-i)\n");
    return 1; // failure
  }

  int retval;
  enum sparse_matrix_file_format_t ext;
  if( (retval = _identify_format_from_extension(n, &ext)) != 0)
    return retval;

  *A = load_sparse_matrix(ext, n);
  if(*A == NULL) {
    fprintf(stderr,"input error: Failed to load matrix\n");
    return 1;
  }
  return 0; // success
}

// returns number of rows in matrix A
unsigned int matrix_rows(struct sparse_matrix_t* A) {
  unsigned int r = 0;
  if(A == NULL)
    return 0;
  int ierr = sparse_matrix_convert(A, COO); assert(ierr == 0);
  struct coo_matrix_t* Acoo = A->repr;
  return Acoo->m;
}
    // identify the file format from the extension
    // return 0: success, 1: failure
int _identify_format_from_extension(char* n, enum sparse_matrix_file_format_t* ext) {

  size_t s = strnlen(n,100);
  if(s > 3) {
    char *e = n + s - 3;
    // strcmp returned match
    if(strcmp(e,".mm") == 0) {
      *ext = MATRIX_MARKET;
      return 0; // success
    }
    else if(strcmp(e,".hb") == 0) {
      *ext = HARWELL_BOEING;
      fprintf(stderr,"input error: Sorry Harwell-Boeing reader is broken\n");
      return 1; // failure
    }
    else if(strcmp(e,".rb") == 0) {
      *ext = HARWELL_BOEING;
      fprintf(stderr,"input error: Sorry Rutherford-Boeing reader is broken\n");
      return 1; // failure
    }
    else if((s > 4) && (strcmp(e-1,".mat") == 0)) {
      *ext = MATLAB;
      fprintf(stderr,"input error: Sorry Matlab reader is broken\n");
      return 1; // failure
    } // TODO test if the matlab reader is actually busted
    else{
      fprintf(stderr,"input error: Unrecognized file extension\n");
      return 1; // failure
    }
  }
}
