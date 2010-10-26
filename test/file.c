/* C Example */
#include <stdio.h>
#include <assert.h>

#include <bebop/util/init.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h>

int main (argc, argv)
     int argc;
     char *argv[];
{
  /* 
  * Do default library initialization.  In particular, this reads the
  * BEBOP_DEBUG_LEVEL environment variable.  If its value is 0, no 
  * debug output is produced, and increasing integral values cause more
  * and more debug output to be produced.  Debug output goes to stderr;
  * you can change this by using a nondefault initialization seqeuence.
  * The documentation gives details.
  */
  int errcode = 0;
  bebop_default_initialize (argc, argv, &errcode);
  assert (errcode == 0);


  const char* out = "test-out.mm";
  //struct sparse_matrix_t* A = load_sparse_matrix (HARWELL_BOEING, "test.hb");
  //struct sparse_matrix_t* A = load_sparse_matrix (HARWELL_BOEING, "test.rb");
  struct sparse_matrix_t* A = load_sparse_matrix (MATRIX_MARKET, "test.mm");
  assert (A != NULL);
  save_sparse_matrix (out, A, MATRIX_MARKET);
  destroy_sparse_matrix (A);

  bebop_exit(errcode);
  assert(errcode == 0);
  return 0;
}
