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

  struct sparse_matrix_t* A;

  printf("harwell-boeing\n");
  A = load_sparse_matrix (HARWELL_BOEING, "test.hb");
  assert (A != NULL);
  //save_sparse_matrix ("out-test.hb2mm", A, MATRIX_MARKET); // TODO broken output (all zeros)
  //save_sparse_matrix ("out-test.hb2hb", A, HARWELL_BOEING); TODO broken output (all zeros)
  destroy_sparse_matrix (A);
  
  printf("rutherford-boeing\n");
  A = load_sparse_matrix (HARWELL_BOEING, "test.rb");
  assert (A != NULL);
  //save_sparse_matrix ("out-test.rb2mm", A, MATRIX_MARKET); // TODO segfault
  //save_sparse_matrix ("out-test.rb2hb", A, HARWELL_BOEING); // TODO broken output (all zeros)
  destroy_sparse_matrix (A);

  printf("matrix market\n");
  A = load_sparse_matrix (MATRIX_MARKET, "test.mm");
  assert (A != NULL);
  save_sparse_matrix ("out-test.mm2mm", A, MATRIX_MARKET);
  save_sparse_matrix ("out-test.mm2hb", A, HARWELL_BOEING);
  destroy_sparse_matrix (A);
  
  // clean up
  printf("PASS\n");
  return 0;
}
