AC_DEFUN([ACX_MPI], [
  AC_ARG_VAR(MPICC,[MPI C compiler command])
  AC_CHECK_PROGS(MPICC, mpicc, $CC)
  acx_mpi_save_CC="$CC"
  CC="$MPICC"
  AC_SUBST(MPICC)
])
