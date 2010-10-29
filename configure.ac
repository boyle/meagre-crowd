#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.65)
AC_INIT([meagre-crowd],[0.0.1],[alistair.js.boyle@gmail.com])

AC_COPYRIGHT([Copyright (C) 2010 Alistair Boyle])

AC_CONFIG_SRCDIR([src/meagre-crowd.c])
AC_CONFIG_HEADERS([src/config.h])

# makefile
AM_INIT_AUTOMAKE
AC_CONFIG_FILES([Makefile])

# Checks for programs.
AC_PROG_CC
#AC_OPENMP
#AC_C_INLINE

# Checks for libraries.
AC_SEARCH_LIBS([dmumps_c], [dmumps],,AC_MSG_ERROR([could not find required double precision MUMPS library]))
AC_SEARCH_LIBS([load_sparse_matrix], [sparse_matrix_converter],,AC_MSG_ERROR([could not find required BeBOP sparse matrix I/O library]))

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])
AC_HEADER_ASSERT

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_OUTPUT
