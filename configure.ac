#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.65)
AC_INIT([meagre-crowd],[0.0.5],[<alistair.js.boyle@gmail.com>])

AC_COPYRIGHT([Copyright (C) 2010 Alistair Boyle])

AC_CONFIG_AUX_DIR([build-aux])

AC_CONFIG_SRCDIR([src/meagre-crowd.c])
AC_CONFIG_HEADERS([src/config.h])

AC_USE_SYSTEM_EXTENSIONS



# debugging

AC_MSG_CHECKING(whether to enable debugging)
default="no"
DEBUG=
AC_ARG_ENABLE(debug, [  --enable-debug=[no/yes] turn on debugging
                       [default=$default]],, enable_debug=$default)
if test "x$enable_debug" = "xyes"; then
  CFLAGS="$CFLAGS -g -DDEBUG -Wall"
  CXXFLAGS="$CXXFLAGS -g -DDEBUG -Wall "
  DEBUG="${DEBUG} -g";
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi
# profiling
AC_MSG_CHECKING(whether to enable profiling)
default="no"
AC_ARG_ENABLE(debug, [  --enable-profile=[no/yes] turn on profiling
                       [default=$default]],, enable_profile=$default)
if test "x$enable_profile" = "xyes"; then
  CFLAGS="${CFLAGS} -g -DDEBUG -pg"
  CXXFLAGS="${CXXFLAGS} -g -DDEBUG -pg"
  DEBUG="${DEBUG} -pg";
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi
# .. and export so Makefile will use them
AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(DEBUG)


# makefile
AM_INIT_AUTOMAKE([1.11 subdir-objects]) # need at least v1.11 of automake for SILENT_RULES
AM_SILENT_RULES([yes]) # enable kernel like cleaner build
AC_CONFIG_FILES([Makefile])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
#AC_OPENMP

# Checks for libraries.
AC_SEARCH_LIBS([dmumps_c], [dmumps],,AC_MSG_ERROR([could not find required double precision MUMPS library]))
AC_SEARCH_LIBS([load_sparse_matrix], [sparse_matrix_converter],,AC_MSG_ERROR([could not find required BeBOP sparse matrix I/O library]))

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h])
AC_HEADER_ASSERT

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_STRNLEN
AC_CHECK_FUNCS([gettimeofday strndup])

# Setup the test suite
AC_CONFIG_TESTDIR([.],[.:tests])
AC_CONFIG_FILES([atlocal])

AC_MSG_CHECKING(for beer in the fridge)
  AC_MSG_RESULT(all systems go)

AC_OUTPUT
