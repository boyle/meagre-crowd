 # Meagre-Crowd: A sparse distributed matrix solver testbench for performance benchmarking.
 # Copyright (C) 2010 Alistair Boyle <alistair.js.boyle@gmail.com>
 #
 #     This file is part of Meagre-Crowd.
 #
 #     Meagre-Crowd program is free software: you can redistribute it and/or modify
 #     it under the terms of the GNU Lesser General Public License as published by
 #     the Free Software Foundation, either version 3 of the License, or
 #     (at your option) any later version.
 #
 #     This program is distributed in the hope that it will be useful,
 #     but WITHOUT ANY WARRANTY; without even the implied warranty of
 #     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 #     GNU Lesser General Public License for more details.
 #
 #     You should have received a copy of the GNU Lesser General Public License
 #     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([meagre-crowd],[0.4.6],[<meagre-crowd@googlegroups.com>])

AC_COPYRIGHT([Copyright (C) 2010, 2011 Alistair Boyle])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])
m4_include([m4/ax_pthread.m4])

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

# Checks for libraries.

# BeBOP: matrix market file format
AC_SEARCH_LIBS([load_sparse_matrix], [sparse_matrix_converter],,AC_MSG_ERROR([missing BeBOP sparse matrix I/O library]))
#  BeBOP: some OSes don't see bebop_util.so properly even though its marked as a dependency of sparse_matrix_converter.so
AC_SEARCH_LIBS([bebop_default_initialize], [bebop_util],,AC_MSG_ERROR([missing BeBOP utils]))

# MPI checks (meagre-crowd usees this too)
AX_MPI()
AC_SUBST(MPILIBS)


AC_ARG_WITH(all-solvers, AS_HELP_STRING(--with-all-solvers, Force presence of all solvers))

m4_define([MC_WITH_LIB],[

dnl ---------------------------------
dnl $1
dnl arguments parsing: handle enabling or disabling various solver libraries
dnl  - behaviour is to try to detect the library, if it works, great... if it fails
m4_define(MC_FULLNAME, $1)
m4_define(MC_LOWERNAME, m4_translit(MC_FULLNAME, [A-Z], [a-z]))
m4_define(MC_UPPERNAME, m4_translit(MC_FULLNAME, [a-z], [A-Z]))
m4_define(MC_FUNC, $2)
m4_define(MC_LIB, $3)
m4_define(MC_EXTRA_CODE, $4)
m4_define(MC_EXTRA_FLAGS, $5)

AC_ARG_WITH(MC_LOWERNAME,
    AS_HELP_STRING([--without-]MC_LOWERNAME, [Ignore presence of ]MC_FULLNAME[ solver and disable it]))

dnl force value for all solvers if --with-all-solvers
AS_IF([test "x$with_all_solvers" = "xyes"],
      [with_]MC_LOWERNAME[=yes],
      [test "x$with_all_solvers" = "xno"],
      [with_]MC_LOWERNAME[=no]
      )

dnl test for library unless explicitly disabled
MC_EXTRA_CODE
AS_IF([test "x$with_]MC_LOWERNAME[" != "xno"],
      [
        AC_SEARCH_LIBS(]MC_FUNC[,]MC_LIB[,have_]MC_LOWERNAME[=yes,have_]MC_LOWERNAME[=no,]MC_EXTRA_FLAGS[)
      ],
      [have_]MC_LOWERNAME[=no
      AC_MSG_CHECKING(for library containing ]MC_FUNC[)
      AC_MSG_RESULT(<skipped>)]
      )

dnl if its found, then all good
dnl if it wasn't found and it was required, then error
AS_IF([test "x$have_]MC_LOWERNAME[" = "xyes"],
      [AC_DEFINE(HAVE_]MC_UPPERNAME[,1,]MC_FULLNAME[ solver is available)],
      [test "x$with_]MC_LOWERNAME[" = "xyes"],
      [AC_MSG_ERROR(]MC_FULLNAME[ solver requested but not found)]
      )

AM_CONDITIONAL([HAVE_]MC_UPPERNAME,[test "x$have_]MC_LOWERNAME[" = "xyes"])
AC_SUBST([have_]MC_LOWERNAME)

m4_undefine([MC_FULLNAME])
m4_undefine([MC_LOWERNAME])
m4_undefine([MC_UPPERNAME])
m4_undefine([MC_FUNC])
m4_undefine([MC_LIB])
m4_undefine([MC_EXTRA_CODE])
m4_undefine([MC_EXTRA_FLAGS])
])

MC_WITH_LIB([MUMPS],   [dmumps_c],      [dmumps])
MC_WITH_LIB([UMFPACK], [umfpack_di_symbolic], [umfpack])
# might not be needed, if already included through another solver
MC_WITH_LIB([CHOLMOD], [cholmod_solve], [cholmod])
# has a 30 day license!
# TODO TAUCS depends on CILK?
MC_WITH_LIB([TAUCS],   [taucs_linsolve],[taucs])
MC_WITH_LIB([SuperLU_DIST], [superlu_gridinit],[superlu_dist])
# TODO (others: gfortran, pthread, m?)     AX_OPENMP([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
# TODO FIXME any MC_WITH_LIBS that come after pardiso is detected will fail if they don't have $OPENMP_CFLAGS since the .so doesn't have the dependency it needs
MC_WITH_LIB([Pardiso], [pardisoinit],   [pardiso], [[AC_OPENMP AX_PTHREAD AC_F77_LIBRARY_LDFLAGS]], $FLIBS $OPENMP_CFLAGS $LIB_PTHREAD)
MC_WITH_LIB([WSMP],    [pwgsmp_],       [pwsmp],   [[AC_OPENMP AX_PTHREAD AC_F77_LIBRARY_LDFLAGS]], $FLIBS $OPENMP_CFLAGS $LIB_PTHREAD)

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h sys/time.h])
AC_HEADER_ASSERT

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
# TODO AC_FUNC_MALLOC -- need replacement rpl_malloc or this breaks...
AC_FUNC_STRNLEN
AC_CHECK_FUNCS([gettimeofday strndup])

# Setup the test suite
AC_CONFIG_TESTDIR([.],[.:tests])
AC_CONFIG_FILES([atlocal])

AC_MSG_CHECKING(for beer in the fridge)
  AC_MSG_RESULT(all systems go)

AC_OUTPUT
