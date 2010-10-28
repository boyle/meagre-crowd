#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([meagre-crowd], [0.0.1], [alistair.js.boyle@gmail.com])

AC_COPYRIGHT([Copyright Alistair Boyle (C) 2010])

AC_CONFIG_SRCDIR([src/meagre-crowd.c])
AC_CONFIG_HEADERS([src/config.h])

# makefile
AM_INIT_AUTOMAKE
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([stamp-h], [echo timestamp > stamp-h])  # rebuild when src changes

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])
AC_HEADER_ASSERT

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_OUTPUT
