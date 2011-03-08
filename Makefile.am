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

# build to GNU standards
AUTOMAKE_OPTIONS = gnu


# the actual program to be installed at the end
bin_PROGRAMS = meagre-crowd
meagre_crowd_SOURCES = src/meagre-crowd.c \
    src/perftimer.c src/matrix.c src/args.c src/file.c \
    src/solvers.c \
    src/solver_mumps.c \
    src/solver_umfpack.c \
    src/solver_cholmod.c \
    src/solver_taucs.c \
    src/solver_pardiso.c

# TODO make openMP optional... OPENMP_CFLAGS, PTHREADS and FLIBS are all required only for Paradiso!
AM_LIBS = $(PTHREAD_LIBS) $(LIBS) $(FLIBS)
AM_CFLAGS = $(CFLAGS) $(OPENMP_CFLAGS) $(PTHREAD_CFLAGS) -Werror #-Wextra
# this is required for Sun systems, but will break when we want to do the same for MPI... (mpicc)
#meagre_crowd_CC = $(PTHREAD_CC)

# header files that need to be distributed but not installed
noinst_HEADERS = src/perftimer.h src/matrix.h src/args.h src/file.h \
                 src/solvers.h src/solver_lookup.h \
		 src/solver_mumps.h src/solver_umfpack.h src/solver_cholmod.h src/solver_pardiso.h src/solver_taucs.h

# want these scripts installed in the end
bin_SCRIPTS = launch_jobs collect_results
# some random files we want to distribute with the package tar
EXTRA_DIST = get_dependencies $(bin_SCRIPTS) \
	     tests/unsym.mm tests/sym.mm tests/sym_posdef.mm \
	     tests/rhs1.mm \
	     tests/unsym-default-ans.mm tests/sym-default-ans.mm tests/sym_posdef-default-ans.mm \
	     tests/unsym-rhs1-ans.mm tests/sym-rhs1-ans.mm tests/sym_posdef-rhs1-ans.mm \
	     tests/test.hb tests/test.rb \
             ChangeLog README NEWS INSTALL LICENSE COPYING COPYING.LESSER

# the tests
BUILT_TESTS = tests/file tests/helloworld-mpi tests/helloworld-mumps tests/hello-pardiso tests/mumps-file tests/unit-perftimer tests/unit-matrix tests/hello-umfpack
#TESTS = $(BUILT_TESTS)
check_PROGRAMS = $(BUILT_TESTS) 
tests_unit_perftimer_SOURCES = tests/unit-perftimer.c src/perftimer.c
tests_unit_perftimer_CPPFLAGS = -I$(srcdir)/src
tests_unit_matrix_SOURCES = tests/unit-matrix.c src/matrix.c
tests_unit_matrix_CPPFLAGS = -I$(srcdir)/src

clean-local: clean-local-check
maintainer-clean-local:
	-rm -f COPYING COPYING.LESSER aclocal.m4 configure ChangeLog Makefile.in src/config.h.in*
	-rm -rf build-aux .deps

# regenerate the ChangeLog from git's history via ./autogen (updates configure, etc too)
# don't fail if autogen doesn't exist though: its not in the .tar-ed releases
dist-hook: $(top_srcdir)/autogen
	[[ -x $(top_srcdir)/autogen ]] && $(top_srcdir)/autogen $(distdir) $(top_srcdir) || true


# apply astyle to source files
style:
	astyle --suffix=none --options=./astylerc {src,tests}/*.{c,h}

######
## vodoo incantations to get the test suite automation to build (autotest)

# The `:;' works around a Bash 3.2 bug when the output is not writeable.
$(srcdir)/package.m4: $(top_srcdir)/configure.ac
	:;{ \
	  echo '# Signature of the current package.' && \
	  echo 'm4_define([AT_PACKAGE_NAME],' && \
	  echo '  [$(PACKAGE_NAME)])' && \
	  echo 'm4_define([AT_PACKAGE_TARNAME],' && \
	  echo '  [$(PACKAGE_TARNAME)])' && \
	  echo 'm4_define([AT_PACKAGE_VERSION],' && \
	  echo '  [$(PACKAGE_VERSION)])' && \
	  echo 'm4_define([AT_PACKAGE_STRING],' && \
	  echo '  [$(PACKAGE_STRING)])' && \
	  echo 'm4_define([AT_PACKAGE_BUGREPORT],' && \
	  echo '  [$(PACKAGE_BUGREPORT)])'; \
	  echo 'm4_define([AT_PACKAGE_URL],' && \
	  echo '  [$(PACKAGE_URL)])'; \
	} >'$(srcdir)/package.m4'


EXTRA_DIST += testsuite.at $(srcdir)/package.m4 $(TESTSUITE) atlocal.in
TESTSUITE = $(srcdir)/testsuite
 
check-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' $(TESTSUITEFLAGS)
 
installcheck-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' AUTOTEST_PATH='$(bindir)' \
	  $(TESTSUITEFLAGS) -k '!hello,!unit'
# don't run dependancy and unit tests

clean-local-check:
	test ! -f '$(TESTSUITE)' || \
	  $(SHELL) '$(TESTSUITE)' --clean
	-rm -f $(TESTSUITE)
	-rm -f $(srcdir)/package.m4
	-rm -f atconfig

atconfig: $(srcdir)/config.status
	$(srcdir)/config.status
 
AUTOM4TE = $(SHELL) $(abs_srcdir)/build-aux/missing --run autom4te
AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/testsuite.at $(srcdir)/package.m4
	$(AUTOTEST) -I '$(abs_srcdir)' -o $@.tmp $@.at
	mv $@.tmp $@


# valgrind configuration
MPI_PREFIX=/usr
MPIWRAP_DEBUG_ARGS=debug
MPI_PLATFORM=amd64-linux
VALGRIND_BIN=valgrind
VALGRIND_ENV=MPIWRAP_DEBUG=${MPIWRAP_DEBUG_ARGS} LD_PRELOAD=${MPI_PREFIX}/lib/valgrind/libmpiwrap-${MPI_PLATFORM}.so
VALGRIND_ARGS=--trace-children=yes --track-origins=yes --leak-check=full
#--db-attach=yes 
MPI_BIN=mpirun
MPI_ARGS=-n 4 
VALGRIND=${VALGRIND_ENV} ${VALGRIND_BIN} ${VALGRIND_ARGS} ${MPI_BIN} ${MPI_ARGS}
valgrind: tests/mumps-file
	cd test && ${VALGRIND} mumps-file

.PHONY: help
help:
	@echo "Useful make targets: make [target]"
	@echo "  <none> == all"
	@echo "  all              - build software"
	@echo "  check            - run testsuite"
	@echo "  install"
	@echo "  uninstall"
	@echo "  clean            - remove built files"
	@echo "For developers:"
	@echo "  distclean        - clean to distributed state"
	@echo "  maintainer-clean - delete all generated files, requires ./autogen to rebuild"
	@echo "  dist             - create a .tar.gz for distribution"
	@echo "  distcheck        - .tar, then untar in new dir and test"
	@echo "  tags             - generate ctags for editing source files"
	@echo "  style            - run astyle on source files to tidy formatting"
