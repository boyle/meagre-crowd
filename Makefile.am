# build to GNU standards
AUTOMAKE_OPTIONS = gnu


# the actual program to be installed at the end
bin_PROGRAMS = meagre-crowd
meagre_crowd_SOURCES = src/meagre-crowd.c

# the tests
TEST_SRC_DIR=test
check_PROGRAMS = ${TEST_SRC_DIR}/file ${TEST_SRC_DIR}/helloworld-mpi ${TEST_SRC_DIR}/helloworld-mumps ${TEST_SRC_DIR}/mumps-file ${TEST_SRC_DIR}/unit_perftimer

__TEST_SRC_DIR__unit_perftimer_SOURCES = ${TEST_SRC_DIR}/unit_perftimer.c src/perftimer.c
__TEST_SRC_DIR__unit_perftimer_CFLAGS = ${CFLAGS} -Werror
#__TEST_SRC_DIR__unit_perftimer_CFLAGS = ${CFLAGS} -Wextra -Werror

clean-local: clean-local-check

######
## vodoo incantations to get the test suite automation to build

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


EXTRA_DIST = testsuite.at $(srcdir)/package.m4 $(TESTSUITE) atlocal.in
TESTSUITE = $(srcdir)/testsuite
 
check-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' $(TESTSUITEFLAGS)
 
installcheck-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' AUTOTEST_PATH='$(bindir)' \
	  $(TESTSUITEFLAGS) -k '!hello'
 
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
VALGRIND_ARGS=--trace-children=yes --track-origins=yes
#--db-attach=yes 
MPI_BIN=mpirun
MPI_ARGS=-n 4 
VALGRIND=${VALGRIND_ENV} ${VALGRIND_BIN} ${VALGRIND_ARGS} ${MPI_BIN} ${MPI_ARGS}
valgrind: test/mumps-file
	cd test && ${VALGRIND} mumps-file
