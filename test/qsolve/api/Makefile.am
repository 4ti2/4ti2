check_PROGRAMS = qsolve_main test_input_conversion_64_to_gmp test_input_conversion_gmp_to_64

# Allow accessing groebner's header files as "groebner/HEADER.h"
AM_CXXFLAGS = -I$(top_srcdir)/src/

qsolve_main_SOURCES = qsolve_main.cpp
qsolve_main_LDADD = $(top_builddir)/src/groebner/lib4ti2int64.la $(top_builddir)/src/4ti2/lib4ti2common.la 

test_input_conversion_64_to_gmp_SOURCES = test_input_conversion_64.cpp
test_input_conversion_64_to_gmp_LDADD = $(top_builddir)/src/groebner/lib4ti2gmp.la $(top_builddir)/src/4ti2/lib4ti2common.la

test_input_conversion_gmp_to_64_SOURCES = test_input_conversion_gmp.cpp
test_input_conversion_gmp_to_64_LDADD = $(top_builddir)/src/groebner/lib4ti2int64.la $(top_builddir)/src/4ti2/lib4ti2common.la

TESTS = qsolve_main test_input_conversion_64_to_gmp.sh test_input_conversion_gmp_to_64.sh

EXTRA_DIST = test_input_conversion_64_to_gmp.sh test_input_conversion_gmp_to_64.sh

$(TESTS): $(top_builddir)/src/groebner/lib4ti2int64.la 
