check_PROGRAMS = 
check_PROGRAMS += qsolve_main rays_main circuits_main \
	test_input_conversion_64_to_gmp test_input_conversion_gmp_to_64

AM_DEFAULT_SOURCE_EXT = .cpp

# Allow accessing groebner's header files as "groebner/HEADER.h"
AM_CXXFLAGS = -I$(top_srcdir)/src/

# 
if HAVE_GMP_WITH_CXX
AM_CXXFLAGS += $(GMP_CFLAGS)
endif

qsolve_main_LDADD = $(top_builddir)/src/groebner/lib4ti2int64.la $(top_builddir)/src/4ti2/lib4ti2common.la 
rays_main_LDADD = $(top_builddir)/src/groebner/lib4ti2int64.la $(top_builddir)/src/4ti2/lib4ti2common.la 
circuits_main_LDADD = $(top_builddir)/src/groebner/lib4ti2int64.la $(top_builddir)/src/4ti2/lib4ti2common.la 

test_input_conversion_64_to_gmp_SOURCES = test_input_conversion_64.cpp
test_input_conversion_64_to_gmp_LDADD = $(top_builddir)/src/groebner/lib4ti2gmp.la $(top_builddir)/src/4ti2/lib4ti2common.la -lgmpxx -lgmp

test_input_conversion_gmp_to_64_SOURCES = test_input_conversion_gmp.cpp
test_input_conversion_gmp_to_64_LDADD = $(top_builddir)/src/groebner/lib4ti2int64.la $(top_builddir)/src/4ti2/lib4ti2common.la -lgmpxx -lgmp

TESTS = 
TESTS += qsolve_main rays_main circuits_main
#test_input_conversion_64_to_gmp.sh test_input_conversion_gmp_to_64.sh

EXTRA_DIST = test_input_conversion_64_to_gmp.sh test_input_conversion_gmp_to_64.sh

$(TESTS): $(top_builddir)/src/groebner/lib4ti2int64.la 
