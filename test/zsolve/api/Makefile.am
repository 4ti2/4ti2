check_PROGRAMS = 
check_PROGRAMS += test_zsolve_api

test_zsolve_api_LDADD = $(top_builddir)/src/zsolve/libzsolve.la $(top_builddir)/src/4ti2/lib4ti2common.la -L$(top_builddir)/src/4ti2/.libs

# Allow accessing groebner's header files as "groebner/HEADER.h"
AM_CXXFLAGS = -I$(top_srcdir)/src/

# 
if HAVE_GMP
AM_CXXFLAGS += $(GMP_CFLAGS)
endif

test_zsolve_api_SOURCES = test_zsolve_api.cpp

TESTS =
TESTS += test_zsolve_api

$(TESTS): $(top_builddir)/src/zsolve/libzsolve.la
$(TESTS): $(top_builddir)/src/4ti2/lib4ti2common.la
