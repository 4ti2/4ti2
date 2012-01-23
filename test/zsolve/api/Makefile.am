check_PROGRAMS = 
if HAVE_INT32_INT64
check_PROGRAMS += test_zsolve_api
endif

LIBS = $(top_builddir)/src/zsolve/libzsolve.la 

# Allow accessing groebner's header files as "groebner/HEADER.h"
AM_CXXFLAGS = -I$(top_srcdir)/src/

test_zsolve_api_SOURCES = test_zsolve_api.cpp

TESTS =
if HAVE_INT32_INT64
TESTS += test_zsolve_api
endif

$(TESTS): $(top_builddir)/src/zsolve/libzsolve.la
