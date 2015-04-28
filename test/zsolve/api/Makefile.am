check_PROGRAMS = 
check_PROGRAMS += test_zsolve_api

LIBS = $(top_builddir)/src/zsolve/libzsolve.la 

# Allow accessing groebner's header files as "groebner/HEADER.h"
AM_CXXFLAGS = -I$(top_srcdir)/src/

# 
if HAVE_GMP_WITH_CXX
AM_CXXFLAGS += $(GMP_CFLAGS)
endif

test_zsolve_api_SOURCES = test_zsolve_api.cpp

TESTS =
TESTS += test_zsolve_api

$(TESTS): $(top_builddir)/src/zsolve/libzsolve.la
