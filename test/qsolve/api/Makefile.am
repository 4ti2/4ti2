check_PROGRAMS = qsolve_main

LIBS = $(top_builddir)/src/groebner/lib4ti2int64.la 

# Allow accessing groebner's header files as "groebner/HEADER.h"
AM_CXXFLAGS = -I$(top_srcdir)/src/

qsolve_main_SOURCES = qsolve_main.cpp

TESTS = qsolve_main

$(TESTS): $(top_builddir)/src/groebner/lib4ti2int64.la 
