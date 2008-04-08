# 4ti2 -- A software package for algebraic, geometric and combinatorial
# problems on linear spaces.
# 
# Copyright (C) 2006 4ti2 team.
# Main author(s): Matthias Walter.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 

bin_PROGRAMS = zsolve 

lib_LTLIBRARIES = libzsolve.la

# Allow accessing zsolve's header files as "zsolve/HEADER.h"
AM_CXXFLAGS = -I$(srcdir)/..

# 4ti2 needs the next option to access maximum values of integer types.
AM_CXXFLAGS += -D__STDC_LIMIT_MACROS
# The optimisation flags.
#AM_CXXFLAGS += -DNDEBUG
# Uncomment the next line for debugging.
#AM_CXXFLAGS += -g
# Uncomment the next line for profiling.
#AM_CXXFLAGS += -pg

libzsolvesources = \
	BitSet.cpp \
	Timer.cpp \
	Options.cpp \
	Library.cpp

zsolveincludedir = $(includedir)/zsolve

zsolveinclude_HEADERS = \
	Integer.h \
	BitSet.h \
	Vector.hpp \
	VectorArray.hpp \
	Variables.hpp \
	Relation.hpp \
	LinearSystem.hpp \
	Lattice.hpp \
	Algorithm.hpp \
	Controller.hpp \
	DefaultController.hpp \
	Norms.hpp \
	Options.h \
	Timer.h \
	Exception.h \
	Heuristics.hpp \
	zsolve.h 

zsolvesources =	\
	main.cpp

#noinst_HEADERS = 

# 32 bit precision flags.
libzsolve_la_SOURCES = $(libzsolvesources)
zsolve_SOURCES = $(zsolvesources)
zsolve_LDADD = libzsolve.la
libzsolve_la_CPPFLAGS =
zsolve_CPPFLAGS = 
if HAVE_GMP_WITH_CXX
libzsolve_la_CPPFLAGS += -D_4ti2_GMP_ ${GMP_CFLAGS}
zsolve_CPPFLAGS += -D_4ti2_GMP_ ${GMP_CFLAGS}
zsolve_LDADD += ${GMP_LIBS}
endif

bin_SCRIPTS = hilbert graver
DISTCLEANFILES = hilbert graver

EXTRA_DIST = hilbert.template graver.template

hilbert: hilbert.template
	$(install_sh) -c -m 755 $< $@

graver: graver.template
	$(install_sh) -c -m 755 $< $@
