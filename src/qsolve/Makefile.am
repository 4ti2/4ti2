# 4ti2 -- A software package for algebraic, geometric and combinatorial
# problems on linear spaces.
# 
# Copyright (C) 2006 4ti2 team.
# Main author(s): Peter Malkin
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

bin_PROGRAMS = 4ti2 

lib_LTLIBRARIES = lib4ti2.la

# Allow accessing 4ti2's header files as "4ti2/HEADER.h"
AM_CXXFLAGS = -I$(srcdir)/..

# 4ti2 needs the next option to access maximum values of integer types.
AM_CXXFLAGS += -D__STDC_LIMIT_MACROS
# The optimisation flags.
AM_CXXFLAGS += -DNDEBUG
# The threaded flags.
AM_CXXFLAGS += -D_PTHREADS -D_REENTRANT
# Uncomment the next line for debugging.
#AM_CXXFLAGS += -g
# Uncomment the next line for profiling.
#AM_CXXFLAGS += -pg
# Uncomment the next lines for multi-threaded running.
#AM_CXXFLAGS += -D_REENTRANT
#4ti2_LDADD += -lpthread

lib4ti2sources = \
	Globals.cpp \
	IndexSetD.cpp \
	IndexSetDS.cpp \
	IndexSetStream.cpp \
	Options.cpp \
	Timer.cpp

4ti2includedir = $(includedir)/4ti2

4ti2include_HEADERS = \
	QSolveAlgorithm.h \
	QSolveAlgorithm.hpp \
	QSolveAPI.h \
	QSolveAPI.hpp \
	QSolveConsOrder.h \
	QSolveVariant.h \
	TypeConversion.h \
	Cone.h \
	DataType.h \
	Euclidean.h \
	Globals.h \
	HermiteAlgorithm.h \
	Index.h \
	IndexSet.h \
	IndexSetD.h \
	IndexSetDS.h \
	IndexSetStream.h \
	Options.h \
	Size.h \
	Timer.h \
	Vector.h \
	VectorD.h \
	VectorS.h \
	VectorStream.h \
	VectorArray.h \
	VectorArrayStream.h

4ti2sources =	\
	main.cpp

#noinst_HEADERS = 

# 32 bit precision flags.
lib4ti2_la_SOURCES = $(lib4ti2sources)
4ti2_SOURCES = $(4ti2sources)
4ti2_LDADD = lib4ti2.la
lib4ti2_la_CPPFLAGS =
4ti2_CPPFLAGS = 
if HAVE_GMP_WITH_CXX
lib4ti2_la_CPPFLAGS += -D_4ti2_GMP_ ${GMP_CFLAGS}
4ti2_CPPFLAGS += -D_4ti2_GMP_ ${GMP_CFLAGS}
4ti2_LDADD += ${GMP_LIBS}
endif
4ti2_LDADD += -lpthread
4ti2_LDADD += -lm
