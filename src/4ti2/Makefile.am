# 4ti2 -- A software package for algebraic, geometric and combinatorial
# problems on linear spaces.
# 
# Copyright (C) 2006 4ti2 team.
# Main author(s): Peter Malkin.
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

fortytwoincludedir = $(includedir)/4ti2

fortytwoinclude_HEADERS = 4ti2.h 4ti2xx.h 4ti2_config.h

BUILT_SOURCES = 4ti2_config.h

lib_LTLIBRARIES = lib4ti2common.la

# Allow accessing zsolve's header files as "zsolve/HEADER.h"
AM_CXXFLAGS = -I$(srcdir)/..

lib4ti2common_la_SOURCES = 4ti2API.cpp

if HAVE_GMP_WITH_CXX
AM_CPPFLAGS = -D_4ti2_GMP_ ${GMP_CFLAGS}
AM_LDFLAGS = ${GMP_LIBS}
endif
