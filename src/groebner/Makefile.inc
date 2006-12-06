# 4ti2 -- A software package for algebraic, geometric and combinatorial
# problems on linear spaces.
# 
# Copyright (C) 2006 4ti2 team.
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

# 4ti2 always uses g++ as the c++ compiler.
COMPILER=g++

# Turn on all warnings.
CXXFLAGS=-Wall
# The libraries.
LIBS=
# The include path for source files.
INCLUDE=
# The include path for libraries.
LIBINCLUDE=

# 4ti2 uses the glpk linear programming library (GLPL).
LIBS+=/usr/lib/libglpk.a
# Add include path to the file glpk.h if necessary.
INCLUDE+=
# Add library include path for the library libglpk if necessary.
LIBINCLUDE+=

# 4ti2 needs the next option to access maximum values of integer types.
CXXFLAGS = -D__STDC_LIMIT_MACROS
# The optimisation flags.
CXXFLAGS += -O3 -DNDEBUG
# Uncomment the next line for debugging.
#CXXFLAGS += -g
# Uncomment the next line for profiling.
#CXXFLAGS += -pg

# The arithmetic precision. The default is 64 bit precision.
# This value should be overridden on the command line to compile with different
# precisions.
INT=int64

# 64 bit precision flags.
ifeq ($(INT),int64)
CXXFLAGS+= -D_4ti2_INT64_
# The flag -ftrapv means that arithmetic overflow checking is enabled.
# This only works with the g++ compiler version 3.4 and newer versions.
TRAPV=-ftrapv
CXXFLAGS+=$(TRAPV)
endif

# 32 bit precision flags.
ifeq ($(INT),int32)
CXXFLAGS+= -D_4ti2_INT32_
endif

# 16 bit precision flags.
ifeq ($(INT),int16)
CXXFLAGS+= -D_4ti2_INT16_
endif

# Arbitrary precision flags.
# 4ti2 uses GMP (GLPL), an arbitrary precision arithmetic library.
ifeq ($(INT),gmp)
LIBS+=/usr/lib/libgmpxx.a /usr/lib/libgmp.a 
CXXFLAGS+= -D_4ti2_GMP_
# Add include path to gmp.h and gmpxx.h files if necessary.
INCLUDE+=
# Add include path for the libraries libgmp and libgmpxx if necessary.
LIBINCLUDE+=
endif
