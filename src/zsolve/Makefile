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
CC = gcc
LD = ld
BINDIR = ../../bin
BINARY_NAME = ${BINDIR}/zsolve
LIBDIR = ../../lib
LIBARY_NAME = ${LIBDIR}/libzsolve.so

SCRIPTS = ${BINDIR}/hilbert ${BINDIR}/graver

OBJDIR = ../../obj/zsolve
OBJECTS = ${OBJDIR}/defs.o ${OBJDIR}/vector.o ${OBJDIR}/vectorarray.o ${OBJDIR}/linearsystem.o ${OBJDIR}/indexarray.o ${OBJDIR}/lattice.o ${OBJDIR}/matrix.o ${OBJDIR}/varproperties.o ${OBJDIR}/valuetrees.o ${OBJDIR}/cputime.o ${OBJDIR}/libzsolve.o
SOLVE_OBJECTS = ${OBJDIR}/solve.o ${OBJDIR}/opts.o

# release?
#DEFINES += -DNDEBUG

CFLAGS = ${DEFINES}
ifeq (${HOSTTYPE},sun4)
CC=cc
CFLAGS += -fast -xarch=generic
LDFLAGS = -G
else
CFLAGS += -Wall -O3 -fomit-frame-pointer -fpic
LDFLAGS = -shared
endif

all: ${LIBARY_NAME} ${BINARY_NAME} ${SCRIPTS}

${OBJDIR}/%.o : %.c
	@mkdir -p ${OBJDIR}
	${CC} ${CFLAGS} -c $< -o $@

${LIBARY_NAME}: ${OBJECTS}
	@mkdir -p ${LIBDIR}
	${CC} ${LDFLAGS} -o ${LIBARY_NAME} ${OBJECTS}

${BINARY_NAME}: ${OBJECTS} ${SOLVE_OBJECTS}
	@mkdir -p ${BINDIR}
	${CC} ${CFLAGS} -o ${BINARY_NAME} ${OBJECTS} ${SOLVE_OBJECTS}

${SCRIPTS}: ${BINDIR}/% : %
	@mkdir -p ${BINDIR}
	cp -p $^ $@
	@touch $@

clean:
	rm -rf ${OBJDIR} ${LIBDIR} ${BINARY_NAME}
