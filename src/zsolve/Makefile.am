#						-*- Makefile -*-
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

AM_CFLAGS = -DNDEBUG

lib_LTLIBRARIES = libzsolve.la

libzsolve_la_SOURCES =				\
	defs.c 					\
	vector.c 				\
	vectorarray.c 				\
	linearsystem.c 				\
	indexarray.c				\
	lattice.c 				\
	matrix.c 				\
	varproperties.c 			\
	valuetrees.c 				\
	cputime.c 				\
	libzsolve.c

include_HEADERS =				\
	cputime.h				\
	defs.h					\
	indexarray.h				\
	lattice.h				\
	libzsolve.h				\
	linearsystem.h				\
	matrix.h				\
	opts.h					\
	valuetrees.h				\
	varproperties.h				\
	vectorarray.h				\
	vector.h

bin_PROGRAMS = zsolve

zsolve_SOURCES = 				\
	solve.c					\
	opts.c

zsolve_LDADD = libzsolve.la

dist_bin_SCRIPTS = hilbert graver
