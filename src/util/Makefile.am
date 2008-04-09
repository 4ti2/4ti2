#						-*- Makefile -*-
# 4ti2 -- A software package for algebraic, geometric and combinatorial
# problems on linear spaces.
# 
# Copyright (C) 2006 4ti2 team.
# Main author(s): Raymond Hemmecke, Ralf Hemmecke, Peter Malkin.
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

bin_PROGRAMS =					\
	genmodel				\
	gensymm					\
	normaliz_wrapper			\
	output

noinst_LIBRARIES = lib4ti2util.a

# The output.c file contains C++-isms (or C99-isms), 
# so compile with the C++ compiler. --mkoeppe
CC = $(CXX)
CFLAGS = $(CXXFLAGS)

lib4ti2util_a_SOURCES =				\
	genmodel.c				\
	gensymm.c				\
	orbit.c					\
	output.c				\
	print.c					\
	vector.c

noinst_HEADERS = 				\
	genmodel.h				\
	gensymm.h				\
	myheader.h				\
	orbit.h					\
	output.h				\
	print.h					\
	vector.h

genmodel_SOURCES = genmodel_main.c
genmodel_LDADD = lib4ti2util.a

gensymm_SOURCES = gensymm_main.c
gensymm_LDADD = lib4ti2util.a

output_SOURCES = output_main.c
output_LDADD = lib4ti2util.a

normaliz_wrapper_SOURCES = normaliz_wrapper.c
normaliz_wrapper_LDADD = lib4ti2util.a

