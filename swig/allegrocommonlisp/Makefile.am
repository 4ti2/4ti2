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

# 4ti2 Makefile

lib_LTLIBRARIES=

# Allegro Common Lisp interface
commonlispdir=$(datadir)/common-lisp/source/4ti2
commonlisp_DATA=4ti2.asd 4ti2-package.lisp common-wrapper-macros.lisp 4ti2-base.lisp 4ti2.cl

lib_LTLIBRARIES += liballegrocl-4ti2.la
liballegrocl_4ti2_la_SOURCES = allegrocl-4ti2.c
#liballegrocl_4ti2_la_LIBADD =
liballegrocl_4ti2_la_LDFLAGS = -L$(top_builddir)/../src/groebner -R$(libdir) -l4ti2int64

AM_CFLAGS = -I$(top_srcdir)/../src 

allegrocl-4ti2.c 4ti2.cl: allegrocl-4ti2.i $(top_srcdir)/../src/4ti2/4ti2.h
	@SWIG@ $(SWIGDEFS) -I$(top_srcdir)/../src \
	-allegrocl -outdir . -cwrap -o allegrocl-4ti2.c \
	-identifier-converter lispify \
	$<
