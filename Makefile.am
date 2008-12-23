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

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src m4 test
if ENABLE_SWIG
SUBDIRS += swig
endif
## We cannot expect "make dist" to create a distribution including the swig
## directory if that directory (containing an independent tree) is not configured!
DIST_SUBDIRS = $(SUBDIRS)
if !ENABLE_SWIG
dist-hook:
	@echo "**************************************************************"
	@echo "WARNING: Distribution will not contain the swig/ subdirectory"
	@echo "since this build of 4ti2 was configured without --enable-swig."
	@echo "**************************************************************"
endif

EXTRA_DIST = autogen.sh TODO NEWS README README.cvs
