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

# By default compile 4ti2 with 32 bit precision, 64 bit precision, and arbitrary precision.
default:
	cd src && $(MAKE)

# The check target runs automatic checks to verify that 4ti2 was compiled
# correctly and is running correctly.
check:
	cd test && ./checkall

# Remove automatically generated files.
clean:
	cd src && $(MAKE) clean
	rm -rf obj bin lib
