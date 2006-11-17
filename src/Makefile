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
default:  util_target zsolve_target ppi_target groebner_target

# The following target compiles everything in the util directory.
util_target:
	cd util && $(MAKE)

# The following target compiles everything in the util directory.
zsolve_target:
	cd zsolve && $(MAKE)

# The following target compiles everything in the groebner directory.
groebner_target: groebner_int32 groebner_int64 groebner_gmp
# The following target compiles 4ti2 groebner with 16 bit precision.
groebner_int16:
	cd groebner && $(MAKE) INT=int16
# The following target compiles 4ti2 groebner with 32 bit precision.
groebner_int32:
	cd groebner && $(MAKE) INT=int32
# The following target compiles 4ti2 groebner with 64 bit precision.
groebner_int64:
	cd groebner && $(MAKE) INT=int64
# The following target compiles 4ti2 groebner with arbitrary precision.
groebner_gmp:
	cd groebner && $(MAKE) INT=gmp

# The following target compiles everything in the ppi directory.
ppi_target:
	cd ppi && $(MAKE)

clean:
	cd util && $(MAKE) clean
	cd groebner && $(MAKE) clean
	cd zsolve && $(MAKE) clean
