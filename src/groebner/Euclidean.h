/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 
*/

#ifndef _4ti2__Euclidean_
#define _4ti2__Euclidean_

#include "DataType.h"

namespace _4ti2_ {

// g0 is the gcd of a and b (g0 > 0).
// p0*a + q0*b = g0.
// p1*a + q1*b = 0 (p1 > 0).

void
euclidean(      IntegerType a,
                IntegerType b,
                IntegerType& g0);

void
euclidean(      IntegerType a,
                IntegerType b,
                IntegerType& g0,
                IntegerType& p0,
                IntegerType& q0);

void
euclidean(      IntegerType a,
                IntegerType b,
                IntegerType& g0,
                IntegerType& p0,
                IntegerType& q0,
                IntegerType& p1,
                IntegerType& q1);

void
lcm(            IntegerType a,
                IntegerType b,
                IntegerType& lcm);

} // namespace _4ti2_

#endif
