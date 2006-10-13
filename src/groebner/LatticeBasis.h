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

#ifndef _4ti2__LatticeBasis_
#define _4ti2__LatticeBasis_

#include "VectorArray.h"

namespace _4ti2_ {

// Computes a lattice basis.
void lattice_basis(const VectorArray& matrix, VectorArray& basis);

IntegerType solve(const VectorArray& matrix, const Vector& rhs, Vector& solution);
void solve(const VectorArray& matrix, const VectorArray& rhs, VectorArray& solutions);

inline
void solve(const VectorArray& matrix, const VectorArray& rhs, VectorArray& solutions)
{
    solutions.renumber(rhs.get_number());
    for (int i = 0; i < rhs.get_number(); ++i)
    {
        IntegerType status = solve(matrix, rhs[i], solutions[i]);
        if (status != 1) { solutions[i].mul(0); }
    }
}

} // namespace _4ti2_

#endif
