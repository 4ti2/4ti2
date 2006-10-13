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

#include "DiagonalAlgorithm.h"
#include "Euclidean.h"
#include "HermiteAlgorithm.h"
#include "VectorArrayStream.h"

using namespace _4ti2_;

Index
_4ti2_::diagonal(VectorArray& vs, Index num_cols)
{
    assert(num_cols <= vs.get_size());
    hermite(vs, num_cols);
    Index pivot_col = 0;
    Index pivot_row = 0;
    while (pivot_col < num_cols && pivot_row < vs.get_number())
    {
        if (vs[pivot_row][pivot_col] != 0)
        {
            for (Index r = 0; r < pivot_row; ++r)
            {
                if (vs[r][pivot_col] != 0)
                {
                    IntegerType g0,p0,q0,p1,q1;
                    euclidean(vs[r][pivot_col],vs[pivot_row][pivot_col],
                                        g0,p0,q0,p1,q1);
                    Vector::add(vs[r], p1, vs[pivot_row], q1, vs[r]);
                }
            }
            ++pivot_row;
        }
        ++pivot_col;
    }
    vs.normalise();
    return pivot_row;
}
