/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Peter Malkin.

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

#ifndef _4ti2__HermiteAlgorithm_
#define _4ti2__HermiteAlgorithm_

#include "VectorArray.h"
#include "Index.h"

namespace _4ti2_ {

// Computes hermite normal form.
Index hermite(VectorArray& ps);

Index hermite(VectorArray& ps, Index num_cols);

template <class ColumnSet>
Index hermite(VectorArray& ps, const ColumnSet& proj);

template <class ColumnSet>
Index hermite(VectorArray& ps, const ColumnSet& proj, int row);

Index upper_triangle(VectorArray& ps);

Index upper_triangle(VectorArray& ps, Index num_cols);

Index upper_triangle(VectorArray& ps, Index num_rows, Index num_cols);

template <class ColumnSet>
Index upper_triangle(VectorArray& ps, const ColumnSet& proj);

template <class ColumnSet>
Index upper_triangle(VectorArray& ps, const ColumnSet& proj, int row);

bool matrix_rank_check(VectorArray& ps);

inline
Index
hermite(VectorArray& vs)
{
    return _4ti2_::hermite(vs, vs.get_size());
}

inline
Index
upper_triangle(VectorArray& vs, Index num_cols)
{
    return upper_triangle(vs, vs.get_number(), num_cols);
}

inline
Index
upper_triangle(VectorArray& vs)
{
    return upper_triangle(vs, vs.get_number(), vs.get_size());
}

template <class ColumnSet>
inline
Index
hermite(VectorArray& vs, const ColumnSet& cols)
{
    return _4ti2_::hermite(vs, cols, 0);
}

template <class ColumnSet>
inline
Index
upper_triangle(VectorArray& vs, const ColumnSet& cols)
{
    return _4ti2_::upper_triangle(vs, cols, 0);
}

} // namespace _4ti2_

// Include template functions.
#include "HermiteAlgorithm.tpp"

#endif
