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

#ifndef _4ti2_groebner__DiagonalAlgorithm_
#define _4ti2_groebner__DiagonalAlgorithm_

#include "groebner/VectorArray.h"
#include "groebner/BitSet.h"

namespace _4ti2_
{

Index diagonal(VectorArray& vs);

Index diagonal(VectorArray& vs, Index num_cols);

template <class ColumnSet>
Index diagonal(VectorArray& vs, const ColumnSet& cols);

template <class ColumnSet>
Index diagonal(VectorArray& vs, const ColumnSet& cols, int row);

inline
Index
diagonal(VectorArray& vs)
{
    return _4ti2_::diagonal(vs, vs.get_size());
}

} // namespace _4ti2_

// Include template functions.
#include "groebner/DiagonalAlgorithm.tpp"

#endif
