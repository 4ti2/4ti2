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

#ifndef _4ti2_qsolve__DiagonalAlgorithm_
#define _4ti2_qsolve__DiagonalAlgorithm_

#include "qsolve/HermiteAlgorithm.h"
#include "qsolve/IndexSet.h"

namespace _4ti2_
{

template <class Matrix>
Index diagonal(Matrix& vs);

template <class Matrix, class ColumnSet>
Index diagonal(Matrix& vs, ColumnSet& pivots);

template <class Matrix, class RowIter, class ColIter>
Index diagonal(Matrix& ps, RowIter rmin, RowIter rmax, ColIter cmin, ColIter cmax);

template <class Matrix, class RowIter, class ColIter, class ColumnSet>
Index diagonal(Matrix& ps, RowIter rmin, RowIter rmax, ColIter cmin, ColIter cmax, ColumnSet& pivots);

// Template function definitions.

template <class Matrix>
inline
Index
diagonal(Matrix& vs)
{
    return diagonal(vs, 0, vs.get_number(), 0, vs.get_size());
}

template <class Matrix, class ColumnSet>
Index
diagonal(Matrix& vs, ColumnSet& pivots)
{
    return diagonal(vs, 0, vs.get_number(), 0, vs.get_size(), pivots);
}

template <class Matrix, class RowIter, class ColIter>
Index
diagonal(Matrix& vs, RowIter rmin, RowIter rmax, ColIter cmin, ColIter cmax)
{
    ColIter pivot_col = cmin;
    RowIter pivot_row = rmin;
    while (pivot_col != cmax && pivot_row != rmax) {
        if (pivot(vs, pivot_row, rmax, pivot_col)) {
            // TODO: Should we start from 0?
            for (Index r = 0; r != pivot_row; ++r) {
                if (vs[r][pivot_col] != 0) {
                    typename Matrix::DataType g0,p0,q0,p1,q1;
                    euclidean(vs[r][pivot_col],vs[pivot_row][pivot_col],
                                        g0,p0,q0,p1,q1);
                    vs[r].add(vs[r], p1, vs[pivot_row], q1);
                }
            }
            ++pivot_row;
        }
        ++pivot_col;
    }
    vs.normalise();
    return pivot_row;
}

template <class Matrix, class RowIter, class ColIter, class ColumnSet>
Index
diagonal(Matrix& vs, RowIter rmin, RowIter rmax, ColIter cmin, ColIter cmax, ColumnSet& pivots)
{
    ColIter pivot_col = cmin;
    RowIter pivot_row = rmin;
    while (pivot_col != cmax && pivot_row != rmax) {
        if (pivot(vs, pivot_row, rmax, pivot_col)) {
            // TODO: Should we start from 0?
            for (Index r = 0; r != pivot_row; ++r) {
                if (vs[r][pivot_col] != 0) {
                    typename Matrix::DataType g0,p0,q0,p1,q1;
                    euclidean(vs[r][pivot_col],vs[pivot_row][pivot_col],
                                        g0,p0,q0,p1,q1);
                    vs[r].add(vs[r], p1, vs[pivot_row], q1);
                }
            }
            ++pivot_row;
            pivots.set(pivot_col);
        }
        ++pivot_col;
    }
    vs.normalise();
    return pivot_row;
}

} // namespace _4ti2_

#endif
