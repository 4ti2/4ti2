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

#ifndef _4ti2_qsolve__HermiteAlgorithm_
#define _4ti2_qsolve__HermiteAlgorithm_

#include "qsolve/Index.h"

namespace _4ti2_ {

template <class Matrix>
Index upper_triangle(Matrix& ps);

template <class Matrix, class RowIter, class ColIter>
Index upper_triangle(Matrix& mat, RowIter rmin, RowIter rmax, ColIter cmin, ColIter cmax);

template <class Matrix, class RowIter, class ColIter, class ColumnSet>
Index upper_triangle(Matrix& mat, RowIter rmin, RowIter rmax, ColIter cmin, ColIter cmax, ColumnSet& pivots);

template <class Matrix, class RowIter>
bool pivot(Matrix& mat, RowIter rmin, RowIter rmax, Index pivot_col);

template <class Matrix>
inline
Index
upper_triangle(Matrix& mat)
{
    return upper_triangle(mat, 0, mat.get_number(), 0, mat.get_size());
}

template <class Matrix, class RowIter>
inline
bool
pivot(Matrix& mat, RowIter rmin, RowIter rmax, Index pivot_col) 
{
    Index pivot_row = rmin;
    int index = -1;
    for (Index r = rmin; r != rmax; ++r) {
        if (mat[r][pivot_col] < 0) { mat[r].muleq(-1); }
        if (index == -1 && mat[r][pivot_col] != 0) { index = r; }
    }
    if (index == -1) { return false; }
    mat.swap_rows(pivot_row, index);
    while (true) {
        assert(mat[pivot_row][pivot_col] > 0);
        Index min_row = pivot_row;
        bool all_zeros = true;
        for (Index r = rmin+1; r != rmax; ++r) {
            if (mat[r][pivot_col] > 0) {
                if (mat[r][pivot_col] < mat[min_row][pivot_col]) { min_row = r; }
                all_zeros = false;
            }
        }
        if (all_zeros) { return true; }
        mat.swap_rows(pivot_row, min_row);
        for (Index r = rmin+1; r != rmax; ++r) {
            if (mat[r][pivot_col] != 0) {
                mat[r].addeq(mat[pivot_row], 
                            -mat[r][pivot_col]/mat[pivot_row][pivot_col]);
            }
        }
    }
}

template <class Matrix, class RowIter, class ColIter>
Index
upper_triangle(Matrix& mat, RowIter rmin, RowIter rmax, ColIter cmin, ColIter cmax)
{
    RowIter pivot_row = rmin;
    ColIter pivot_col = cmin;
    while (pivot_col != cmax && pivot_row != rmax) {
        if (pivot(mat, pivot_row, rmax, pivot_col)) { ++pivot_row; }
        ++pivot_col;
    }
    return pivot_row;
}

template <class Matrix, class RowIter, class ColIter, class ColumnSet>
Index
upper_triangle(Matrix& mat, RowIter rmin, RowIter rmax, ColIter cmin, ColIter cmax, ColumnSet& pivots)
{
    RowIter pivot_row = rmin;
    ColIter pivot_col = cmin;
    while (pivot_col != cmax && pivot_row != rmax) {
        if (pivot(mat, pivot_row, rmax, pivot_col)) { ++pivot_row; pivots.set(pivot_col); }
        ++pivot_col;
    }
    return pivot_row;
}

} // namespace _4ti2_

#endif
