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

#include "HermiteAlgorithm.h"
#include "BitSetStream.h"
#include "VectorArrayStream.h"
#include "Globals.h"

//#define DEBUG_4ti2(X) X
#include "Debug.h"

using namespace _4ti2_;

Index
_4ti2_::hermite(VectorArray& vs, Index num_cols)
{
    assert(num_cols <= vs.get_size());
    Index pivot_col = 0;
    Index pivot_row = 0;
    while (pivot_col < num_cols && pivot_row < vs.get_number())
    {
        DEBUG_4ti2(*out << "VS:\n" << vs << "\n";)
        int index = -1;
        for (Index r = pivot_row; r < vs.get_number(); ++r)
        {
            if (vs[r][pivot_col] < 0) { vs[r].mul(-1); }
            if (index == -1 && vs[r][pivot_col] != 0) { index = r; }
        }
        if (index != -1)
        {
            vs.swap_vectors(pivot_row, index);
            while (true)
            {
                assert(vs[pivot_row][pivot_col] > 0);
                int min_row = pivot_row;
                bool all_zeros = true;
                for (Index r = pivot_row+1; r < vs.get_number(); ++r)
                {
                    if (vs[r][pivot_col] > 0)
                    {
                        if (vs[r][pivot_col] < vs[min_row][pivot_col])
                        {
                            min_row = r;
                        }
                        all_zeros = false;
                    }
                }
                if (all_zeros) { break; }
                vs.swap_vectors(pivot_row, min_row);
                for (Index r = pivot_row+1; r < vs.get_number(); ++r)
                {
                    if (vs[r][pivot_col] != 0)
                    {
                       IntegerType mul = vs[r][pivot_col]/vs[pivot_row][pivot_col];
                       vs[r].sub(vs[pivot_row], mul);
                    }
                }
            }
            for (Index r = 0; r < pivot_row; ++r)
            {
                if (vs[r][pivot_col] != 0)
                {
                    IntegerType mul = vs[r][pivot_col]/vs[pivot_row][pivot_col];
                    Vector::sub(vs[r], 1, vs[pivot_row], mul, vs[r]);
                    //if (vs[r][pivot_col] < 0)
                    //{
                    //    Vector::add(vs[r], vs[pivot_row], vs[r]);
                    //}
                    if (vs[r][pivot_col] > 0)
                    {
                        Vector::sub(vs[r], vs[pivot_row], vs[r]);
                    }
                }
            }
            ++pivot_row;
        }
        ++pivot_col;
    }
    return pivot_row;
}

Index
_4ti2_::upper_triangle(VectorArray& vs, int num_rows, int num_cols)
{
    Index pivot_col = 0;
    Index pivot_row = 0;
    while (pivot_col < num_cols && pivot_row < num_rows)
    {
        int index = -1;
        for (Index r = pivot_row; r < num_rows; ++r)
        {
            if (vs[r][pivot_col] < 0) { vs[r].mul(-1); }
            if (index == -1 && vs[r][pivot_col] != 0) { index = r; }
        }
        if (index != -1)
        {
            vs.swap_vectors(pivot_row, index);
            while (true)
            {
                assert(vs[pivot_row][pivot_col] > 0);
                int min_row = pivot_row;
                bool all_zeros = true;
                for (Index r = pivot_row+1; r < num_rows; ++r)
                {
                    if (vs[r][pivot_col] > 0)
                    {
                        if (vs[r][pivot_col] < vs[min_row][pivot_col])
                        {
                            min_row = r;
                        }
                        all_zeros = false;
                    }
                }
                if (all_zeros) { break; }
                vs.swap_vectors(pivot_row, min_row);
                for (Index r = pivot_row+1; r < num_rows; ++r)
                {
                    if (vs[r][pivot_col] != 0)
                    {
                        IntegerType mul = vs[r][pivot_col]/vs[pivot_row][pivot_col];
                        vs[r].sub(vs[pivot_row], mul);
                    }
                }
            }
            ++pivot_row;
            //*out << "Pivot Row = " << pivot_row << std::endl;
        }
        ++pivot_col;
    }
    return pivot_row;
}

#if 0
bool
_4ti2_::matrix_rank_check(VectorArray& vs)
{
    bool check = false;
    Index pivot_col = 0;
    Index pivot_row = 0;
    while (pivot_col < vs.get_size() && pivot_row < vs.get_number())
    {
        Index i = pivot_row;
        while (i < vs.get_number() && vs[i][pivot_col] == 0) { ++i; }
        if (i < vs.get_number())
        {
            vs.swap_vectors(pivot_row, i);
            for (Index r = pivot_row+1; r < vs.get_number(); ++r)
            {
                if (vs[r][pivot_col] != 0)
                {
                    Vector::sub(vs[r],vs[pivot_row][pivot_col],
                                vs[pivot_row], vs[r][pivot_col], vs[r]);
                }
            }
            ++pivot_row;
        }
        else
        {
            if (check) { return false; }
            else { check = true; }
        }
        ++pivot_col;
    }
    return true;
}
#endif
