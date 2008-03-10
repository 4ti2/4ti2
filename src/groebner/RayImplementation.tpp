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

#include "groebner/RayImplementation.h"
#include "groebner/CircuitOptions.h"
#include "groebner/Globals.h"
#include "groebner/HermiteAlgorithm.h"
#include "groebner/VectorStream.h"

#include "groebner/ShortDenseIndexSet.h"
#include "groebner/LongDenseIndexSet.h"

using namespace _4ti2_;

template <class IndexSet>
RayImplementation<IndexSet>::RayImplementation()
{
    switch (CircuitOptions::instance()->next_column)
    {
    case CircuitOptions::MININDEX:
        compare = &minindex;
        break;
    case CircuitOptions::MAXCUTOFF:
        compare = &maxcutoff;
        break;
    case CircuitOptions::MINCUTOFF:
        compare = &mincutoff;
        break;
    case CircuitOptions::MAXINTER:
        compare = &maxinter;
        break;
    default:
        //compare = &maxinter;
        compare = &maxcutoff;
        break;
    }
}

template <class IndexSet>
RayImplementation<IndexSet>::~RayImplementation()
{
}

template <class IndexSet>
int
RayImplementation<IndexSet>::next_column(
                const VectorArray& vs,
                const IndexSet& remaining,
                int& next_positive_count,
                int& next_negative_count,
                int& next_zero_count)
{
    // Sanity Checks.
    assert(vs.get_size() == remaining.get_size());

    int num_cols = vs.get_size();
    Index c = 0;
    while (c < num_cols && !remaining[c]) { ++c; }
    assert(c < num_cols);
    int next_col = c;
    column_count(vs, next_col, next_positive_count, next_negative_count, next_zero_count);
    while (c < num_cols)
    {
        if (remaining[c])
        {
            int positive_count = 0;
            int negative_count = 0;
            int zero_count = 0;
            column_count(vs, c, positive_count, negative_count, zero_count);
            if ((*compare)(next_positive_count, next_negative_count, next_zero_count,
                        positive_count, negative_count, zero_count))
            {
                next_col = c;
                next_positive_count = positive_count;
                next_negative_count = negative_count;
                next_zero_count = zero_count;
            }
        }
        ++c;
    }
    return next_col;
}

template <class IndexSet>
bool
RayImplementation<IndexSet>::minindex(
                int next_positive_count,
                int next_negative_count,
                int next_zero_count,
                int positive_count,
                int negative_count,
                int zero_count)
{
    return false;
}

template <class IndexSet>
bool
RayImplementation<IndexSet>::maxinter(
                int next_positive_count,
                int next_negative_count,
                int next_zero_count,
                int positive_count,
                int negative_count,
                int zero_count)
{
    return (zero_count > next_zero_count);
}

template <class IndexSet>
bool
RayImplementation<IndexSet>::maxcutoff(
                int next_positive_count,
                int next_negative_count,
                int next_zero_count,
                int positive_count,
                int negative_count,
                int zero_count)
{
    return (negative_count > next_negative_count);
}

template <class IndexSet>
bool
RayImplementation<IndexSet>::mincutoff(
                int next_positive_count,
                int next_negative_count,
                int next_zero_count,
                int positive_count,
                int negative_count,
                int zero_count)
{
    return (negative_count < next_negative_count);
}

// Count how many zero, positive and negative entries there are in a column.
template <class IndexSet>
void
RayImplementation<IndexSet>::column_count(
        const VectorArray& vs,
        Index c,
        int& positive_count,
        int& negative_count,
        int& zero_count)
{
    zero_count = 0;
    positive_count = 0;
    negative_count = 0;
    for (Index r = 0; r < vs.get_number(); ++r)
    {
        if (vs[r][c] == 0) { ++zero_count; }
        else if (vs[r][c] > 0) { ++positive_count; }
        else { ++negative_count; }
    }
}

template <class IndexSet>
void
RayImplementation<IndexSet>::sort(
                VectorArray& vs,
                std::vector<IndexSet>& supports,
                int next_col,
                int next_zero_count,
                int next_positive_count, 
                int next_negative_count)
{
        int zero_index = 0;
        for (int i = 0; i < vs.get_number(); ++i)
        {
            if (vs[i][next_col] == 0)
            {
                vs.swap_vectors(i,zero_index);
                IndexSet::swap(supports[i], supports[zero_index]);
                ++zero_index;
            }
        }
        int positive_index = next_zero_count;
        for (int i = positive_index; i < vs.get_number(); ++i)
        {
            if (vs[i][next_col] > 0)
            {
                vs.swap_vectors(i,positive_index);
                IndexSet::swap(supports[i], supports[positive_index]);
                ++positive_index;
            }
        }
}

template <class IndexSet>
void
RayImplementation<IndexSet>::sort(
                VectorArray& vs,
                std::vector<IndexSet>& supports,
                Fathers& fathers,
                std::vector<IndexSet>& zeros,
                int next_col,
                int next_zero_count,
                int next_positive_count, 
                int next_negative_count)
{
    std::vector<int> perm(vs.get_number());
    for (int i = 0; i < vs.get_number(); ++i) { perm[i] = i; }

    int zero_index = 0;
    for (int i = 0; i < vs.get_number(); ++i)
    {
        if (vs[i][next_col] == 0)
        {
            vs.swap_vectors(i,zero_index);
            IndexSet::swap(supports[i], supports[zero_index]);
            IndexSet::swap(zeros[i], zeros[zero_index]);
            int tmp = fathers[i]; fathers[i] = fathers[zero_index]; fathers[zero_index] = tmp;
            tmp = perm[i]; perm[i] = perm[zero_index]; perm[zero_index] = tmp;
            ++zero_index;
        }
    }
    int positive_index = next_zero_count;
    for (int i = positive_index; i < vs.get_number(); ++i)
    {
        if (vs[i][next_col] > 0)
        {
            vs.swap_vectors(i,positive_index);
            IndexSet::swap(supports[i], supports[positive_index]);
            IndexSet::swap(zeros[i], zeros[positive_index]);
            int tmp = fathers[i]; fathers[i] = fathers[positive_index]; fathers[positive_index] = tmp;
            tmp = perm[i]; perm[i] = perm[positive_index]; perm[positive_index] = tmp;
            ++positive_index;
        }
    }
    //*out << "\nPermutation:\n";
    //for (int i = 0; i < vs.get_number(); ++i) { *out << perm[i] << " "; }
    //*out << "\n";
    std::vector<int> back(vs.get_number());
    for (int i = 0; i < vs.get_number(); ++i) { back[perm[i]] = i; }
    //for (int i = 0; i < vs.get_number(); ++i) { *out << back[i] << " "; }
    //*out << "\n";
    for (int i = 0; i < vs.get_number(); ++i)
    { 
        if (fathers[i] >= 0) { fathers[i] = back[fathers[i]]; }
    }
}

template <class IndexSet>
void
RayImplementation<IndexSet>::create_new_vector(
                VectorArray& vs,
                std::vector<IndexSet>& supports,
                int r1, int r2, int next_col,
                int next_positive_count, int next_negative_count,
                Vector& temp, IndexSet& temp_supp)
{
    if (next_positive_count <= next_negative_count)
    {
        Vector::sub(vs[r2],vs[r1][next_col],
                    vs[r1],vs[r2][next_col],temp);
    }
    else
    {
        Vector::sub(vs[r1],vs[r2][next_col],
                    vs[r2],vs[r1][next_col],temp);
    }
    temp.normalise();
    vs.insert(temp);
    IndexSet::set_union(supports[r1],supports[r2],temp_supp);
    supports.push_back(temp_supp);
#if 0
    *out << "\nADDING VECTOR.\n";
    *out << "R1: " << r1 << "\n";
    *out << vs[r1] << "\n";
    *out << "R2: " << r2 << "\n";
    *out << vs[r2] << "\n";
    *out << temp << "\n";
#endif
}

template <class IndexSet>
void
RayImplementation<IndexSet>::create_new_vector(
                VectorArray& vs,
                std::vector<IndexSet>& supports,
                Fathers& fathers,
                int r1, int r2, int next_col,
                int next_positive_count, int next_negative_count,
                Vector& temp, IndexSet& temp_supp)
{
    if (next_positive_count <= next_negative_count)
    {
        Vector::sub(vs[r2],vs[r1][next_col],
                    vs[r1],vs[r2][next_col],temp);
        fathers.push_back(r1);
    }
    else
    {
        Vector::sub(vs[r1],vs[r2][next_col],
                    vs[r2],vs[r1][next_col],temp);
        fathers.push_back(r2);
    }
    temp.normalise();
    vs.insert(temp);
    IndexSet::set_union(supports[r1],supports[r2],temp_supp);
    supports.push_back(temp_supp);
#if 0
    *out << "\nADDING VECTOR.\n";
    *out << "R1[" << r1 << "]:\n";
    *out << vs[r1] << "\n";
    *out << "R2[" << r2 << "]:\n";
    *out << vs[r2] << "\n";
    *out << "R[" << vs.get_number()-1 << "]: " << fathers.back() << "\n";
    *out << temp << "\n";
#endif
}

