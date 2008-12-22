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

#include "groebner/CircuitImplementation.h"
#include "groebner/Debug.h"
#include "groebner/Globals.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/HermiteAlgorithm.h"

#include <iostream>

using namespace _4ti2_;

template <class IndexSet>
void
CircuitImplementation<IndexSet>::compute(
                const VectorArray& matrix,
                VectorArray& vs,
                const IndexSet& rs,
                const IndexSet& cirs)
{
    // TODO: Complete this implementation.
#if 0
    VectorArray full_matrix(matrix.get_number(), matrix.get_size()+cirs.count());
    VectorArray cirs_matrix(matrix.get_number(), cirs.count());
    VectorArray::project(matrix, cirs, cirs_matrix);
    cirs_matrix.mul(-1);
    VectorArray::concat(matrix, cirs_matrix, full_matrix);

    DEBUG_4ti2(*out << "MATRIX:\n" << matrix;)
    DEBUG_4ti2(*out << "FULL MATRIX:\n" << full_matrix;)

    VectorArray full_vs(vs.get_number(), vs.get_size()+cirs.count());
    VectorArray zero_vs(vs.get_number(), cirs.count(), 0);
    VectorArray::concat(vs, zero_vs, full_vs);
    int index = 0;
    for (int i = 0; i < vs.get_size(); ++i)
    {
        if (cirs[i])
        {
            Vector v(full_vs.get_size(), 0);
            v[i] = 1;
            v[index+vs.get_size()] = 1;
            ++index;
            full_vs.insert(v);
        }
    }
    DEBUG_4ti2(*out << "VS:\n" << vs;)
    DEBUG_4ti2(*out << "FULL VS:\n" << full_vs;)

    IndexSet full_rs(full_vs.get_size(), true);
    for (int i = 0; i < vs.get_size(); ++i) 
    {
        if (!rs[i] && !cirs[i]) { full_rs.unset(i); }
    }
    DEBUG_4ti2(*out << "RS:\n" << rs << "\n";)
    DEBUG_4ti2(*out << "CIRS:\n" << cirs << "\n";)
    DEBUG_4ti2(*out << "FULL RS:\n" << full_rs << "\n";)

    RayAlgorithm algorithm;
    //algorithm.compute(full_matrix, full_vs, full_rs);

    DEBUG_4ti2(*out << "FULL CIRCUITS BEFORE:\n" << full_vs;)
    for (int i = 0; i < full_vs.get_number(); ++i)
    {
        int index = 0;
        for (int j = 0; j < vs.get_size(); ++j)
        {
            if (cirs[j])
            {
                full_vs[i][j] -= full_vs[i][index+vs.get_size()];
                ++index;
            }
        }
    }
    DEBUG_4ti2(*out << "FULL CIRCUITS AFTER:\n" << full_vs;)
    vs.renumber(full_vs.get_number());
    VectorArray::project(full_vs, 0, vs.get_size(), vs);
    Vector zero(vs.get_size(), 0);
    // Remove zero vectors and duplicates.
    // A vector is a duplicate if it is zero on the ray components and the first
    // non-zero entry is negative.
    for (int i = vs.get_number()-1; i >= 0; --i)
    {
        if (vs[i].is_zero(rs) && vs[i] <= zero) { vs.remove(i); }
    }
#endif
}

template <class IndexSet>
int
CircuitImplementation<IndexSet>::next_column(
                const VectorArray& vs,
                const IndexSet& remaining)
{
    // Sanity Check
    assert(vs.get_size() == remaining.get_size());
    int next_pos_count;
    int next_neg_count;
    int next_zero_count;

    int num_cols = vs.get_size();
    Index c = 0;
    while (c < num_cols && !remaining[c]) { ++c; }
    assert(c < num_cols);
    int next_col = c;
    column_count(vs, next_col, next_pos_count, next_neg_count, next_zero_count);
    while (c < num_cols)
    {
        if (remaining[c])
        {
            int pos_count = 0;
            int neg_count = 0;
            int zero_count = 0;
            column_count(vs, c, pos_count, neg_count, zero_count);
            if (zero_count > next_zero_count)
            {
                next_col = c;
                next_pos_count = pos_count;
                next_neg_count = neg_count;
                next_zero_count = zero_count;
            }
        }
        ++c;
    }
    return next_col;
}

// Count the number of rows with positive, negative, and zero entries in column
// c.
template <class IndexSet>
void
CircuitImplementation<IndexSet>::column_count(
        const VectorArray& vs,
        Index c,
        int& pos_count,
        int& neg_count,
        int& zero_count)
{
    zero_count = 0;
    pos_count = 0;
    neg_count = 0;
    for (Index r = 0; r < vs.get_number(); ++r)
    {
        if (vs[r][c] == 0) { ++zero_count; }
        else if (vs[r][c] > 0) { ++pos_count; }
        else { ++neg_count; }
    }
}

template <class IndexSet>
void
CircuitImplementation<IndexSet>::switch_supports(
                int start, int end,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps)
{
    for (int i = start; i < end; ++i) { IndexSet::swap(pos_supps[i], neg_supps[i]); }
}


// Pushes zeros to the beginning.
template <class IndexSet>
void
CircuitImplementation<IndexSet>::sort_nonzeros(
                VectorArray& vs,
                int start, int end,
                std::vector<bool>& rays,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps,
                int next_col,
                int& middle)
{
    assert(start >= 0 && start <= end && end <= vs.get_number());
    int index = start;
    for (int i = start; i < end; ++i)
    {
        if (vs[i][next_col] != 0)
        {
            vs.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            IndexSet::swap(pos_supps[i], pos_supps[index]);
            IndexSet::swap(neg_supps[i], neg_supps[index]);
            bool tmp = rays[i];
            rays[i] = rays[index];
            rays[index] = tmp;
            ++index;
        }
    }
    middle = index;
}

// Pushes positives to the beginning.
template <class IndexSet>
void
CircuitImplementation<IndexSet>::sort_positives(
                VectorArray& vs,
                int start, int end,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps,
                int next_col,
                int& middle)
{
    assert(start >= 0 && start <= end && end <= vs.get_number());
    int index = start;
    for (int i = start; i < end; ++i)
    {
        if (vs[i][next_col] > 0)
        {
            vs.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            IndexSet::swap(pos_supps[i], pos_supps[index]);
            IndexSet::swap(neg_supps[i], neg_supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes negatives to the beginning.
template <class IndexSet>
void
CircuitImplementation<IndexSet>::sort_negatives(
                VectorArray& vs,
                int start, int end,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps,
                int next_col,
                int& middle)
{
    assert(start >= 0 && start <= end && end <= vs.get_number());
    int index = start;
    for (int i = start; i < end; ++i)
    {
        if (vs[i][next_col] < 0)
        {
            vs.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            IndexSet::swap(pos_supps[i], pos_supps[index]);
            IndexSet::swap(neg_supps[i], neg_supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes rays to the beginning.
template <class IndexSet>
void
CircuitImplementation<IndexSet>::sort_rays(
                VectorArray& vs,
                int start, int end,
                std::vector<bool>& rays,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps,
                int& middle)
{
    int index = start;
    for (int i = start; i < end; ++i)
    {
        if (rays[i])
        {
            vs.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            IndexSet::swap(pos_supps[i], pos_supps[index]);
            IndexSet::swap(neg_supps[i], neg_supps[index]);
            rays[i] = rays[index];
            rays[index] = true;
            ++index;
        }
    }
    middle = index;
}

// Splits vs into rays and circuits.
template <class IndexSet>
void
CircuitImplementation<IndexSet>::split_rays(
                VectorArray& vs,
                std::vector<bool>& rays,
                VectorArray& circuits)
{
    int index = 0;
    for (int i = 0; i < vs.get_number(); ++i)
    {
        if (rays[i])
        {
            vs.swap_vectors(i,index);
            ++index;
        }
    }
    VectorArray::transfer(vs, index, vs.get_number(), circuits, 0);
}


template <class IndexSet>
void
CircuitImplementation<IndexSet>::update_supports(
                std::vector<IndexSet>& supps,
                int index, int start, int end)
{
    for (int i = start; i < end; ++i) { supps[i].set(index); }
}

template <class IndexSet>
void
CircuitImplementation<IndexSet>::check(
                const Vector& v,
                const IndexSet& supp,
                const IndexSet& pos_supp,
                const IndexSet& neg_supp,
                const IndexSet& remaining)
{
    for (int i = 0; i < v.get_size(); ++i)
    {
        if (!remaining[i])
        {
            if ((v[i] != 0) != supp[i]) { *out << "Support Check failed.\n"; }
            if ((v[i] > 0) != pos_supp[i]) { *out << "Positive Support Check failed.\n"; }
            if ((v[i] < 0) != neg_supp[i]) { *out << "Negative Support Check failed.\n"; }
        }
    }
}
