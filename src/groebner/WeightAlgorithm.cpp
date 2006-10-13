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

#include "WeightAlgorithm.h"

using namespace _4ti2_;

bool
WeightAlgorithm::check_weights(
                const VectorArray& matrix,
                const VectorArray& lattice,
                const BitSet& urs,
                VectorArray& weights)
{
    Vector result(matrix.get_number());
    for (Index i = 0; i < weights.get_number(); ++i)
    {
        for (Index j = 0; j < matrix.get_number(); ++j)
        {
            if (Vector::dot(weights[i], matrix[j]) != 0) { return false; }
        }
    }

    for (Index i = 0; i < weights.get_number(); ++i)
    {
        if (violates_urs(weights[i], urs)) { return false; }
    }

    Vector zero(weights.get_size(),0);
    for (Index i = 0; i < weights.get_number(); ++i)
    {
        if (weights[i] < zero) { return false; }
    }
    return true;
}

bool
WeightAlgorithm::get_weights(
                const VectorArray& matrix,
                const VectorArray& lattice,
                const BitSet& urs,
                VectorArray& weights)
{
    weights.renumber(0);
    // First check to see if all ones vector is in the row span of the matrix.
    Vector all_ones(lattice.get_size());
    for (Index i = 0; i < all_ones.get_size(); ++i)
    {
        if (urs[i] == 1) { all_ones[i] = 0; }
        else { all_ones[i] = 1; }
    }
    Vector result(lattice.get_number());
    VectorArray::dot(lattice, all_ones, result);
    if (result.is_zero())
    {
        weights.insert(all_ones);
        return true;
    }

    // Look for a set of vectors which together imply a weighting.
    BitSet mask(matrix.get_size(),false);
    while (mask.count()<mask.get_size()-urs.count() &&
           get_weights(matrix,urs,mask,weights));
    if (mask.count() == mask.get_size()-urs.count()) { return true; }

    // We have not found a grading.
    weights.insert(all_ones);
    return false;
}

bool
WeightAlgorithm::get_weights(
                const VectorArray& matrix,
                const BitSet& urs,
                BitSet& mask,
                VectorArray& weights)
{
    int max = 0;
    int index = -1;
    for (int i = 0; i < matrix.get_number(); ++i)
    {
        if (is_candidate(matrix[i], urs, mask))
        {
            int count = positive_count(matrix[i], mask);
            if (count > max)
            {
                index = i;
                max = count;
            }
        }
    }
    if (index != -1)
    {
        weights.insert(matrix[index]);
        update_mask(mask, matrix[index]);
        return true;
    }
    return false;
}

bool
WeightAlgorithm::is_candidate(
                const Vector& v,
                const BitSet& urs,
                const BitSet& mask)
{
    assert(v.get_size() == mask.get_size());
    for (int i = 0; i < v.get_size(); ++i)
    {
        if (mask[i] == 0 && v[i] < 0) { return false; }
        if (urs[i] == 1 && v[i] != 0) { return false; }
    }
    return true;
}

bool
WeightAlgorithm::violates_urs(
                const Vector& v,
                const BitSet& urs)
{
    assert(v.get_size() == urs.get_size());
    for (int i = 0; i < v.get_size(); ++i)
    {
        if (urs[i] == 1 && v[i] != 0) { return true; }
    }
    return false;
}

int
WeightAlgorithm::positive_count(const Vector& v, const BitSet& mask)
{
    assert(v.get_size() == mask.get_size());
    int count = 0;
    for (int i = 0; i < v.get_size(); ++i)
    {
        if (mask[i] == 0 && v[i] > 0) { ++count; }
    }
    return count;
}

void
WeightAlgorithm::update_mask(BitSet& mask, const Vector& v)
{
    assert(v.get_size() == mask.get_size());
    for (int i = 0; i < v.get_size(); ++i)
    {
        if (v[i] > 0) { mask.set(i); }
    }
}

void
WeightAlgorithm::strip_weights(
                VectorArray* weights,
                Weight* max_weights,
                const BitSet& urs)
{
    if (max_weights == 0 || weights == 0) { return; }
    if (weights->get_number() == 0) { return; }
    assert(weights->get_number() == max_weights->get_size());
    BitSet good(max_weights->get_size(), true);
    Vector zero(weights->get_size(), 0);
    for (Index i = weights->get_number()-1; i >= 0; --i)
    {
        if ((*weights)[i] < zero || violates_urs((*weights)[i], urs))
        {
            weights->remove(i);
            good.unset(i);
        }
    }
    max_weights->project(good);
}
