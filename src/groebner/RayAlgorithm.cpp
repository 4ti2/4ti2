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

#include "RayAlgorithm.h"
#include "RayImplementation.h"
#include "RayMatrixAlgorithm.h"
#include "RaySupportAlgorithm.h"
#include "CircuitOptions.h"
#include "LongDenseIndexSet.h"
#include "ShortDenseIndexSet.h"
#include "IndexSetConverter.h"
#include "Globals.h"

#include "Debug.h"

using namespace _4ti2_;

RayAlgorithm::RayAlgorithm()
{
}

RayAlgorithm::~RayAlgorithm()
{
}

BitSet
RayAlgorithm::compute(
                VectorArray& matrix,
                VectorArray& vs,
                VectorArray& subspace,
                const BitSet& rs)
{
    // We remove the linear subspace if there is one.
    linear_subspace(matrix, vs, rs, subspace);
    BitSet result(rs.get_size());
    if (CircuitOptions::instance()->algorithm == CircuitOptions::SUPPORT)
    {
        if  (rs.get_size() <= ShortDenseIndexSet::max_size)
        {
            DEBUG_4ti2(*out << "Using Short BitSet.\n";)
            ShortDenseIndexSet rs_tmp(rs.get_size());
            convert(rs, rs_tmp);
            RaySupportAlgorithm<ShortDenseIndexSet> algorithm;
            ShortDenseIndexSet tmp_result = algorithm.compute(matrix, vs, rs_tmp);
            convert(tmp_result, result);
        }
        else
        {
            DEBUG_4ti2(*out << "Using Long BitSet.\n";)
            RaySupportAlgorithm<BitSet> algorithm;
            result = algorithm.compute(matrix, vs, rs);
        }
    }
    else //if (CircuitOptions::instance()->algorithm == CircuitOptions::MATRIX)
    {
        if (rs.get_size() <= ShortDenseIndexSet::max_size)
        {
            DEBUG_4ti2(*out << "Using Short BitSet.\n";)
            ShortDenseIndexSet rs_tmp(rs.get_size());
            convert(rs, rs_tmp);
            RayMatrixAlgorithm<ShortDenseIndexSet> algorithm;
            ShortDenseIndexSet tmp_result = algorithm.compute(matrix, vs, rs_tmp);
            convert(tmp_result, result);
        }
        else
        {
            DEBUG_4ti2(*out << "Using Long BitSet.\n";)
            RayMatrixAlgorithm<BitSet> algorithm;
            result = algorithm.compute(matrix, vs, rs);
        }
    }

    return result;
}

void
RayAlgorithm::linear_subspace(
                    VectorArray& matrix,
                    VectorArray& vs,
                    const BitSet& rs,
                    VectorArray& subspace)
{
    subspace.renumber(0);
    if (rs.count() == matrix.get_size()) { return; }
    Index rows = upper_triangle(vs, rs);
    VectorArray::transfer(vs, rows, vs.get_number(), subspace, 0);
    rows = upper_triangle(subspace);
    if (rows != 0)
    {
        *out << "Cone is not pointed.\n";
        subspace.remove(rows, subspace.get_number());
        // We insert the linear subspace into the matrix to make the cone
        // pointed.
        matrix.insert(subspace);
    }
}
