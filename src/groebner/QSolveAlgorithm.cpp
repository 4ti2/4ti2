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

#include "groebner/QSolveAlgorithm.h"
#include "groebner/CircuitMatrixAlgorithm.h"
#include "groebner/CircuitSupportAlgorithm.h"
#include "groebner/RayMatrixAlgorithm.h"
#include "groebner/RaySupportAlgorithm.h"
#include "groebner/LongDenseIndexSet.h"
#include "groebner/ShortDenseIndexSet.h"
#include "groebner/IndexSetConverter.h"
#include "groebner/Globals.h"

#include "groebner/Debug.h"

using namespace _4ti2_;

QSolveAlgorithm::QSolveAlgorithm()
{
    order = MAXCUTOFF;
#ifdef _4ti2_GMP_
    variant = SUPPORT;
#else
    variant = MATRIX; 
#endif
}

QSolveAlgorithm::QSolveAlgorithm(QSolveConsOrder o, QSolveVariant v)
    : order(o), variant(v)
{
}

QSolveAlgorithm::~QSolveAlgorithm()
{
}

BitSet
QSolveAlgorithm::compute(
                VectorArray& matrix,
                VectorArray& vs,
                VectorArray& subspace,
                const BitSet& rs)
{
    // We remove the linear subspace if there is one.
    linear_subspace(matrix, vs, rs, subspace);
    BitSet result(rs.get_size());
    if (variant == SUPPORT) {
        if  (rs.get_size() <= ShortDenseIndexSet::max_size) {
            DEBUG_4ti2(*out << "Using Short BitSet.\n";)
            ShortDenseIndexSet rs_tmp(rs.get_size());
            convert(rs, rs_tmp);
            RaySupportAlgorithm<ShortDenseIndexSet> algorithm(order);
            ShortDenseIndexSet tmp_result = algorithm.compute(matrix, vs, rs_tmp);
            convert(tmp_result, result);
        }
        else {
            DEBUG_4ti2(*out << "Using Long BitSet.\n";)
            RaySupportAlgorithm<BitSet> algorithm(order);
            result = algorithm.compute(matrix, vs, rs);
        }
    }
    else { //if (variant == MATRIX)
        if (rs.get_size() <= ShortDenseIndexSet::max_size) {
            DEBUG_4ti2(*out << "Using Short BitSet.\n";)
            ShortDenseIndexSet rs_tmp(rs.get_size());
            convert(rs, rs_tmp);
            RayMatrixAlgorithm<ShortDenseIndexSet> algorithm(order);
            ShortDenseIndexSet tmp_result = algorithm.compute(matrix, vs, rs_tmp);
            convert(tmp_result, result);
        }
        else {
            DEBUG_4ti2(*out << "Using Long BitSet.\n";)
            RayMatrixAlgorithm<BitSet> algorithm(order);
            result = algorithm.compute(matrix, vs, rs);
        }
    }
    return result;
}

void
QSolveAlgorithm::linear_subspace(
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

void
QSolveAlgorithm::compute(
                VectorArray& matrix,
                VectorArray& vs,
                VectorArray& circuits,
                VectorArray& subspace,
                const BitSet& rs,
                const BitSet& cirs)
{
    // Remove the linear subspace if there is one.
    linear_subspace(matrix, vs, rs, cirs, subspace);
    if (variant == SUPPORT) {
        if (cirs.get_size()+cirs.count() <= ShortDenseIndexSet::max_size) {
            DEBUG_4ti2(*out << "Using Short BitSet.\n";)
            ShortDenseIndexSet cirs_tmp(cirs.get_size());
            convert(cirs, cirs_tmp);
            ShortDenseIndexSet rs_tmp(cirs.get_size());
            convert(rs, rs_tmp);
            CircuitSupportAlgorithm<ShortDenseIndexSet> algorithm;
            algorithm.compute(matrix, vs, circuits, rs_tmp, cirs_tmp);
        }
        else {
            DEBUG_4ti2(*out << "Using Long BitSet.\n";)
            CircuitSupportAlgorithm<BitSet> algorithm;
            algorithm.compute(matrix, vs, circuits, rs, cirs);
        }
    }
    else { // if (variant == MATRIX)
        if (cirs.get_size() <= ShortDenseIndexSet::max_size) {
            DEBUG_4ti2(*out << "Using Short BitSet.\n";)
            ShortDenseIndexSet cirs_tmp(cirs.get_size());
            convert(cirs, cirs_tmp);
            ShortDenseIndexSet rs_tmp(cirs.get_size());
            convert(rs, rs_tmp);
            CircuitMatrixAlgorithm<ShortDenseIndexSet> algorithm;
            algorithm.compute(matrix, vs, circuits, rs_tmp, cirs_tmp);
        }
        else {
            DEBUG_4ti2(*out << "Using Long BitSet.\n";)
            CircuitMatrixAlgorithm<BitSet> algorithm;
            algorithm.compute(matrix, vs, circuits, rs, cirs);
        }
    }
}

void
QSolveAlgorithm::linear_subspace(
                    VectorArray& matrix,
                    VectorArray& vs,
                    const BitSet& rs,
                    const BitSet& cirs,
                    VectorArray& subspace)
{
    assert(BitSet::set_disjoint(rs,cirs));
    if (rs.count()+cirs.count() == matrix.get_size()) { return; }
    Index rs_rows = upper_triangle(vs, rs);
    Index cirs_rows = upper_triangle(vs, cirs, rs_rows);
    subspace.renumber(0);
    VectorArray::transfer(vs, cirs_rows, vs.get_number(), subspace, 0);
    Index rows = upper_triangle(subspace);
    if (rows != 0)
    {
        *out << "Cone is not pointed.\n";
        subspace.remove(rows, subspace.get_number());
        // We insert the linear subspace into the matrix to make the cone
        // pointed.
        matrix.insert(subspace);
    }
}


