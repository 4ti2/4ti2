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
#include "groebner/LatticeBasis.h"
#include "groebner/Globals.h"

#include "groebner/Debug.h"

using namespace _4ti2_;

QSolveAlgorithm::QSolveAlgorithm()
{
#ifdef _4ti2_GMP_
    variant = SUPPORT;
#else
    variant = MATRIX; 
#endif
    order = MAXCUTOFF;
}

QSolveAlgorithm::QSolveAlgorithm(QSolveVariant v, QSolveConsOrder o)
    : variant(v), order(o)
{
}

QSolveAlgorithm::~QSolveAlgorithm()
{
}

BitSet
QSolveAlgorithm::compute(
                const VectorArray& matrix,
                VectorArray& vs,
                VectorArray& subspace,
                const Vector& rel,
                const Vector& sign)
{
    Index extra_cols = 0;
    for (Index i = 0; i < rel.get_size(); ++i) { if (rel[i] != 0 && rel[i] != 3) { ++extra_cols; } }
    if (extra_cols == 0) {
        BitSet rs(sign.get_size(), false);
        BitSet cirs(sign.get_size(), false);
        convert_sign(sign, rs, cirs);
        if (!cirs.empty()) {
            std::cerr << "ERROR: Circuits components not supported.\n;";
            exit(1);
        }
        lattice_basis(matrix, vs);
        return compute(matrix, vs, subspace, rs);      
    }
    else {
        VectorArray ext_matrix(matrix.get_number(),matrix.get_size()+extra_cols,0);
        VectorArray ext_vs(0,vs.get_size()+extra_cols,0);
        VectorArray ext_subspace(0,subspace.get_size()+extra_cols,0);
        Vector ext_sign(matrix.get_size()+extra_cols,0);

        VectorArray::lift(matrix, 0, matrix.get_size(), ext_matrix);
        Vector::lift(sign, 0, sign.get_size(), ext_sign);
        Index ext = matrix.get_size();
        for (Index i = 0; i < matrix.get_number(); ++i) {
            if (rel[i] == 1) {
                ext_matrix[i][ext] = -1;
                ext_sign[ext] = 1;
                ++ext;
            }
            else if (rel[i] == -1) {
                ext_matrix[i][ext] = 1;
                ext_sign[ext] = 1;
                ++ext;
            }
            else if (rel[i] == 2) {
                std::cerr << "ERROR: Circuit components not supported.\n";
                exit(1);
            }
        }

        lattice_basis(ext_matrix, ext_vs);

        BitSet ext_rs(ext_sign.get_size(), false);
        BitSet ext_cirs(ext_sign.get_size(), false);
        convert_sign(ext_sign, ext_rs, ext_cirs);
        if (!ext_cirs.empty()) {
                std::cerr << "ERROR: Circuit components not supported.\n";
                exit(1);
        }
        BitSet ext_result(ext_matrix.get_size(), false);
        ext_result = compute(ext_matrix, ext_vs, ext_subspace, ext_rs);
        BitSet result(matrix.get_size(), false);
        BitSet::shrink(ext_result, result);

        vs.renumber(ext_vs.get_number());
        VectorArray::project(ext_vs, 0, vs.get_size(), vs);
        subspace.renumber(ext_subspace.get_number());
        VectorArray::project(ext_subspace, 0, subspace.get_size(), subspace);

        return result;
    }
}


// TODO: This function is essentially a copy of the one above. 
// The code should be refactored.
void
QSolveAlgorithm::compute(
                const VectorArray& matrix,
                VectorArray& vs,
                VectorArray& circuits,
                VectorArray& subspace,
                const Vector& rel,
                const Vector& sign)
{
    Index extra_cols = 0;
    for (Index i = 0; i < rel.get_size(); ++i) { if (rel[i] != 0 && rel[i] != 3) { ++extra_cols; } }
    if (extra_cols == 0) {
        BitSet rs(sign.get_size(), false);
        BitSet cirs(sign.get_size(), false);
        convert_sign(sign, rs, cirs);
        lattice_basis(matrix, vs);
        compute(matrix, vs, circuits, subspace, rs, cirs);      
    }
    else {
        VectorArray ext_matrix(matrix.get_number(),matrix.get_size()+extra_cols,0);
        VectorArray ext_vs(0,vs.get_size()+extra_cols,0);
        VectorArray ext_circuits(0,circuits.get_size()+extra_cols,0);
        VectorArray ext_subspace(0,subspace.get_size()+extra_cols,0);
        Vector ext_sign(matrix.get_size()+extra_cols,0);

        VectorArray::lift(matrix, 0, matrix.get_size(), ext_matrix);
        Vector::lift(sign, 0, sign.get_size(), ext_sign);
        Index ext = matrix.get_size();
        for (Index i = 0; i < matrix.get_number(); ++i) {
            if (rel[i] == 1) {
                ext_matrix[i][ext] = -1;
                ext_sign[ext] = 1;
                ++ext;
            }
            else if (rel[i] == 2) {
                ext_matrix[i][ext] = -1;
                ext_sign[ext] = 2;
                ++ext;
            }
            else if (rel[i] == -1) {
                ext_matrix[i][ext] = 1;
                ext_sign[ext] = 1;
                ++ext;
            }
        }

        lattice_basis(ext_matrix, ext_vs);

        BitSet ext_rs(ext_sign.get_size(), false);
        BitSet ext_cirs(ext_sign.get_size(), false);
        convert_sign(ext_sign, ext_rs, ext_cirs);

        DEBUG_4ti2(*out << "MATRIX:\n" << matrix << "\n";)
        DEBUG_4ti2(*out std::cout << "EXT MATRIX:\n" << ext_matrix << "\n";)

        compute(ext_matrix, ext_vs, ext_circuits, ext_subspace, ext_rs, ext_cirs);

        vs.renumber(ext_vs.get_number());
        VectorArray::project(ext_vs, 0, vs.get_size(), vs);
        subspace.renumber(ext_subspace.get_number());
        VectorArray::project(ext_subspace, 0, subspace.get_size(), subspace);
        circuits.renumber(ext_circuits.get_number());
        VectorArray::project(ext_circuits, 0, circuits.get_size(), circuits);
    }
}

void
QSolveAlgorithm::convert_sign(const Vector& sign, BitSet& rs, BitSet& cirs)
{
    assert(sign.get_size() == rs.get_size() && sign.get_size() == cirs.get_size());
    for (Index i = 0; i < sign.get_size(); ++i) {
        if (sign[i] == 1) { rs.set(i); }
        else if (sign[i] == 2) { cirs.set(i); }
        else if (sign[i] == -1) {
            std::cerr << "ERROR: non-positive variables not yet supported.\n";
            exit(1);
        }
    }
}

BitSet
QSolveAlgorithm::compute(
                const VectorArray& matrix,
                VectorArray& vs,
                VectorArray& subspace,
                const BitSet& rs)
{
    // Remove the linear subspace if there is one.
    linear_subspace(matrix, vs, rs, subspace);
    if (!subspace.get_number()) {
        return compute(matrix, vs, rs);
    }
    else {
        // We insert the linear subspace into the matrix to make the cone
        // pointed.
        VectorArray ext_matrix(matrix);
        ext_matrix.insert(subspace);
        return compute(ext_matrix, vs, rs);
    }
}

BitSet
QSolveAlgorithm::compute(
                const VectorArray& matrix,
                VectorArray& vs,
                const BitSet& rs)
{
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
                    const VectorArray& matrix,
                    VectorArray& vs,
                    const BitSet& rs,
                    VectorArray& subspace)
{
    subspace.renumber(0);
    if (rs.count() == matrix.get_size()) { return; }
    Index rows = upper_triangle(vs, rs);
    VectorArray::transfer(vs, rows, vs.get_number(), subspace, 0);
    rows = upper_triangle(subspace);
    if (rows != 0) {
        *out << "Cone is not pointed.\n";
        subspace.remove(rows, subspace.get_number());
    }
}

void
QSolveAlgorithm::compute(
                const VectorArray& matrix,
                VectorArray& vs,
                VectorArray& circuits,
                VectorArray& subspace,
                const BitSet& rs,
                const BitSet& cirs)
{
    // Remove the linear subspace if there is one.
    linear_subspace(matrix, vs, rs, cirs, subspace);
    if (!subspace.get_number()) {
        compute(matrix, vs, circuits, rs, cirs);
    }
    else {
        // We insert the linear subspace into the matrix to make the cone
        // pointed.
        VectorArray ext_matrix(matrix);
        ext_matrix.insert(subspace);
        compute(ext_matrix, vs, circuits, rs, cirs);
    }
}

void
QSolveAlgorithm::compute(
                const VectorArray& matrix,
                VectorArray& vs,
                VectorArray& circuits,
                const BitSet& rs,
                const BitSet& cirs)
{
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
                    const VectorArray& matrix,
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
    if (rows != 0) {
        *out << "Cone is not pointed.\n";
        subspace.remove(rows, subspace.get_number());
    }
}


