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

#include "LatticeBasis.h"
#include "HermiteAlgorithm.h"
#include "BitSet.h"
#include "Globals.h"

#include "VectorArrayStream.h"
#include "VectorStream.h"

//#define DEBUG_4ti2(X) X
#include "Debug.h"

// Inputs a matrix and outputs a lattice basis in the same set.
void
_4ti2_::lattice_basis(const VectorArray& matrix, VectorArray& temp)
{
    DEBUG_4ti2(*out << "Computing Lattice basis ...\n";)
    VectorArray trans(matrix.get_size(), matrix.get_number());
    VectorArray::transpose(matrix, trans);
    DEBUG_4ti2(*out << "Transpose:\n" << trans << "\n";)

    VectorArray basis(matrix.get_size(), matrix.get_size(), 0);
    for (Index i = 0; i < basis.get_number(); ++i) { basis[i][i] = 1; }
    DEBUG_4ti2(*out << "Identity:\n" << basis << "\n";)

    VectorArray trans_basis(trans.get_number(),
                        trans.get_size()+basis.get_size());
    VectorArray::concat(trans, basis, trans_basis);
    DEBUG_4ti2(
        *out << "Transpose + Identity:\n" << trans_basis << "\n";)

    Index rows = upper_triangle(trans_basis, trans_basis.get_number(), trans.get_size());
    DEBUG_4ti2(
        *out << "Full Hermite:\n" << trans_basis << "\n";)

    VectorArray::project(trans_basis, trans.get_size(),
                    trans_basis.get_size(), basis);
    basis.remove(0,rows);
    DEBUG_4ti2(*out << "Basis:\n" << basis << "\n";)

    temp = basis;
}

// Inputs a matrix and outputs a lattice basis in the same set.
IntegerType
_4ti2_::solve(const VectorArray& matrix, const Vector& rhs, Vector& solution)
{
    VectorArray trans(matrix.get_size(), matrix.get_number());
    VectorArray::transpose(matrix, trans);
    Vector rhs_negative(rhs);
    rhs_negative.mul(-1);
    trans.insert(rhs_negative);
    DEBUG_4ti2(*out << "Transpose:\n" << trans << "\n";)

    VectorArray basis(matrix.get_size()+1, matrix.get_size()+1, 0);
    for (Index i = 0; i < basis.get_number(); ++i) { basis[i][i] = 1; }
    DEBUG_4ti2(*out << "Identity:\n" << basis << "\n";)

    VectorArray trans_basis(trans.get_number(),
                        trans.get_size()+basis.get_size());
    VectorArray::concat(trans, basis, trans_basis);
    DEBUG_4ti2(
        *out << "Transpose + Identity:\n" << trans_basis << "\n";)

    Index rows = upper_triangle(trans_basis, trans_basis.get_number(), trans.get_size());
    DEBUG_4ti2(
        *out << "Full Hermite:\n" << trans_basis << "\n";)

    VectorArray::project(trans_basis, trans.get_size(),
                    trans_basis.get_size(), basis);
    DEBUG_4ti2(
        *out << "Projected Hermite:\n" << basis << "\n";)
    basis.remove(0,rows);
    DEBUG_4ti2(*out << "Basis:\n" << basis << "\n";)

    BitSet proj(basis.get_size(), false);
    proj.set(basis.get_size()-1);
    upper_triangle(basis, proj);
    DEBUG_4ti2(
        *out << "Solution Hermite:\n" << basis << "\n";)
    if (basis.get_number() == 0)
    {
        solution.mul(0);
        return 0;
    }

    proj.set_complement();
    Vector::project(basis[0], proj, solution);
    DEBUG_4ti2(*out << "Solution:\n" << solution << "\n";)

    return basis[0][basis.get_size()-1];
}
