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

#include "Feasible.h"
#include "Bounded.h"
#include "LatticeBasis.h"
#include "WeightAlgorithm.h"
#include "VectorArrayStream.h"
#include "Globals.h"

//#define DEBUG_4ti2(X) X
#include "Debug.h"

using namespace _4ti2_;

Feasible::Feasible(
                const VectorArray* _basis,
                const VectorArray* _matrix,
                const BitSet* _urs,
                const Vector* _rhs,
                const VectorArray* _weights,
                const Vector* _max_weights)
{
    assert(_matrix != 0 || _basis != 0);
    // Set the dimension of the problem.
    if (_matrix != 0) { dim = _matrix->get_size(); }
    else { dim = _basis->get_size(); }

    assert(_matrix == 0 || _matrix->get_size() == dim);
    assert(_basis == 0 || _basis->get_size() == dim);
    assert(_urs == 0 || _urs->get_size() == dim);
    assert(_rhs == 0 || _rhs->get_size() == dim);
    assert(_weights == 0 || _weights->get_size() == dim);
    assert(_max_weights == 0 || _max_weights->get_size() == dim);

    basis = new VectorArray(0, dim);
    matrix = new VectorArray(0, dim);
    urs = new BitSet(dim);

    if (_basis != 0) { *basis = *_basis; }
    else { lattice_basis(*_matrix, *basis); }

    if (_matrix != 0) { *matrix = *_matrix; }
    else { lattice_basis(*_basis, *matrix); }
    DEBUG_4ti2(*out << "Matrix:\n" << *matrix << "\n";)
    DEBUG_4ti2(*out << "Basis:\n" << *basis << "\n";)

    if (_urs != 0) { *urs = *_urs; }

    rhs = 0;
    weights = 0;
    max_weights = 0;
    if (_rhs != 0) { rhs = new Vector(*_rhs); }
    if (_weights != 0) { weights = new VectorArray(*_weights); }
    if (_max_weights != 0) { max_weights = new Vector(*_max_weights); }
    WeightAlgorithm::strip_weights(weights, max_weights, *urs);

    computed_bounded = false;
    bnd = 0;
    unbnd = 0;
    grading = 0;
    ray = 0;
}

Feasible::Feasible(const Feasible& feasible)
{
    *this = feasible;
}

Feasible::Feasible(const Feasible& feasible, const BitSet& _urs)
{
    *this = feasible;

    if (*urs == _urs) { return; }
    else { computed_bounded = false; }

    // If some of the bounded variables have become urs, then the set of bounded
    // variables is invalid.
    if (bnd != 0 && !BitSet::set_disjoint(*bnd, _urs))
    {
        delete bnd; bnd = 0;
        delete grading; grading = 0;
    }

    if (unbnd != 0) {
        // If some of the urs variables have become non-negative, then the set
        // of unbounded variables is invalid.
        if (!BitSet::set_subset(*urs, _urs)) {
            delete unbnd; unbnd = 0;
            delete ray; ray = 0;
        }
        else { unbnd->set_difference(_urs); }
    }
    *urs = _urs;
    WeightAlgorithm::strip_weights(weights, max_weights, *urs);
}

Feasible&
Feasible::operator=(const Feasible& feasible)
{
    dim = feasible.dim;
    basis = new VectorArray(*feasible.basis);
    matrix = new VectorArray(*feasible.matrix);
    urs = new BitSet(*feasible.urs);

    rhs = 0;
    weights = 0;
    max_weights = 0;
    if (feasible.rhs != 0) {
        rhs = new Vector(*feasible.rhs);
    }
    if (feasible.weights != 0) {
        weights = new VectorArray(*feasible.weights);
    }
    if (feasible.max_weights != 0) {
        max_weights = new Vector(*feasible.max_weights);
    }

    computed_bounded = feasible.computed_bounded;

    bnd = 0;
    unbnd = 0;
    grading = 0;
    ray = 0;
    if (feasible.bnd != 0) {
        bnd = new BitSet(*feasible.bnd);
    }
    if (feasible.unbnd != 0) {
        unbnd = new BitSet(*feasible.unbnd);
    }
    if (feasible.grading != 0) {
        grading = new Vector(*feasible.grading);
    }
    if (feasible.ray != 0) {
        ray = new Vector(*feasible.ray);
    }
    return *this;
}

Feasible::~Feasible()
{
    delete basis;
    delete matrix;
    delete urs;
    delete rhs;
    delete weights;
    delete max_weights;
    delete bnd;
    delete unbnd;
    delete grading;
    delete ray;
}

void
Feasible::compute_bounded()
{
    if (computed_bounded == true) { return; }
    if (bnd == 0) { bnd = new BitSet(dim); }
    if (unbnd == 0) { unbnd = new BitSet(dim); }
    if (grading == 0) { grading = new Vector(dim,0); }
    if (ray == 0) { ray = new Vector(dim,0); }

    _4ti2_::bounded(*matrix, *basis, *urs, *bnd, *grading, *unbnd, *ray);
    computed_bounded = true;
}

bool
Feasible::bounded(
                const VectorArray& cost,
                BitSet& cost_unbnd)
{
    compute_bounded();
    return _4ti2_::bounded(*matrix, *basis,
                    *urs, cost,
                    *bnd, *unbnd, cost_unbnd);
}

#if 0
void
Feasible::add_inequalities(const VectorArray& constr, const Vector& rhs)
{
    assert(constr.get_size() == dim);
    assert(constr.get_number() == rhs.get_size());

    // Update the dimension.
    dim = matrix->get_size()+constr.get_number();

    // Add the constraints to the matrix.
    assert(matrix != 0);
    int oldm = matrix->get_number();
    int oldn = matrix->get_size();
    int m = matrix->get_number()+constr.get_number();
    int n = matrix->get_size()+constr.get_number();

    matrix->insert(constr);
    VectorArray* ext_matrix = new VectorArray(m,n,0);
    VectorArray::lift(*matrix, 0, matrix->get_size(), *ext_matrix);

    for (int i = 0; i < constr.get_number(); ++i)
    {
        (*ext_matrix)[i+oldm][i+oldn] = 1;
    }
    delete matrix;
    matrix = ext_matrix;
    DEBUG_4ti2(*out << "Full Matrix:\n" << *matrix << "\n";)

    // Add the constraints to the basis.
    VectorArray* ext_basis = new VectorArray(basis->get_number(),n);

    VectorArray::lift(*basis, 0, basis->get_size(), *ext_basis);
    VectorArray slacks(lattice.get_number());
    VectorArray::dot(lattice, cost, cost_col);
    for (int i = 0; i < ext_lattice.get_number(); ++i)
    {   ext_lattice[i][lattice.get_size()] = -cost_col[i]; }
    DEBUG_4ti2(*out << "Full Lattice:\n" << ext_lattice << "\n";)

    // Extend the urs.
    const BitSet& urs = feasible.get_urs();
    BitSet ext_urs(urs.get_size()+1);
    BitSet::extend(urs, ext_urs);
    DEBUG_4ti2(*out << "Extended URS:\n" << ext_urs << "\n";)
}
#endif
