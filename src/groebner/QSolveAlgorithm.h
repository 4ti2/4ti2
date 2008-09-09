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

#ifndef _4ti2_groebner__QSolveAlgorithm_
#define _4ti2_groebner__QSolveAlgorithm_

#include "groebner/VectorArray.h"
#include "groebner/BitSet.h"
#include "groebner/QSolveConsOrder.h"
#include "groebner/QSolveVariant.h"

namespace _4ti2_
{

class QSolveAlgorithm
{
public:
    QSolveAlgorithm();
    QSolveAlgorithm(QSolveVariant v, QSolveConsOrder o);
    ~QSolveAlgorithm();

    void set_constraint_order(QSolveConsOrder o);
    void set_variant(QSolveVariant v);

    BitSet compute(
                    VectorArray& matrix,
                    VectorArray& vs,
                    VectorArray& subspace,
                    Vector& rels,
                    Vector& sign);

    void compute(
                    VectorArray& matrix,
                    VectorArray& vs,
                    VectorArray& circuits,
                    VectorArray& subspace,
                    Vector& rels,
                    Vector& sign);

    BitSet compute(
                    VectorArray& matrix,
                    VectorArray& vs,
                    VectorArray& subspace,
                    const BitSet& rs);

    void compute(
                    VectorArray& matrix,
                    VectorArray& vs,
                    VectorArray& circuits,
                    VectorArray& subspace,
                    const BitSet& rs,
                    const BitSet& cirs);

protected:
    void linear_subspace(
                    VectorArray& matrix,
                    VectorArray& vs,
                    const BitSet& rs,
                    const BitSet& cirs,
                    VectorArray& subspace);

    void linear_subspace(
                    VectorArray& matrix,
                    VectorArray& vs,
                    const BitSet& rs,
                    VectorArray& subspace);

    void convert_sign(const Vector& sign, BitSet& rs, BitSet& cirs);

    QSolveVariant variant;
    QSolveConsOrder order;
};

inline
void
QSolveAlgorithm::set_constraint_order(QSolveConsOrder o) 
{
    order = o;
}

inline
void
QSolveAlgorithm::set_variant(QSolveVariant v)
{
    variant = v;
}

} // namespace _4ti2_

#endif
