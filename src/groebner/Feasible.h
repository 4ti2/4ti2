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

#ifndef _4ti2__Feasible_
#define _4ti2__Feasible_

#include "VectorArray.h"
#include "BitSet.h"

namespace _4ti2_
{

class Feasible
{
public:
    Feasible(   const VectorArray* basis,
                const VectorArray* matrix = 0,
                const BitSet* urs = 0,
                const Vector* rhs = 0,
                const VectorArray* weights = 0,
                const Vector* max_weights = 0);
    Feasible(const Feasible& feasible);
    Feasible(const Feasible& feasible, const BitSet& urs);
    Feasible& operator=(const Feasible& feasible);
    ~Feasible();

    int get_dimension();

    const VectorArray& get_matrix();
    const VectorArray& get_basis();
    const BitSet& get_urs();
    const Vector* get_rhs();
    const VectorArray* get_weights();
    const Vector* get_max_weights();

    // Is the cost function bounded?
    bool bounded(const VectorArray& cost, BitSet& cost_unbnd);

    void add_inequalities(const VectorArray& constr, const Vector& rhs);

    // The properties of the feasible sets.
    const BitSet& get_bnd();
    const BitSet& get_unbnd();
    const Vector& get_grading();
    const Vector& get_ray();

protected:
    void compute_bounded();
    void check_weights();

    int dim;

    VectorArray* basis;
    VectorArray* matrix;
    BitSet* urs;
    Vector* rhs;
    VectorArray* weights;
    Vector* max_weights;

    bool computed_bounded;
    BitSet* bnd;
    BitSet* unbnd;
    Vector* grading;
    Vector* ray;
};

inline
int
Feasible::get_dimension()
{
    return dim;
}

inline
const VectorArray&
Feasible::get_matrix()
{
    assert(matrix != 0);
    return *matrix;
}

inline
const VectorArray&
Feasible::get_basis()
{
    assert(basis != 0);
    return *basis;
}

inline
const BitSet&
Feasible::get_urs()
{
    assert(urs != 0);
    return *urs;
}

inline
const Vector*
Feasible::get_rhs()
{
    return rhs;
}

inline
const VectorArray*
Feasible::get_weights()
{
    return weights;
}

inline
const Vector*
Feasible::get_max_weights()
{
    return max_weights;
}

inline
const BitSet&
Feasible::get_bnd()
{
    compute_bounded();
    return *bnd;
}

inline
const BitSet&
Feasible::get_unbnd()
{
    compute_bounded();
    return *unbnd;
}

inline
const Vector&
Feasible::get_grading()
{
    compute_bounded();
    return *grading;
}

inline
const Vector&
Feasible::get_ray()
{
    compute_bounded();
    return *ray;
}

} // namespace _4ti2_

#endif
