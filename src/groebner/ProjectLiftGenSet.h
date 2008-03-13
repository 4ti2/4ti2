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

#ifndef _4ti2_groebner__ProjectLiftGenSet_
#define _4ti2_groebner__ProjectLiftGenSet_

#include "groebner/VectorArray.h"
#include "groebner/Weight.h"
#include "groebner/BitSet.h"
#include "groebner/Feasible.h"

namespace _4ti2_
{

class ProjectLiftGenSet
{
public:
    ProjectLiftGenSet();
    virtual ~ProjectLiftGenSet();

    virtual void compute(
                    Feasible& feasible,
                    VectorArray& gens,
                    bool minimal = true);
    virtual void compute(
                    Feasible& feasible,
                    VectorArray& gens,
                    VectorArray& feasibles,
                    bool minimal = true);

protected:
    virtual void compute_bounded(
                    Feasible& feasible,
                    VectorArray& gens,
                    VectorArray& feasibles,
                    bool minimal = true);
    virtual void compute_unbounded(
                    Feasible& feasible,
                    VectorArray& gens,
                    VectorArray& feasibles,
                    bool minimal = true);

    int next_support(
                    const VectorArray& gens,
                    BitSet& fin);

    int add_support(
                    const VectorArray& gens,
                    BitSet& fin);

    int positive_count(
                    const VectorArray& gens,
                    int c);

    void make_feasible(
                    VectorArray& feasibles,
                    const Vector& ray);
};

inline
void
ProjectLiftGenSet::compute(
                    Feasible& feasible,
                    VectorArray& gens,
                    bool minimal)
{
    VectorArray feasibles(0,feasible.get_dimension());
    compute(feasible, gens, feasibles, minimal);
}

} // namespace _4ti2_

#endif
