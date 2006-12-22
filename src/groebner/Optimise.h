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

#ifndef _4ti2__Optimise_
#define _4ti2__Optimise_

#include "VectorArray.h"
#include "Weight.h"
#include "BitSet.h"
#include "Feasible.h"

namespace _4ti2_
{

class Optimise
{
public:
    Optimise();
    virtual ~Optimise();

    virtual int compute(
                    Feasible& feasible,
                    Vector& cost,
                    Vector& sol);

protected:
    int compute_infeasible(
                    Feasible& feasible,
                    Vector& cost,
                    Vector& sol);
    int compute_feasible(
                    Feasible& feasible,
                    Vector& cost,
                    Vector& sol);
    int compute_feasible(
                    Feasible& feasible,
                    int cost_col,
                    IntegerType cost_offset,
                    Vector& sol);
    int compute_bounded(
                    Feasible& feasible,
                    int cost_col,
                    IntegerType cost_offset,
                    Vector& sol,
                    IntegerType upper_bound);

    int next_support(
                    const VectorArray& gens,
                    const BitSet& fin,
                    const Vector& sol);

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

} // namespace _4ti2_

#endif
