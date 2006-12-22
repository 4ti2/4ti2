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

#ifndef _4ti2__Completion_
#define _4ti2__Completion_

#include "VectorArray.h"
#include "BitSet.h"
#include "Weight.h"
#include "Timer.h"
#include "Feasible.h"

namespace _4ti2_
{

class Algorithm;

class Completion
{
public:
    Completion();
    ~Completion();

    // Compute test set.
    void compute(
                    Feasible& feasible,
                    const VectorArray& cost,
                    VectorArray& vs);
    // Compute test set and optimize feasibles.
    void compute(
                    Feasible& feasible,
                    const VectorArray& cost,
                    VectorArray& vs,
                    VectorArray& feasibles);

    // Compute test set.
    void compute(
                    Feasible& feasible,
                    const VectorArray& cost,
                    const BitSet& sat,
                    VectorArray& vs);
    // Compute test set and optimize feasibles.
    void compute(
                    Feasible& feasible,
                    const VectorArray& cost,
                    const BitSet& sat,
                    VectorArray& vs,
                    VectorArray& feasibles);

    void set_algorithm(Algorithm* _alg);

protected:
    Timer t;
    Algorithm* alg;
};

inline
void
Completion::compute(
                Feasible& feasible,
                const VectorArray& cost,
                VectorArray& vs)
{
    VectorArray feasibles(0,feasible.get_dimension());
    compute(feasible, cost, vs, feasibles);
}

inline
void
Completion::compute(
                Feasible& feasible,
                const VectorArray& cost,
                const BitSet& sat,
                VectorArray& vs)
{
    VectorArray feasibles(0,feasible.get_dimension());
    compute(feasible, cost, sat, vs, feasibles);
}

} // namespace _4ti2_

#endif
