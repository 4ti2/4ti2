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

#include "GroebnerBasis.h"
#include "Completion.h"
#include "ProjectLiftGenSet.h"
#include "SaturationGenSet.h"
#include "Globals.h"

#include <iostream>
#include <iomanip>
#include "VectorStream.h"
#include "VectorArrayStream.h"
#include "BitSetStream.h"

#define DEBUG_4ti2_GroebnerBasis(X) //X

using namespace _4ti2_;

GroebnerBasis::GroebnerBasis(
                Feasible& _feasible,
                VectorArray* _cost,
                VectorArray* gb)
{
    feasible = &_feasible;
    if (_cost != 0) { cost = new VectorArray(*_cost); }
    else { cost = new VectorArray(0,feasible->get_dimension()); }

    if (gb != 0) { gens = new VectorArray(*gb); }
    else 
    {
        gb = new VectorArray(feasible->get_basis());
        GeneratingSet::compute();
        compute();
    }
}

GroebnerBasis::GroebnerBasis(
                GeneratingSet& gs,
                VectorArray* _cost)
{
    feasible = &gs.get_feasible();
    gens = new VectorArray(gs.get_generating_set());

    if (_cost != 0) { cost = new VectorArray(*_cost); }
    else { cost = new VectorArray(0,feasible->get_dimension()); }

    compute();
}

GroebnerBasis::GroebnerBasis(
                GroebnerBasis& gb,
                VectorArray* _cost)
{
    feasible = &gb.get_feasible();
    gens = new VectorArray(gb.get_groebner_basis());

    if (_cost != 0) { cost = new VectorArray(*_cost); }
    else { cost = new VectorArray(0,feasible->get_dimension()); }

    compute();
}

GroebnerBasis::~GroebnerBasis()
{
    delete cost;
}

const VectorArray&
GroebnerBasis::get_groebner_basis()
{
    assert(gens != 0);
    return *gens;
}

const VectorArray&
GroebnerBasis::get_cost()
{
    assert(cost != 0);
    return *cost;
}

void
GroebnerBasis::compute()
{
    assert(cost != 0 && gens != 0 && feasible != 0);
    Completion algorithm;
    algorithm.compute(*feasible, *cost, *gens);
    gens->sort();
}
