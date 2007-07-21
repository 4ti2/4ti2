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

#include "GeneratingSet.h"
#include "ProjectLiftGenSet.h"
#include "SaturationGenSet.h"
#include "HybridGenSet.h"
#include "Markov.h"
#include "Options.h"
#include "Globals.h"
#include "VectorArrayStream.h"
#include "Debug.h"

using namespace _4ti2_;

GeneratingSet::GeneratingSet()
{
    feasible = 0;
    gens = 0;
}

GeneratingSet::~GeneratingSet()
{
    delete gens;
}

GeneratingSet::GeneratingSet(
                Feasible& _feasible,
                VectorArray* _gens)
        : feasible(&_feasible)
{
    gens = _gens;
    if (gens == 0)
    {
        gens = new VectorArray(0, feasible->get_dimension());
        compute();
    }
}

void
GeneratingSet::compute()
{
    assert(gens != 0 && feasible != 0);
    if (Globals::generation == Globals::SATURATION)
    {
        SaturationGenSet algorithm;
        BitSet sat(feasible->get_dimension(), false);
        algorithm.compute(*feasible, *gens, sat, Globals::minimal);
    }
    else if (Globals::generation == Globals::PROJECT_AND_LIFT)
    {
#if 0
        // TODO: The following code can be used to compute feasible points.
        std::string vs_file(Options::instance()->filename + ".vs");
        VectorArray* vs = input_VectorArray(vs_file.c_str());
        if (vs == 0)
        {
            vs = new VectorArray(0,feasible.get_dimension());
        }
        ProjectLiftGenSet algorithm;
        algorithm.compute(*feasible, *gens, *vs, Globals::minimal);
        DEBUG_4ti2(*out << "Feasibles\n" << *vs << "\n";)
#endif
        ProjectLiftGenSet algorithm;
        algorithm.compute(*feasible, *gens, Globals::minimal);
    }
    else //if (Globals::generation == Globals::HYBRID)
    {
        HybridGenSet algorithm;
        algorithm.compute(*feasible, *gens, Globals::minimal);
    }
}

void
GeneratingSet::minimal()
{
    Markov algorithm;
    algorithm.compute(*feasible, *gens);
}

void
GeneratingSet::standardise()
{
    assert(gens != 0 && feasible != 0);
    Vector zero(feasible->get_dimension(), 0);
    for (int i = 0; i < gens->get_number(); ++i)
    {
        if ((*gens)[i] < zero) { (*gens)[i].mul(-1); }
    }
    gens->sort();
}

const VectorArray&
GeneratingSet::get_generating_set()
{
    assert(gens != 0);
    return *gens;
}

Feasible&
GeneratingSet::get_feasible()
{
    assert(feasible != 0);
    return *feasible;
}
