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

#include "groebner/HybridGenSet.h"
#include "groebner/SaturationGenSet.h"
#include "groebner/Statistics.h"
#include "groebner/Bounded.h"
#include "groebner/Globals.h"
#include "groebner/BasicCompletion.h"
#include "groebner/SyzygyCompletion.h"
#include "groebner/HermiteAlgorithm.h"
#include "groebner/Completion.h"
#include "groebner/Timer.h"
#include "groebner/Markov.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

#include "groebner/BitSetStream.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/FeasibleStream.h"

//#define DEBUG_4ti2(X) X
#include "Debug.h"

using namespace _4ti2_;

HybridGenSet::HybridGenSet()
{
}

HybridGenSet::~HybridGenSet()
{
}

void
HybridGenSet::compute(
                Feasible& feasible,
                VectorArray& gens,
                bool minimal)
{
    *out << "Computing generating set (Hybrid) ...\n";
    DEBUG_4ti2(*out << "Compute Unbounded.\n";)
    DEBUG_4ti2(*out << feasible << "\n";)

    if (!feasible.get_bnd().empty())
    {
        BitSet proj(feasible.get_urs());
        proj.set_union(feasible.get_unbnd());
        Feasible bounded(feasible, proj);
        compute_bounded(bounded, gens, minimal);
    }

    if (!feasible.get_unbnd().empty())
    {
        VectorArray basis(feasible.get_basis());
        int row = upper_triangle(basis, feasible.get_bnd());
        basis.remove(0,row);
        gens.insert(basis);
        gens.insert(feasible.get_ray());
    }
}

void
HybridGenSet::compute_bounded(
                Feasible& feasible,
                VectorArray& gens,
                bool minimal)
{
    int dim = feasible.get_dimension();
    const BitSet& urs = feasible.get_urs();
    const BitSet& unbnd = feasible.get_unbnd();

    if (!unbnd.empty())
    {
        std::cerr << "ERROR: Expected fully bounded problem.\n";
        exit(1);
    }

    BitSet fin(dim);
    Vector rhs(dim,1);
    if (feasible.get_rhs() != 0) { rhs = *feasible.get_rhs(); }
    bounded_projection(feasible.get_matrix(), feasible.get_basis(), urs, rhs, fin);

    // Use the saturation algorithm to compute the first generating set.
    assert(BitSet::set_disjoint(fin, urs));
    BitSet fin_union_urs(fin.get_size());
    BitSet::set_union(fin, urs, fin_union_urs);
    *out << "Phase 1:\n";
    Feasible sat_feasible(feasible, fin_union_urs);
    SaturationGenSet saturation_algorithm;
    BitSet sat(feasible.get_dimension());
    saturation_algorithm.compute(sat_feasible, gens, sat, false);

    Timer t;
    *out << "Phase 2:\n";
    *out << "Lifting " << fin.count() << " variable(s).\n";
    add_support(gens, fin);
    char buffer[250];
    int column = -1;
    while (!fin.empty())
    {
        column = next_support(gens, fin);
        VectorArray cost(1,dim,0);
        cost[0][column] = -1;
        sprintf(buffer, "  Lift %3d: Col: %3d ", fin.count(), column); 
        Globals::context = buffer;
        BitSet::set_union(fin, urs, fin_union_urs);
        Feasible projection(feasible, fin_union_urs);

        Completion algorithm;
        algorithm.compute(projection, cost, gens);

        fin.unset(column);
        add_support(gens, fin);
    }
    Globals::context = "";
    *out << "Done. ";
    *out << "Size: " << std::setw(6) << gens.get_number();
    *out << ", Time: " << t.get_elapsed_time() << " / ";
    *out << Timer::global << " secs" << std::endl;

    if (minimal)
    {
        Markov markov;
        if (column != -1)
        {
            VectorArray cost(1,dim,0);
            cost[0][column] = -1;
            markov.compute(feasible, cost, gens);
        }
        else
        {
            markov.compute(feasible, gens);
        }
    }
}

int
HybridGenSet::add_support(
                const VectorArray& gens,
                BitSet& fin)
{
    int num_lifts = 0;
    for (Index c = 0; c < gens.get_size(); ++c)
    {
        if (fin[c] == 1 && positive_count(gens, c) == 0)
        {
            fin.unset(c);
            ++num_lifts;
        }
    }
    if (num_lifts != 0) 
    {
        *out << "  Lifted already on " << num_lifts << " variable(s)";
        *out << std::endl;
    }
    return num_lifts;
}

int
HybridGenSet::positive_count(
                const VectorArray& gens,
                int c)
{
    int count = 0;
    for (Index r = 0; r < gens.get_number(); ++r)
    {
        if (gens[r][c] > 0) ++count;
    }
    return count;
}

int
HybridGenSet::next_support(
                const VectorArray& gens,
                BitSet& fin)
{
    int min_count = gens.get_number()+1;
    int col = -1;
    for (Index c = 0; c < gens.get_size(); ++c)
    {
        if (fin[c] == 1)
        {
            int count = positive_count(gens, c);
            if (count < min_count)
            {
                col = c;
                min_count = count;
            }
        }
    }
    assert(col != -1);
    assert(min_count != (int) gens.get_number()+1);
    return col;
}
