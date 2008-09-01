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

#include "groebner/ProjectLiftGenSet.h"
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
#include "groebner/BitSetStream.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/FeasibleStream.h"

//#define DEBUG_4ti2(X) X
#include "groebner/Debug.h"

using namespace _4ti2_;

ProjectLiftGenSet::ProjectLiftGenSet()
{
}

ProjectLiftGenSet::~ProjectLiftGenSet()
{
}

void
ProjectLiftGenSet::compute(
                Feasible& feasible,
                VectorArray& gens,
                VectorArray& feasibles,
                bool minimal)
{
    *out << "Computing generating set (Project-and-Lift) ...\n";
    compute_unbounded(feasible, gens, feasibles, minimal);
}


void
ProjectLiftGenSet::compute_unbounded(
                Feasible& feasible,
                VectorArray& gens,
                VectorArray& feasibles,
                bool minimal)
{
    DEBUG_4ti2(*out << "Compute Unbounded.\n";)
    DEBUG_4ti2(*out << feasible << "\n";)

    if (!feasible.get_bnd().empty())
    {
        BitSet proj(feasible.get_urs());
        proj.set_union(feasible.get_unbnd());
        Feasible bounded(feasible, proj);
        compute_bounded(bounded, gens, feasibles, minimal);
    }

    if (!feasible.get_unbnd().empty())
    {
        VectorArray basis(feasible.get_basis());
        int row = upper_triangle(basis, feasible.get_bnd());
        basis.remove(0,row);
        gens.insert(basis);
        gens.insert(feasible.get_ray());
        make_feasible(feasibles, feasible.get_ray());
        *out << "  Lifting " << feasible.get_unbnd().count() << " unbounded.\n";
    }
}

void
ProjectLiftGenSet::compute_bounded(
                Feasible& feasible,
                VectorArray& gens,
                VectorArray& feasibles,
                bool minimal)
{
    DEBUG_4ti2(*out << "Compute Bounded.\n";)
    DEBUG_4ti2(*out << feasible << "\n";)
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

    BitSet fin_union_urs(fin.get_size());
    BitSet::set_union(fin, urs, fin_union_urs);
    int index = 0;
    while (index < dim && fin_union_urs[index]) { ++index; }
    assert(index != dim);
    DEBUG_4ti2(*out<<"Column = "<<index<<"\n";)
    fin_union_urs.set(index);
    Feasible projection(feasible, fin_union_urs);
    compute_unbounded(projection, gens, feasibles, false);
    VectorArray cost(1,dim,0);
    cost[0][index] = -1;
    char buffer[250];
    sprintf(buffer, "  Lift %3d: Col: %3d ", fin.count()+1, index); 
    Globals::context = buffer;
    Completion algorithm;
    DEBUG_4ti2(*out << "Column:" << index << "\n";)
    algorithm.compute(projection, cost, gens, feasibles);
    DEBUG_4ti2(*out << "GENS:\n" << gens << "\n";)

    Timer t;
    add_support(gens, fin);
    int column = index;
    while (!fin.empty())
    {
        column = next_support(gens, fin);
        VectorArray cost(1,dim,0);
        cost[0][column] = -1;
        sprintf(buffer, "  Lift %3d: Col: %3d ", fin.count(), column); 
        Globals::context = buffer;
        BitSet fin_union_urs(fin.get_size());
        BitSet::set_union(fin, urs, fin_union_urs);
        Feasible projection(feasible, fin_union_urs);
        DEBUG_4ti2(*out << projection << "\n";)

        Completion algorithm;
        algorithm.compute(projection, cost, gens, feasibles);
        DEBUG_4ti2(*out << "GENS:\n" << gens << "\n";)

        fin.unset(column);
        add_support(gens, fin);
    }
    Globals::context = "";
    *out << "Done. ";
    *out << "Size: " << std::setw(6) << gens.get_number();
    *out << ", Time: " << t << " / ";
    *out << Timer::global << " secs" << std::endl;

    if (minimal)
    {
        Markov markov;
        VectorArray cost(1,dim,0);
        cost[0][column] = -1;
        markov.compute(feasible, cost, gens);
        //markov.compute(feasible, gens);
    }
}

int
ProjectLiftGenSet::add_support(
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
        *out << "  Lifted already on " << num_lifts << " variable(s).";
        *out << std::endl;
    }
    return num_lifts;
}

int
ProjectLiftGenSet::positive_count(
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
ProjectLiftGenSet::next_support(
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

void
ProjectLiftGenSet::make_feasible(
                VectorArray& feasibles,
                const Vector& ray)
{
    DEBUG_4ti2(*out << "Making feasible.\n";)
    DEBUG_4ti2(*out << "Feasibles:\n" << feasibles << "\n";)
    DEBUG_4ti2(*out << "Ray:\n" << ray << "\n";)
    IntegerType factor = 0;
    for (int i = 0; i < feasibles.get_number(); ++i)
    {
        for (int j = 0; j < ray.get_size(); ++j)
        {
            if (feasibles[i][j] < 0 && ray[j] > 0)
            {
                IntegerType ratio = -feasibles[i][j]/ray[j] + 1;
                if (ratio > factor) { factor = ratio; }
            }
        }
        if (factor != 0) { feasibles[i].add(ray, factor); } 
    }
}
