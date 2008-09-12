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

#include "Feasible.h"
#include "FeasibleStream.h"
#include "MaxMinGenSet.h"
#include "Completion.h"
#include "SyzygyCompletion.h"
#include "OrderedCompletion.h"
#include "Timer.h"
#include "Bounded.h"
#include "HermiteAlgorithm.h"
#include "VectorArrayStream.h"
#include "Globals.h"
#include "Markov.h"

#include <iostream>
#include <iomanip>
#include "BitSetStream.h"

#define DEBUG_4ti2(X) //X
#include "Debug.h"

using namespace _4ti2_;

MaxMinGenSet::MaxMinGenSet()
{
}

MaxMinGenSet::~MaxMinGenSet()
{
}

void
MaxMinGenSet::compute(
                Feasible& feasible,
                VectorArray& gens,
                BitSet& sat,
                bool minimal)
{
    *out << "Computing generating set (MaxMin) ...\n";
    DEBUG_4ti2(*out << "Compute Unbounded.\n";)
    DEBUG_4ti2(*out << feasible << "\n";)

    if (!feasible.get_bnd().empty())
    {
        BitSet proj(feasible.get_urs());
        proj.set_union(feasible.get_unbnd());
        Feasible bounded(feasible, proj);
        compute_bounded(bounded, gens, sat, minimal);
    }

    if (!feasible.get_unbnd().empty())
    {
        VectorArray basis(feasible.get_basis());
        int row = hermite(basis, feasible.get_bnd());
        basis.remove(0,row);
        gens.insert(basis);
        gens.insert(feasible.get_ray());
    }
}

void
MaxMinGenSet::compute_bounded(
                Feasible& feasible,
                VectorArray& gens,
                BitSet& sat,
                bool minimal)
{
    assert(sat.get_size() == feasible.get_dimension());
    if (!feasible.get_unbnd().empty())
    {
        std::cerr << "Attempting saturation when not fully bounded.\n";
        exit(1);
    }

    int dim = feasible.get_dimension();
    const BitSet& urs = feasible.get_urs();

    DEBUG_4ti2(*out << feasible;)
    Timer t;

    gens.insert(feasible.get_basis());

    BitSet todo(dim, false);
    compute_saturations(gens, sat, urs, todo);
    DEBUG_4ti2(*out << "Saturations ToDo:\n" << todo << "\n";)

    // Extend the gens.
    VectorArray ext_gens(gens.get_number(), gens.get_size()+1, 0);
    VectorArray::lift(gens, 0, gens.get_size(), ext_gens);
    Vector ext_vec(gens.get_size()+1,0);
    for (int i = 0; i < gens.get_size(); ++i)
    {  if (todo[i]) { ext_vec[i] = 1; } }
    ext_vec[gens.get_size()] = -1;
    //ext_vec[gens.get_size()] = 1;
    ext_gens.insert(ext_vec, 0);
    DEBUG_4ti2(*out << "Extended Gens:\n" << ext_gens << "\n";)

    // Extend the matrix.
    const VectorArray& matrix = feasible.get_matrix();
    VectorArray ext_matrix(matrix.get_number(),matrix.get_size()+1,0);
    VectorArray::lift(matrix, 0, matrix.get_size(), ext_matrix);
    for (int i = 0; i < matrix.get_number(); ++i)
    {
        IntegerType acc = 0;
        for (int j = 0; j < matrix.get_size(); ++j)
        { if (todo[j]) { acc += matrix[i][j]; } }
        ext_matrix[i][matrix.get_size()] = acc;
    }
    DEBUG_4ti2(*out << "Extended Matrix:\n" << ext_matrix << "\n";)

    // Extend the urs.
    BitSet ext_urs(urs.get_size()+1, false);
    BitSet::extend(urs, ext_urs);
    //ext_urs.set(urs.get_size());
    DEBUG_4ti2(*out << "Extended URS:\n" << ext_urs << "\n";)

    // Extend the sat.
    BitSet ext_sat(sat.get_size()+1, false);
    BitSet::extend(sat, ext_sat);
    DEBUG_4ti2(*out << "Extended SAT:\n" << ext_sat << "\n";)

    // Extend the feasibles.
    Feasible ext_feasible(&ext_gens, &ext_matrix, &ext_urs);
    DEBUG_4ti2(*out << "Extended Feasible:\n" << ext_feasible << "\n";)

    // We perform one iteration of the saturation algorithm.
    DEBUG_4ti2(*out << "MaxMin:\n" << todo << "\n";)
    VectorArray ext_cost(1,dim+1,0);
    ext_cost[0][dim] = -1;
    //ext_cost[0][dim] = 1;
    for (int i = 0; i < dim; ++i)
    {
        if (todo[i])
        {
            Vector new_cost(dim+1,0);
            new_cost[i] = -1;
            ext_cost.insert(new_cost);
        }
    }
    DEBUG_4ti2(*out << "Ext Cost:\n" << ext_cost << "\n";)

    Completion algorithm;
    algorithm.compute(ext_feasible, ext_cost, ext_sat, ext_gens);
    DEBUG_4ti2(*out << "Ext Gens:\n" << ext_gens << "\n";)

    std::cout << "Ext vector:\n" << ext_vec << "\n";
    for (int i = ext_gens.get_number()-1; i >= 0 ; --i)
    {
        if (ext_gens[i][dim] != 0)
        {
            ext_gens[i].add(ext_vec, ext_gens[i][dim]);
            if (ext_gens[i].is_zero()) { ext_gens.remove(i); }
        }
    }
    DEBUG_4ti2(*out << "Ext Gens Modified:\n" << ext_gens << "\n";)

    gens.renumber(ext_gens.get_number());
    VectorArray::project(ext_gens, 0, gens.get_size(), gens);

    Globals::context = "";
    *out << "Done. ";
    *out << "Size: " << std::setw(6) << gens.get_number();
    *out << ", Time: " << t.get_elapsed_time() << " / ";
    *out << Timer::global << " secs" << std::endl;
    DEBUG_4ti2(*out << "GENS:\n" << gens << "\n";)

    if (minimal)
    {
        Markov markov;
        markov.compute(feasible, gens);
    }
}

int
MaxMinGenSet::compute_saturations(
                    const VectorArray& gens,
                    const BitSet& sat,
                    const BitSet& urs,
                    BitSet& todo)
{
    BitSet tmp_sat(sat);
    DEBUG_4ti2(*out << "Saturation:\n" << tmp_sat << "\n";)

    // We look for any columns that are all zeros and set them as saturated.
    saturate_zero_columns(gens, tmp_sat, urs);

    int index;
    int count = 0;
    while (!is_saturated(tmp_sat, urs))
    {
        index = next_saturation(gens, tmp_sat, urs);
        ++count;
        tmp_sat.set(index);
        todo.set(index);
        DEBUG_4ti2(*out << "Saturation:\n" << tmp_sat << "\n";)
        saturate(gens, tmp_sat, urs);
    }
    DEBUG_4ti2(*out << "ToDo:\n" << todo << "\n";)
    return count;
}

void
MaxMinGenSet::saturate_zero_columns(
                const VectorArray& gens,
                BitSet& sat,
                const BitSet& urs)
{
    int num_zero_cols = 0;
    for (Index c = 0; c < gens.get_size(); ++c)
    {
        if (urs[c] == 0 && sat[c] == 0 && is_column_zero(gens, c))
        {
            sat.set(c);
            ++num_zero_cols;
        }
    }
}

bool
MaxMinGenSet::is_saturated(
                const BitSet& sat,
                const BitSet& urs)
{
    assert(sat.get_size() == urs.get_size());
    for (Index i = 0; i < sat.get_size(); ++i)
    {
        if (sat[i] == 0 && urs[i] == 0) { return false; }
    }
    return true;
}

int
MaxMinGenSet::saturate(
                const VectorArray& gens,
                BitSet& sat,
                const BitSet& urs)
{
    assert(gens.get_size() == sat.get_size());
    assert(sat.get_size() == urs.get_size());
    bool changed = true;
    int num_sats = 0;
    while (changed)
    {
        changed = false;
        for (Index i = 0; i < gens.get_number(); ++i)
        {
            int pos, neg;
            support_count(gens[i], sat, urs, pos, neg);
            if ((pos == 0 && neg != 0) || (neg == 0 && pos != 0))
            {
                DEBUG_4ti2(*out << "Adding support:\n" << gens[i] << "\n";)
                num_sats += add_support(gens[i], sat, urs);
                changed = true;
            }
        }
    }
    return num_sats;
}

int
MaxMinGenSet::next_saturation(
                const VectorArray& gens,
                BitSet& sat, 
                const BitSet& urs)
{
    assert(gens.get_size() == sat.get_size());
    assert(sat.get_size() == urs.get_size());
    assert(!is_saturated(sat, urs));
    int min = gens.get_size();
    int index = -1;
    int pos_neg = 0;
    for (Index i = 0; i < gens.get_number(); ++i)
    {
        int pos, neg;
        support_count(gens[i], sat, urs, pos, neg);
        if (pos != 0 && pos < min)
        {
            min = pos;
            index = i;
            pos_neg = 1;
        }
        if (neg != 0 && neg < min)
        {
            min = neg;
            index = i;
            pos_neg = -1;
        }
    }
    assert(min != (int) gens.get_size());
    assert(index != -1);
    assert(pos_neg != 0);

    for (Index i = 0; i < gens.get_size(); ++i)
    {
        if (sat[i] == 0 && urs[i] == 0 && gens[index][i]*pos_neg > 0)
        {
            return i;
        }
    }
    assert(false);
    return 0;
}

int
MaxMinGenSet::add_support(
                const Vector& p,
                BitSet& sat,
                const BitSet& urs)
{
    assert(p.get_size() == sat.get_size());
    assert(sat.get_size() == urs.get_size());
    int num_sats = 0;
    for (Index i = 0; i < p.get_size(); ++i)
    {
        if (sat[i] == 0 && urs[i] == 0 && p[i] != 0)
        {
            sat.set(i);
            ++num_sats;
        }
    }
    return num_sats;
}

void
MaxMinGenSet::support_count(
                const Vector& p,
                BitSet& sat,
                const BitSet& urs,
                int& pos,
                int& neg)
{
    assert(p.get_size() == sat.get_size());
    assert(sat.get_size() == urs.get_size());
    pos = 0;
    neg = 0;
    for (Index i = 0; i < p.get_size(); ++i)
    {
        if (sat[i] == 0 && urs[i] == 0)
        {
            if (p[i] > 0) { ++pos; }
            if (p[i] < 0) { ++neg; }
        }
    }
}

bool
MaxMinGenSet::is_column_zero(
                const VectorArray& gens,
                int col)
{
    for (Index i = 0; i < gens.get_number(); ++i)
    {
        if (gens[i][col] != 0) { return false; }
    }
    return true;
}
