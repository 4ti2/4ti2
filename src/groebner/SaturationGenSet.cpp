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
#include "SaturationGenSet.h"
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

//#define DEBUG_4ti2(X) X
#include "Debug.h"

using namespace _4ti2_;

SaturationGenSet::SaturationGenSet()
{
}

SaturationGenSet::~SaturationGenSet()
{
}

void
SaturationGenSet::compute(
                Feasible& feasible,
                VectorArray& gens,
                BitSet& sat,
                bool minimal)
{
    *out << "Computing generating set (Saturation) ...\n";
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
SaturationGenSet::compute_bounded(
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

    *out << "Saturating " << urs.get_size()-urs.count() << " variable(s).\n";
    int index = -1;
    // We look for any columns that are all zeros and set them as saturated.
    saturate_zero_columns(gens, sat, urs);
    // We see if we can infer that any variables are saturated.
    saturate(gens, sat, urs);
    // We perform one iteration of the saturation algorithm.
    char buffer[250];
    if (!is_saturated(sat, urs) && gens.get_number() != 0)
    {
        DEBUG_4ti2(*out << "Saturation:\n" << sat << "\n";)
        index = next_saturation(gens, sat, urs);
        VectorArray cost(1,dim,0);
        cost[0][index] = 0;
        sprintf(buffer, "  Sat %3d: Col: %3d ",
                        urs.get_size()-urs.count()-sat.count(), index);
        Globals::context = buffer;
        cost[0][index] = -1;

        Completion algorithm;
        algorithm.compute(feasible, cost, sat, gens);
        //DEBUG_4ti2(*out << "Gens:\n" << gens << "\n";)

        sat.set(index);
        // We look for any columns that are all zeros and set them as saturated.
        // This is possible here if truncation is strong.
        saturate_zero_columns(gens, sat, urs);
        // We see if we can infer that any variables are saturated.
        saturate(gens, sat, urs);
    }

    // We compute the saturations that we need to perform using the set of
    // vectors given from the first saturation step.
    // We store the vectors that we used to infer saturation in the vector array
    // useful_gens.
    // TODO: We could probably do this in a cleaner way.
    VectorArray useful_gens(0,gens.get_size());
    compute_saturations(gens, sat, urs, useful_gens);

    while (!is_saturated(sat, urs) && gens.get_number() != 0)
    {
        DEBUG_4ti2(*out << "Saturation:\n" << sat << "\n";)
        index = next_saturation(useful_gens, sat, urs);
        VectorArray cost(1,dim,0);
        cost[0][index] = 0;
        sprintf(buffer, "  Sat %3d: Col: %3d ",
                        urs.get_size()-urs.count()-sat.count(), index);
        Globals::context = buffer;
        cost[0][index] = -1;

        Completion algorithm;
        algorithm.compute(feasible, cost, sat, gens);
        //DEBUG_4ti2(*out << "Gens:\n" << gens << "\n";)

        sat.set(index);
        saturate_zero_columns(gens, sat, urs);
        saturate(useful_gens, sat, urs);
    }

    Globals::context = "";
    *out << "Done. ";
    *out << "Size: " << std::setw(6) << gens.get_number();
    *out << ", Time: " << t.get_elapsed_time() << " / ";
    *out << Timer::global << " secs" << std::endl;
    DEBUG_4ti2(*out << "GENS:\n" << gens << "\n";)

    if (minimal)
    {
        Markov markov;
        // TODO: We should be able to use the fast markov here, but there are
        // problems with determining what the term order should be since
        // components are permutated under projection.
        if (0) //(index != -1)
        {
            VectorArray cost(1,dim,0);
            cost[0][index] = -1;
            markov.compute(feasible, cost, gens);
        }
        else
        {
            markov.compute(feasible, gens);
        }
    }
}

int
SaturationGenSet::compute_saturations(
                    const VectorArray& gens,
                    const BitSet& sat,
                    const BitSet& urs,
                    VectorArray& useful_gens)
{
    BitSet tmp_sat(sat);
    int index;
    int count = 0;
    while (!is_saturated(tmp_sat, urs))
    {
        DEBUG_4ti2(*out << "Saturation:\n" << tmp_sat << "\n";)
        index = next_saturation(gens, tmp_sat, urs);
        ++count;
        tmp_sat.set(index);
        saturate(gens, tmp_sat, urs, useful_gens);
    }
    DEBUG_4ti2(*out << "Useful:\n" << useful_gens << "\n";)
    return count;
}

void
SaturationGenSet::saturate_zero_columns(
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
    if (num_zero_cols != 0)
    {
        *out << "  Saturated already on " << num_zero_cols << " variable(s).";
        *out << std::endl;
    }
}

bool
SaturationGenSet::is_saturated(
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
SaturationGenSet::saturate(
                const VectorArray& gens,
                BitSet& sat,
                const BitSet& urs,
                VectorArray& useful)
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
                useful.insert(gens[i]);
                changed = true;
            }
        }
    }
    return num_sats;
}

int
SaturationGenSet::saturate(
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
    if (num_sats != 0)
    {
        *out << "  Saturated already on " << num_sats << " variable(s).";
        *out << std::endl;
    }
    return num_sats;
}

int
SaturationGenSet::next_saturation(
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
SaturationGenSet::add_support(
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
SaturationGenSet::support_count(
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
SaturationGenSet::is_column_zero(
                const VectorArray& gens,
                int col)
{
    for (Index i = 0; i < gens.get_number(); ++i)
    {
        if (gens[i][col] != 0) { return false; }
    }
    return true;
}
