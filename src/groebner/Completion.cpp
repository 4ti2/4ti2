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

#include "Completion.h"
#include "BasicCompletion.h"
#include "SyzygyCompletion.h"
#include "FlipCompletion.h"
#include "OrderedCompletion.h"
#include "Algorithm.h"
#include "Feasible.h"
#include "FeasibleStream.h"
#include "BinomialFactory.h"
#include "BinomialSetStream.h"

#include <iostream>
#include <iomanip>
#include "VectorStream.h"
#include "VectorArrayStream.h"
#include "BitSetStream.h"
#include "Globals.h"

//#define DEBUG_4ti2(X) X
#include "Debug.h"

using namespace _4ti2_;

Completion::Completion()
{
    switch (Globals::algorithm) {
    case Globals::WEIGHTED:
        alg = new OrderedCompletion;
        break;
    case Globals::GEBAUER_MOELLER:
        alg = new SyzygyCompletion;
        break;
    case Globals::NORMAL:
        alg = new BasicCompletion;
        break;
    // In the default case, the algorithm is chosen according to the properties
    // of the input.
    default: // Globals::DEFAULT
        alg = 0;
        break;
    }
}

Completion::~Completion()
{
    delete alg;
}

void
Completion::compute(
                Feasible& feasible,
                const VectorArray& cost,
                VectorArray& vs,
                VectorArray& feasibles)
{
    // Reset the stopwatch.
    t.reset();

    DEBUG_4ti2(
        *out << "COST:\n" << cost << "\n";
        *out << "Feasible:\n" << feasible << "\n";
        *out << "VS:\n" << vs << "\n";
        *out << std::endl;
    )

    // If no algorithm has been selected, then select one according to the ratio
    // of bounded variables to unbounded variables.
    if (alg == 0)
    {
        if (feasible.get_unbnd().count() / (feasible.get_bnd().count()+1) >= 2)
        {
            //alg = new OrderedCompletion;
            alg = new SyzygyCompletion;
        }
        else
        {
            alg = new BasicCompletion;
        }
    }

    // Initialise the binomials.
    BinomialFactory factory(feasible, cost);
    BinomialSet bs;
    factory.convert(vs, bs);
    DEBUG_4ti2(*out << "Before:\n" << bs;)
    alg->algorithm(bs);
    DEBUG_4ti2(*out << "After:\n" << bs;)
    Binomial b;
    for (int i = 0; i < feasibles.get_number(); ++i)
    {
        factory.convert(feasibles[i], b);
        bs.minimize(b);
        factory.convert(b, feasibles[i]);
    }
    factory.convert(bs, vs);
    bs.clear();

    *out << "\r" << Globals::context << alg->get_name();
    *out << " Size: " << std::setw(6) << vs.get_number();
    *out << ", Time: " << t << " / ";
    *out << Timer::global << " secs.          " << std::endl;
}

void
Completion::compute(
                Feasible& feasible,
                const VectorArray& cost,
                const BitSet& sat,
                VectorArray& vs,
                VectorArray& feasibles)
{
    // Reset the stopwatch.
    t.reset();

    DEBUG_4ti2(
        *out << "COST:\n" << cost << "\n";
        *out << "Feasible:\n" << feasible << "\n";
        *out << "VS:\n" << vs << "\n";
        *out << std::endl;
    )

    // If no algorithm has been selected, then select one according to the ratio
    // of saturated variables to unsaturated variables.
    if (alg == 0)
    {
        if ((feasible.get_dimension()-sat.count()) / (sat.count()+1) >= 3)
        {
            alg = new SyzygyCompletion;
        }
        else
        {
            alg = new BasicCompletion;
        }
    }

    // Initialise the binomials.
    BinomialFactory factory(feasible, cost, sat);
    BinomialSet bs;
    factory.convert(vs, bs);
    DEBUG_4ti2(*out << "Before:\n" << bs;)
    alg->algorithm(bs);
    DEBUG_4ti2(*out << "After:\n" << bs;)
    Binomial b;
    for (int i = 0; i < feasibles.get_number(); ++i)
    {
        factory.convert(feasibles[i], b);
        bs.minimize(b);
        factory.convert(b, feasibles[i]);
    }
    factory.convert(bs, vs);
    bs.clear();

    *out << "\r" << Globals::context << alg->get_name();
    *out << " Size: " << std::setw(6) << vs.get_number();
    *out << ", Time: " << t << " / ";
    *out << Timer::global << " secs.          " << std::endl;

    bs.clear();
}

void
Completion::set_algorithm(Algorithm* _alg)
{
    delete alg;
    alg = _alg;
}
