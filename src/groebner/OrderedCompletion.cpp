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

#include "OrderedCompletion.h"
#include "Generation.h"
#include "Globals.h"
#include "Debug.h"

#include <iostream>
#include <iomanip>

using namespace _4ti2_;

OrderedCompletion::OrderedCompletion()
{
    name = "(W)";
}

OrderedCompletion::~OrderedCompletion()
{
}

bool
OrderedCompletion::algorithm(BinomialSet& bs)
{
    bs.auto_reduce_once();
    WeightedBinomialSet c;
    for (Index i = 0; i < bs.get_number(); ++i) { c.add(bs[i]); }
    bs.clear();
    return algorithm(c, bs);
}

bool
OrderedCompletion::algorithm(
                WeightedBinomialSet& c,
                BinomialSet& bs)
{
    long int num_iterations = 0;
    Binomial b;
    bool homogeneous = true;
    if (Binomial::bnd_end != Binomial::rs_end) { homogeneous = false; }
    while(!c.empty())
    {
        ++num_iterations;
        c.next(b);
        bool zero = false;
        bs.reduce(b, zero);
        if (!zero)
        {
            bs.add(b);
            gen->generate(bs, bs.get_number()-1, c);
        }
        if (num_iterations % Globals::output_freq == 0)
        {
            *out << "\r" << Globals::context << name;
            *out << " Size: " << std::setw(6) << bs.get_number();
            *out << " Degree: " << std::setw(6) << c.min_grade();
            *out << " ToDo: " << std::setw(6) << c.get_size();
            *out << std::flush;
            DEBUG_4ti2(*out << "\n";)
        }
        if (!homogeneous && num_iterations % Globals::auto_reduce_freq == 0)
        {
            int i = bs.get_number();
            bs.auto_reduce_once(i);
            if (i != bs.get_number())
            {
                gen->generate(bs, i, bs.get_number()-1, c);
            }
        }
    }
    // If the lattice is not homogeneous, then minimize the set.
    // If is is homogeneous, then it must be minimal.
    if (!homogeneous) { bs.minimal(); }
    // TODO: There is no need to perform reduction completely.
    // We can do it incrementally in a smarter way.
    bs.reduced();
    return true;
}
