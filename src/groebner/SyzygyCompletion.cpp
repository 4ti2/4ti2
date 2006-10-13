/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.

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

#include "SyzygyCompletion.h"
#include "SyzygyGeneration.h"
#include "BasicGeneration.h"
#include "Globals.h"
#include "Debug.h"

#include <iostream>
#include <iomanip>

using namespace _4ti2_;

SyzygyCompletion::SyzygyCompletion()
{
    name = "(U)";
    delete gen;
    gen = new SyzygyGeneration;
}

SyzygyCompletion::~SyzygyCompletion()
{
}

bool
SyzygyCompletion::algorithm(BinomialSet& bs)
{
    WeightedBinomialSet s;
    bs.auto_reduce_once();

    // long int num_iterations = 0;
    int previous_num = 0;
    int current_num = bs.get_number();

    Binomial b;

    // A heuristic to switch between using a queue or not.
    const int MAX_DIFF=200;

    while (previous_num != current_num)
    {
        *out << "\r" << Globals::context << name;
        *out << " Size: " << std::setw(8) << bs.get_number();
        *out << ", ToDo: " << std::setw(8);
        *out << current_num - previous_num << std::flush;
        DEBUG_4ti2(*out << "\n";)

        // If there are many binomials to process, then we queue them for
        // processing otherwise we just place them straight into the set of
        // binomials.
        if (current_num - previous_num >= MAX_DIFF)
        {
            gen->generate(bs, previous_num, current_num, s);
            while (!s.empty())
            {
                s.next(b);
                bool zero = false;
                bs.reduce(b, zero);
                if (!zero) { bs.add(b); }
            }
        }
        else
        {
            gen->generate(bs, previous_num, current_num, bs);
        }
        bs.auto_reduce(current_num);

        previous_num = current_num;
        current_num = bs.get_number();
    }
    bs.minimal();
    bs.reduced();
    return true;
}
