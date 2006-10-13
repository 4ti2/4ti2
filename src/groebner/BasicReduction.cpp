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

#include "BasicReduction.h"
#include "Statistics.h"

using namespace _4ti2_;

BasicReduction::BasicReduction()
{
}

BasicReduction::~BasicReduction()
{
}

void
BasicReduction::add(const Binomial& b)
{
    binomials.push_back(&b);
}

void
BasicReduction::remove(const Binomial& b)
{
    for (unsigned i = 0; i < binomials.size(); ++i)
    {
        if (binomials[i] == &b)
        {
            binomials.erase(binomials.begin()+i);
            return;
        }
    }
}

void
BasicReduction::clear()
{
    binomials.clear();
}

// Returns a point in the point set which dominates p excluding the point
// itself and b1. Returns 0 if no point exists.
const Binomial*
BasicReduction::reducable(const Binomial& b, const Binomial* b1) const
{
    //Statistics::incr_num_reducable_checks();
    for (unsigned int i = 0; i < binomials.size(); ++i)
    {
        if (Binomial::reduces(*binomials[i], b))
        {
            if (binomials[i] != &b && binomials[i] != b1) return binomials[i];
        }
    }
    return 0;
}


