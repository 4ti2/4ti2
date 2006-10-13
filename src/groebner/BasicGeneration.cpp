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

#include "BasicGeneration.h"
#include "BinomialStream.h"
#include "Globals.h"

//#define DEBUG_4ti2(X) X
#include "Debug.h"

// TODO: Use temporary variables better.
// Do not use a separate b per function call.

using namespace _4ti2_;

void
BasicGeneration::generate(
                const BinomialSet& bs,
                Index i,
                Index j,
                BinomialCollection& bc)
{
    for (Index k = i; k < j; ++k) { generate(bs, k, bc); }
}

void
output_stuff(const Binomial& b1, const Binomial& b2)
{
    Binomial z;
    for (int i = 0; i < Binomial::urs_end; ++i)
    {
        if (b1[i] >= 0 && b1[i] >= b2[i]) { z[i] = b1[i]; }
        else if (b2[i] >= 0 && b2[i] >= b1[i]) { z[i] = b2[i]; }
        else { z[i] = 0; }
    }
    Binomial x;
    for (int i = 0; i < Binomial::urs_end; ++i)
    {
        x[i] = z[i] - b1[i];
    }
    Binomial y;
    for (int i = 0; i < Binomial::urs_end; ++i)
    {
        y[i] = z[i] - b2[i];
    }
    for (int i = Binomial::urs_end; i < Binomial::size; ++i)
    { z[i] = 0; x[i] = 0; y[i] = 0; }

    std::cout << "Z = " << z << "\n";
    std::cout << "X = " << x << "\n";
    std::cout << "Y = " << y << "\n";
}

void
BasicGeneration::generate(
                const BinomialSet& bs,
                Index i,
                BinomialCollection& bc)
{
    Binomial b;
    const Binomial& bi = bs[i];
    for (Index j = 0; j < i; ++j)
    {
        DEBUG_4ti2(*out << "(" << j << "," << i << "):\n";)
        DEBUG_4ti2(*out << "u = " << bs[j] << "\n";)
        DEBUG_4ti2(*out << "v = " << bi << "\n";)
        if (bs.is_negative_disjoint(i,j))
        {
            if (!bs.is_positive_disjoint(i,j))
            {
                Binomial::spair(bi, bs[j], b);
                DEBUG_4ti2(*out << b << "\n";)
                DEBUG_4ti2(output_stuff(bs[j], bi);)
                if (BinomialSet::check(bs, b)) { bc.add(b); }
            }
            else
            {
                DEBUG_4ti2(*out << "Positive disjoint.\n";)
            }
        }
        else
        {
            DEBUG_4ti2(*out << "Negative disjoint.\n";)
        }
    }
}
