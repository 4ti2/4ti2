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

#include "Minimize.h"
#include "BinomialFactory.h"
#include "BinomialSet.h"
#include "Bounded.h"

using namespace _4ti2_;

Minimize::Minimize()
{
}

Minimize::~Minimize()
{
}

// Compute the normal form of point.
void
Minimize::minimize(
                Feasible& feasible,
                const VectorArray& cost,
                const VectorArray& vs,
                Vector& v)
{
    BinomialFactory factory(feasible, cost);
    BinomialSet bs;
    factory.convert(vs, bs);
    Binomial b;
    factory.convert(v, b);
    bs.minimize(b);
    factory.convert(b, v);
    bs.clear();
}

// Compute the normal form of a set of points.
void
Minimize::minimize(
                Feasible& feasible,
                const VectorArray& cost,
                const VectorArray& vs,
                VectorArray& ps)
{
    BinomialFactory factory(feasible, cost);
    BinomialSet bs;
    factory.convert(vs, bs);
    Binomial b;
    for (int i = 0; i < ps.get_number(); ++i)
    {
        factory.convert(ps[i], b);
        bs.minimize(b);
        factory.convert(b, ps[i]);
    }
    bs.clear();
}
