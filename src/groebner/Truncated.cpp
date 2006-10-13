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

#include "Truncate.h"
#include "Binomial.h"
#include "BinomialFactory.h"

using namespace _4ti2_;

void
_4ti2_::truncate(Feasible& feasible, VectorArray& vs)
{
    VectorArray cost(0, feasible.get_dimension());
    BinomialFactory factory(feasible, cost);

    Binomial b;
    for (int i = vs.get_number()-1; i >= 0; --i)
    {
        factory.convert(vs[i], b);
        if (Binomial::overweight(b) || Binomial::truncated(b))
        {
            vs.remove(i);
        }
    }
}
