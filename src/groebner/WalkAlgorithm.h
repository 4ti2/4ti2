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

#ifndef _4ti2__WalkAlgorithm_
#define _4ti2__WalkAlgorithm_

#include "groebner/VectorArray.h"
#include "groebner/BitSet.h"
#include "groebner/Weight.h"
#include "groebner/Timer.h"
#include "groebner/Feasible.h"
#include "groebner/Binomial.h"
#include "groebner/BinomialSet.h"
#include "groebner/TermOrder.h"

namespace _4ti2_
{

class WalkAlgorithm
{
public:
    WalkAlgorithm();
    virtual ~WalkAlgorithm();

    virtual void compute(
                    Feasible& feasible,
                    const VectorArray& costold,
                    VectorArray& vs1,
                    const VectorArray& costnew);

protected:
    IntegerType compare(const Binomial& b1, const Binomial& b2);
    RationalType tvalue(const Binomial& b);
    void tvector(Vector& c1, Vector& c2, Vector& v, Vector& tv);
    bool next(const BinomialSet& bs, const TermOrder& term_order, int& min);

    int costnew_start;
    int costnew_end;
    int costold_start;
    int costold_end;
    Timer t;
};

}

#endif
