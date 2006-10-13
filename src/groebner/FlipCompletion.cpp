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

#include "FlipCompletion.h"
#include "Globals.h"

#include <iostream>
#include <iomanip>

//#define DEBUG_4ti2(X) X
#include "Debug.h"

using namespace _4ti2_;

FlipCompletion::FlipCompletion()
{
}

FlipCompletion::~FlipCompletion()
{
}

bool
FlipCompletion::algorithm(BinomialSet& bs, const Binomial& b)
{
    long int num_iterations = 0;
    Binomial tmp;
    Index j = 0;
    BitSet neg_supp(b.negative_support());
    BitSet pos_supp(b.positive_support());
    bool zero = false;
    while(j < bs.get_number())
    {
        const Binomial& bj = bs[j];
        //Statistics::incr_num_critical_pairs();
        //if (Binomial::is_negative_disjoint(b, bj))
        if (bs.is_negative_disjoint(j, neg_supp))
        {
            //Statistics::incr_num_graded_critical_pairs();
            //if (!Binomial::is_positive_disjoint(b, bj))
            if (!bs.is_positive_disjoint(j, pos_supp))
            {
                Binomial::spair(bj, b, tmp);
                if (!Binomial::overweight(tmp) && !bs.reducable(tmp))
                {
                    bs.reduce_negative(tmp, zero);
                    if (!zero && !Binomial::truncated(tmp)) { bs.add(tmp); }
                }
            }
        }
        ++num_iterations;
        ++j;
    }
    return true;
}
