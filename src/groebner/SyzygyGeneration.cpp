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

#include "SyzygyGeneration.h"
#include "Statistics.h"

#include <iostream>
#include <iomanip>

// TODO: Should we eliminate overweight binomials before performing syzygy
// generation.

using namespace _4ti2_;

void
SyzygyGeneration::generate(
                const BinomialSet& bs,
                Index i,
                BinomialCollection& bc)
{
    Binomial b;
    const Binomial& bi = bs[i];

    std::vector<std::pair<IntegerType,int> > ordering;
    ordering.reserve(i); // Reserve space.
    for (Index j = 0; j < i; ++j)
    {
        //Statistics::incr_num_critical_pairs();
        //if (!Binomial::is_positive_disjoint(bi, bs[j]))
        if (!bs.is_positive_disjoint(i, j))
        {
            //Statistics::incr_num_disjoint_critical_pairs();
            IntegerType normij = Binomial::max_l1_norm(bi,bs[j]);
            ordering.push_back(std::pair<IntegerType,int>(normij,j));
        }
    }
    sort(ordering.begin(), ordering.end());
    std::vector<int> syzergies;
    for (Index k = 0; k < (Index) ordering.size(); ++k)
    {
        Index index = ordering[k].second;
        const Binomial& bk = bs[index];
        if (!dominated(syzergies, bs, bi, bk))
        {
            syzergies.push_back(index);
            //Statistics::incr_num_syzygy_critical_pairs();
            //if (Binomial::is_negative_disjoint(bi, bk))
            if (bs.is_negative_disjoint(i,index))
            {
                //Statistics::incr_num_graded_critical_pairs();
                Binomial::spair(bi, bk, b);
                if (BinomialSet::check(bs, b)) { bc.add(b); }
            }
        }
    }
}

void
SyzygyGeneration::generate(
                const BinomialSet& bs,
                Index start,
                Index end,
                BinomialCollection& bc)
{
    Binomial b;
    std::vector<std::vector<int> > back_syzergies(end-start);
    for (int i = end-1; i >= (int) start; --i)
    {
        const Binomial& bi = bs[i];
        std::vector<std::pair<IntegerType,int> > ordering;
        ordering.reserve(i);
        for (int j = 0; j < i; ++j)
        {
            //Statistics::incr_num_critical_pairs();
            //if (!Binomial::is_positive_disjoint(bi, bs[j]))
            if (!bs.is_positive_disjoint(i, j))
            {
                //Statistics::incr_num_disjoint_critical_pairs();
                IntegerType norm = Binomial::max_l1_norm(bi, bs[j]);
                ordering.push_back(std::pair<IntegerType,int>(norm,j));
            }
        }
        sort(ordering.begin(),ordering.end());
        std::vector<int> syzergies;
        for (Index k = 0; k < (Index) ordering.size(); ++k)
        {
            Index index = ordering[k].second;
            const Binomial& bk = bs[index];
            if (!dominated(syzergies, back_syzergies[i-start], bs, bi, bk))
            {
                syzergies.push_back(index);
                //Statistics::incr_num_syzygy_critical_pairs();
                //if (Binomial::is_negative_disjoint(bi, bk))
                if (bs.is_negative_disjoint(i,index))
                {
                    Binomial::spair(bi, bk, b);
                    if (BinomialSet::check(bs, b)) { bc.add(b); }
                }
            }
        }
        for (Index k = 0; k < (Index) syzergies.size(); ++k)
        {
            int index = syzergies[k]-start;
            if (index >= 0) back_syzergies[index].push_back(i);
        }
        back_syzergies[i-start].clear();
    }
}

bool
SyzygyGeneration::dominated(
                std::vector<int>& syzergies,
                const BinomialSet& bs,
                const Binomial& b0,
                const Binomial& b1)
{
    for (Index i = 0; i < (Index) syzergies.size(); ++i)
    {
         const Binomial& b2 = bs[syzergies[i]];
         if (Binomial::reduces(b0, b2, b1)) return true;
    }
    return false;
}

bool
SyzygyGeneration::dominated(
                std::vector<int>& syzergies,
                std::vector<int>& back_syzergies,
                const BinomialSet& bs,
                const Binomial& b0,
                const Binomial& b1)
{
    for (Index i = 0; i < (Index) syzergies.size(); ++i)
    {
         const Binomial& b2 = bs[syzergies[i]];
         if (Binomial::reduces(b0, b2, b1)) { return true; }
    }

    for (Index i = 0; i < (Index) back_syzergies.size(); ++i)
    {
         const Binomial& b2 = bs[back_syzergies[i]];
         if (Binomial::reduces(b0, b2, b1))
         {
             IntegerType norm01 = Binomial::max_l1_norm(b0, b1);
             IntegerType norm02 = Binomial::max_l1_norm(b0, b2);
             if (norm01 != norm02)
             {
                IntegerType norm12 = Binomial::max_l1_norm(b1, b2);
                if (norm01 != norm12) { return true; }
             }
         }
    }
    return false;
}
