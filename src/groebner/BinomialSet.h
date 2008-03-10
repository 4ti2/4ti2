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

#ifndef _4ti2__BinomialSet_
#define _4ti2__BinomialSet_

#include "groebner/Binomial.h"
#include "groebner/VectorArray.h"
#include "groebner/Reduction.h"
#include "groebner/BitSet.h"
#include "groebner/IndexBinomialSet.h"
#include "groebner/BinomialCollection.h"
#include "groebner/TermOrder.h"

#include <iostream>
#include "groebner/BitSetStream.h"

namespace _4ti2_
{

class BinomialSet : public BinomialCollection
{
public:
    BinomialSet();
    ~BinomialSet();

    const Binomial& operator[](Index) const;

    Size get_number() const;

    void add(const Binomial& b);
    void remove(Index i);
    void clear();

    bool reduce(    Binomial& b,
                    bool& zero,
                    Binomial* ptr = 0) const;
    bool reduce_negative(
                    Binomial& b,
                    bool& zero,
                    Binomial* ptr = 0) const;

    bool minimize(Binomial& b) const;

    bool auto_reduce();
    bool auto_reduce(Index& index);
    bool auto_reduce_once();
    bool auto_reduce_once(Index& index);
    bool auto_reduce_once(int first, int last, Index& index);
    bool minimal();
    bool reduced();
    bool reducable(const Binomial& b);
    bool reducable_negative(const Binomial& b);

    bool is_positive_disjoint(int i, int j) const;
    bool is_positive_disjoint(int i, const BitSet& bs) const;
    bool is_negative_disjoint(int i, int j) const;
    bool is_negative_disjoint(int i, const BitSet& bs) const;

    static bool check(const BinomialSet& bs, Binomial& b);

protected:
    BinomialSet(const BinomialSet&);
    BinomialSet& operator=(const BinomialSet& b);

    Reduction reduction;
    IndexBinomialSet binomials;

    std::vector<BitSet> pos_supps;
    std::vector<BitSet> neg_supps;
};

inline
const Binomial&
BinomialSet::operator[](Index index) const
{
    return *binomials[index];
}

inline
Size
BinomialSet::get_number() const
{
    return binomials.size();
}

inline
bool
BinomialSet::is_negative_disjoint(int i, int j) const
{
    assert(0 <= i && i < (int) binomials.size());
    assert(0 <= j && j < (int) binomials.size());
    return BitSet::set_disjoint(neg_supps[i], neg_supps[j]);
}

inline
bool
BinomialSet::is_positive_disjoint(int i, int j) const
{
    assert(0 <= i && i < (int) binomials.size());
    assert(0 <= j && j < (int) binomials.size());
    return BitSet::set_disjoint(pos_supps[i], pos_supps[j]);
}

inline
bool
BinomialSet::is_negative_disjoint(int i, const BitSet& bs) const
{
    assert(0 <= i && i < (int) binomials.size());
    return BitSet::set_disjoint(neg_supps[i], bs);
}

inline
bool
BinomialSet::is_positive_disjoint(int i, const BitSet& bs) const
{
    assert(0 <= i && i < (int) binomials.size());
    return BitSet::set_disjoint(pos_supps[i], bs);
}


inline
bool
BinomialSet::check(const BinomialSet& bs, Binomial& b)
{
    if (!Binomial::overweight(b))
    {
        b.orientate();
        bool zero = false;
        bs.reduce(b, zero);
        if (!zero && !Binomial::truncated(b)) { return true; }
    }
    return false;
}

} // namespace _4ti2_

#endif

