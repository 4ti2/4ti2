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

#ifndef _4ti2__BinomialArray_
#define _4ti2__BinomialArray_

#include "groebner/Binomial.h"
#include "groebner/VectorArray.h"
#include "groebner/BinomialCollection.h"
#include "groebner/BitSet.h"

#include <iostream>
#include "groebner/BitSetStream.h"

namespace _4ti2_
{

class BinomialArray : public BinomialCollection
{
public:
    BinomialArray();
    ~BinomialArray();

    const Binomial& operator[](Index) const;
    Binomial& operator[](Index);

    const Binomial& front() const;
    Binomial& front();
    const Binomial& back() const;
    Binomial& back();

    Size get_number() const;

    void add(const Binomial& b);
    void remove(Index i);
    void clear();

    static void transfer(BinomialArray& bs1, BinomialArray& bs2);
    static void transfer(
                    BinomialArray& bs1,
                    Index start,
                    Index end,
                    BinomialArray& bs2,
                    Index pos);
protected:
    BinomialArray(const BinomialArray&);
    BinomialArray& operator=(const BinomialArray& b);

    typedef std::vector<Binomial*> Container;
    Container binomials;
};

inline
const Binomial&
BinomialArray::operator[](Index index) const
{
    assert(0 <= index && index <= get_number());
    return *binomials[index];
}

inline
Size
BinomialArray::get_number() const
{
    return binomials.size();
}

inline
const Binomial&
BinomialArray::front() const
{
    assert(!binomials.empty());
    return *binomials.front();
}

inline
Binomial&
BinomialArray::front()
{
    assert(!binomials.empty());
    return *binomials.front();
}

inline
const Binomial&
BinomialArray::back() const
{
    assert(!binomials.empty());
    return *binomials.back();
}

inline
Binomial& 
BinomialArray::back()
{
    assert(!binomials.empty());
    return *binomials.back();
}

} // namespace _4ti2_

#endif

