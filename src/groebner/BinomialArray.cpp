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

#include "BinomialArray.h"
#include "Statistics.h"
#include "Globals.h"

#include <iostream>
#include "BitSetStream.h"

#define DEBUG_4ti2(X) //X
#include "Debug.h"


using namespace _4ti2_;

BinomialArray::BinomialArray()
{
}

BinomialArray::~BinomialArray()
{
    for (Index i = 0; i < (Index) binomials.size(); ++i)
    {
        delete binomials[i];
    }
    binomials.clear();
}

void
BinomialArray::add(const Binomial& b)
{
    Binomial* bptr = new Binomial(b);
    binomials.push_back(bptr);
}

void
BinomialArray::remove(Index i)
{
    delete binomials[i];
    binomials.erase(binomials.begin()+i);
}

void
BinomialArray::clear()
{
    for (Index i = 0; i < (Index) binomials.size(); ++i)
    {
        delete binomials[i];
    }
    binomials.clear();
}

void
BinomialArray::transfer(
                BinomialArray& bs1,
                BinomialArray& bs2)
{
    transfer(bs1, 0, bs1.get_number(), bs2, bs2.get_number());
}

void
BinomialArray::transfer(
                BinomialArray& bs1,
                Index start,
                Index end,
                BinomialArray& bs2,
                Index pos)
{
    assert(start >= 0 && start <= end);
    assert(end <= bs1.get_number());
    assert(pos >= 0 && pos <= bs2.get_number());
    bs2.binomials.insert(bs2.binomials.begin()+pos,
                    bs1.binomials.begin()+start, bs1.binomials.begin()+end);
    bs1.binomials.erase(bs1.binomials.begin()+start, bs1.binomials.begin()+end);
}
