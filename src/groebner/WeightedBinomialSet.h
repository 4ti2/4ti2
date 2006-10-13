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

#ifndef _4ti2__WeightedBinomialSet_
#define _4ti2__WeightedBinomialSet_

#include "Binomial.h"
#include "BinomialCollection.h"
#include "Grading.h"
#include <map>
#include <set>

namespace _4ti2_ 
{

class WeightedBinomialSet : public BinomialCollection
{
public:
    WeightedBinomialSet();
    ~WeightedBinomialSet();

    void add(const Binomial& b);
    void next(Binomial& b);
    void clear();
    bool empty() const;
    Size get_size() const;

    const Grade min_grade() const;
    const Grade max_grade() const;

    void print() const;

protected:
    typedef std::set<std::pair<std::pair<Grade,Grade>,Binomial> > Container;
    //typedef std::multimap<Grade,Binomial> Container;
    Container data;
};

inline
Size
WeightedBinomialSet::get_size() const
{
    return data.size();
}

inline
bool
WeightedBinomialSet::empty() const
{
    return data.empty();
}

inline
const Grade
WeightedBinomialSet::min_grade() const
{
    if (!data.empty()) { return data.begin()->first.first; }
    return 0;
}

inline
const Grade
WeightedBinomialSet::max_grade() const
{
    if (!data.empty()) { return (--data.end())->first.first; }
    return 0;
}

} // namespace _4ti2_

#endif

