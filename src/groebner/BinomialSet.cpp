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

#include "BinomialSet.h"
#include "BinomialStream.h"
#include "Statistics.h"
#include "Globals.h"

#include <iostream>
#include "BitSetStream.h"

//#define DEBUG_4ti2(X) X
#include "Debug.h"

using namespace _4ti2_;

BinomialSet::BinomialSet()
{
}

BinomialSet::~BinomialSet()
{
    for (Index i = 0; i < (Index) binomials.size(); ++i)
    {
        delete binomials[i];
    }
    binomials.clear();
}

void
BinomialSet::add(const Binomial& b)
{
    Binomial* bptr = new Binomial(b);
    binomials.push_back(bptr);
    reduction.add(*bptr);
    pos_supps.push_back(bptr->positive_support());
    neg_supps.push_back(bptr->negative_support());
}

void
BinomialSet::remove(Index i)
{
    reduction.remove(*binomials[i]);
    delete binomials[i];
    binomials.erase(binomials.begin()+i);
    pos_supps.erase(pos_supps.begin()+i);
    neg_supps.erase(neg_supps.begin()+i);
}

void
BinomialSet::clear()
{
    reduction.clear();
    for (Index i = 0; i < (Index) binomials.size(); ++i)
    {
        delete binomials[i];
    }
    binomials.clear();
    neg_supps.clear();
    pos_supps.clear();
}

// Reduce point by point set.
// Returns true if point changed.
// Sets is_zero to true if point reduces to zero.
bool
BinomialSet::reduce(
                Binomial& b,
                bool& zero,
                Binomial* bptr) const
{
    assert(!b.is_non_positive());
    DEBUG_4ti2(*out << "Reducing: " << b << "\n";)
    zero = false;
    bool changed = false;
    bool reducing = true;
    const Binomial* ptr;
    //Statistics::incr_num_reduction_steps();
    DEBUG_4ti2(unsigned int count = 0;)
    while (reducing)
    {
        reducing = false;
        ptr = reduction.reducable(b, bptr);
        if (ptr) // If ptr == 0, then no binomials reduce b.
        {
            DEBUG_4ti2(*out << "+Reducable by: " << *ptr << "\n";)
            if (!Binomial::is_negative_disjoint(b, *ptr))
            {
                DEBUG_4ti2(*out << "Reduces to 0." << "\n";)
                zero = true;
                return true;
            }
            b.reduce(*ptr);
            if (!b.orientate())
            {
                DEBUG_4ti2(*out << "Reduces to 0." << "\n";)
                zero = true;
                return true;
            }
            DEBUG_4ti2(*out << "+Reduced: " << b << "\n";)
            reducing = true;
            changed = true;
            DEBUG_4ti2(++count;)
        }
    }
    DEBUG_4ti2(
    if (count >= 100) { std::cout << "\nNum reductions = " << count << "\n"; })
    // Now reduce the negative part of the Binomial.
    reducing = true;
    while (reducing)
    {
        reducing = false;
        ptr = reduction.reducable_negative(b, bptr);
        if (ptr)
        {
            DEBUG_4ti2(*out<<"-Reducable by: " << *ptr << "\n";)
            if (!Binomial::is_pos_neg_disjoint(b, *ptr))
            {
                DEBUG_4ti2(*out << "Reduces to 0.\n";)
                zero = true;
                return true;
            }
            b.reduce_negative(*ptr);
            DEBUG_4ti2(*out << "-Reduced: " << b << "\n";)
            reducing = true;
            changed = true;
        }
    }
    if (b.is_non_positive())
    {
        // TODO: use proper error handling.
        std::cerr << "Problem is unbounded." << std::endl;
        std::cout << b << "\n";
        exit(1);
    }
    return changed;
}

// Reduce point by point set.
// Returns true if point changed.
// Sets is_zero to true if point reduces to zero.
bool
BinomialSet::reduce_negative(
                Binomial& b,
                bool& zero,
                Binomial* bptr) const
{
    assert(!b.is_non_positive());
    DEBUG_4ti2(*out << "Reducing: " << b << "\n";)
    zero = false;
    bool changed = false;
    bool reducing = true;
    const Binomial* ptr;
    //Statistics::incr_num_reduction_steps();
    // Reduce the negative part of the Binomial.
    reducing = true;
    while (reducing)
    {
        reducing = false;
        ptr = reduction.reducable_negative(b, bptr);
        if (ptr)
        {
            DEBUG_4ti2(*out<<"-Reducable by: " << *ptr << "\n";)
            if (!Binomial::is_pos_neg_disjoint(b, *ptr))
            {
                DEBUG_4ti2(*out << "Reduces to 0.\n";)
                zero = true;
                return true;
            }
            b.reduce_negative(*ptr);
            DEBUG_4ti2(*out << "-Reduced: " << b << "\n";)
            reducing = true;
            changed = true;
        }
    }
    if (b.is_non_positive())
    {
        // TODO: use proper error handling.
        std::cerr << "Problem is unbounded." << std::endl;
        std::cout << b << "\n";
        exit(1);
    }
    return changed;
}

void
BinomialSet::reducers(
                const Binomial& b,
                std::vector<const Binomial*>& reducers) const
{
    reduction.reducable(b, reducers);
}

bool
BinomialSet::reducable(const Binomial& b)
{
    return reduction.reducable(b);
}

bool
BinomialSet::reducable_negative(const Binomial& b)
{
    return reduction.reducable_negative(b);
}

bool
BinomialSet::auto_reduce_once()
{
    bool changed = false;
    Binomial b;
    for (int i = binomials.size()-1; i >= 0; --i)
    {
        b = *binomials[i];
        bool zero = false;
        bool point_changed = reduce(b, zero, binomials[i]);
        if (point_changed)
        {
            changed = true;
            remove(i);
            if (!zero) { add(b); }
        }
    }
    return changed;
}

bool
BinomialSet::auto_reduce_once(Index& index)
{
    bool changed = false;
    Binomial b;
    for (int i = binomials.size()-1; i >= 0; --i)
    {
        b = *binomials[i];
        bool zero = false;
        bool point_changed = reduce(b, zero, binomials[i]);
        if (point_changed)
        {
            if (i < index) { --index; }
            changed = true;
            remove(i);
            if (!zero) { add(b); }
        }
    }
    return changed;
}

bool
BinomialSet::auto_reduce_once(int first, int last, Index& index)
{
    assert(0 <= first && first <= last && last <= (int) binomials.size());
    bool changed = false;
    Binomial b;
    for (int i = last-1; i >= first; --i)
    {
        b = *binomials[i];
        bool zero = false;
        bool point_changed = reduce(b, zero, binomials[i]);
        if (point_changed)
        {
            if (i < index) { --index; }
            changed = true;
            remove(i);
            if (!zero) { add(b); }
        }
    }
    return changed;
}

bool
BinomialSet::auto_reduce()
{
    bool changed = false;
    while (auto_reduce_once()) { changed = true; }
    return changed;
}

bool
BinomialSet::auto_reduce(Index& index)
{
    bool changed = false;
    while (auto_reduce_once(index)) { changed = true; }
    return changed;
}

// Removes elements which are dominated by another element in the point set.
// True if point set changed.
bool
BinomialSet::minimal()
{
    bool changed = false;
    //Statistics::set_size_of_set_before_minimal(binomials.size());
    for (int i = binomials.size()-1; i >= 0; --i)
    {
        if (reduction.reducable(*binomials[i]))
        {
            changed = true;
            remove(i);
        }
    }
    return changed;
}

// Computes the reduced point set by reducing elements in the set by each other.
// Returns true if point set changed.
bool
BinomialSet::reduced()
{
    bool changed = false;
    for (int i = binomials.size()-1; i >= 0; --i)
    {
        bool reducing = true;
        const Binomial* ptr = 0;
        //Statistics::incr_num_reduction_steps();
        while (reducing)
        {
            reducing = false;
            ptr = reduction.reducable_negative(*binomials[i]);
            if (ptr)
            {
                binomials[i]->reduce_negative(*ptr);
                reducing = true;
                changed = true;
            }
        }
    }
    return changed;
}

// Minimizes binomial by set.
// Returns true if point changed.
bool
BinomialSet::minimize(Binomial& b) const
{
    bool changed = false;
    const Binomial* ptr;
    //Statistics::incr_num_reduction_steps();
    // If ptr == 0, then no binomials reduce b.
    while ((ptr = reduction.reducable(b)))
    {
        b.reduce(*ptr);
        changed = true;
    }
    return changed;
}
