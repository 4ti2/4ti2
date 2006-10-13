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

#ifndef _4ti2__Binomial_
#define _4ti2__Binomial_

#include "Vector.h"
#include "VectorArray.h"
#include "DataType.h"
#include "Statistics.h"
#include "BitSet.h"
#include "Weight.h"
#include "Grading.h"
#include "Permutation.h"
#include "VectorStream.h"
#include "Index.h"
#include "Size.h"
#include "Bounded.h"
#include "Filter.h"
#include <iostream>

namespace _4ti2_
{

class Binomial
{
public:
    Binomial();
    Binomial(const Binomial& b);
    Binomial& operator=(const Binomial& b);
    ~Binomial();

    const IntegerType& operator[](Index i) const;
    IntegerType& operator[](Index i);

    static bool reduces(
                    const Binomial& b1,
                    const Binomial& b2);
    static bool reduces(
                    const Binomial& b1,
                    const Filter& filter,
                    const Binomial& b2);
    void reduce(const Binomial& b);
    void reduce_once(const Binomial& b);
    static bool reduces_negative(
                    const Binomial& b1,
                    const Binomial& b2);
    static bool reduces_negative(
                    const Binomial& b1,
                    const Filter& filter,
                    const Binomial& b2);
    void reduce_negative(const Binomial& b);
    void reduce_once_negative(const Binomial& b);

    static bool reduces(
                    const Binomial& b0,
                    const Binomial& b1,
                    const Binomial& b2);

    static bool is_positive_disjoint(
                    const Binomial& b1,
                    const Binomial& b2);
    static bool is_negative_disjoint(
                    const Binomial& b1,
                    const Binomial& b2);
    static bool is_negative_disjoint(
                    const Binomial& b1,
                    const Binomial& b2,
                    const BitSet& is);
    static bool is_pos_neg_disjoint(
                    const Binomial& b1,
                    const Binomial& b2);
    static bool is_neg_pos_disjoint(
                    const Binomial& b1,
                    const Binomial& b2);

    IntegerType l1_norm() const;
    IntegerType l1_norm_negative() const;
    static IntegerType max_l1_norm(
                    const Binomial& b1,
                    const Binomial& b2);

    bool is_non_positive() const;
    bool is_positive(Index i) const;

    bool orientate();
    void flip();

    static void spair(
                    const Binomial& b1,
                    const Binomial& b2,
                    Binomial& b);

    void convert(Vector& v) const;

    friend bool operator==(const Binomial&, const Binomial&);
    friend bool operator!=(const Binomial&, const Binomial&);
    friend bool operator<(const Binomial&, const Binomial&);
    friend bool operator<=(const Binomial&, const Binomial&);

    void positive_support(BitSet& supp) const;
    BitSet positive_support() const;
    void negative_support(BitSet& supp) const;
    BitSet negative_support() const;
    void get_filter(Filter& filter) const;

    IntegerType degree() const;
    IntegerType degree(const Grading& grading) const;

    static bool overweight(const Binomial& b);
    static bool truncated(const Binomial& b);

    void get_positive(Vector& positive);

public:
    static Size get_num_vars();
    static Size get_num_svars();
    static Size get_num_bvars();

    static Index bnd_end;
    static Index rs_end;
    static Index urs_end;
    static Index cost_start;
    static Index cost_end;
    static Size size;

    static VectorArray* weights;
    static Weight* max_weights;
    static Vector* rhs;
    static VectorArray* lattice;
    static Grading* grading; // TODO: lazily evaluate the grading.

protected:
    IntegerType reduction_factor(const Binomial& b);
    IntegerType reduction_negative_factor(const Binomial& b);

    IntegerType* data;

    friend class BinomialFactory;
};

inline
const IntegerType&
Binomial::operator[](Index i) const
{
    assert(i >= 0 && i < size);
    return data[i];
}

inline
IntegerType&
Binomial::operator[](Index i)
{
    assert(i >= 0 && i < size);
    return data[i];
}

inline
bool
Binomial::reduces(
                const Binomial& b1,
                const Binomial& b2)
{
    for (Index i = 0; i < rs_end; ++i)
    {
        if (b1.data[i] > 0 && b1.data[i] > b2.data[i]) { return false; }
    }
    return true;
}

inline
bool
Binomial::reduces(
                const Binomial& b1,
                const Filter& filter,
                const Binomial& b2)
{
    for (Index i = 0; i < (int) filter.size(); ++i)
    {
        if (b1.data[filter[i]] > b2.data[filter[i]]) { return false; }
    }
    return true;
}

inline
void
Binomial::reduce_once(
                const Binomial& b)
{
    //Statistics::incr_num_reductions();
    for (Index i = 0; i < size; ++i) { data[i] -= b.data[i]; }
}

inline
void
Binomial::reduce(
                const Binomial& b)
{
    IntegerType factor = reduction_factor(b);
    if (factor == 1)
    {
        for (Index i = 0; i < size; ++i) { data[i] -= b.data[i]; }
    }
    else
    {
        for (Index i = 0; i < size; ++i) { data[i] -= factor*b.data[i]; }
    }
    //Statistics::incr_num_reductions();
}

inline
bool
Binomial::reduces_negative(
                const Binomial& b1,
                const Binomial& b2)
{
    for (Index i = 0; i < rs_end; ++i)
    {
        if (b1.data[i] > 0 && b1.data[i] > -b2.data[i]) { return false; }
    }
    return true;
}

inline
bool
Binomial::reduces_negative(
                const Binomial& b1,
                const Filter& filter,
                const Binomial& b2)
{
    for (Index i = 0; i < (int) filter.size(); ++i)
    {
        if (b1.data[filter[i]] > -b2.data[filter[i]]) { return false; }
    }
    return true;
}

inline
void
Binomial::reduce_once_negative(
                const Binomial& b)
{
    //Statistics::incr_num_reductions();
    for (Index i = 0; i < size; ++i) { data[i] += b.data[i]; }
}

inline
void
Binomial::reduce_negative(
                const Binomial& b)
{
    IntegerType factor = reduction_negative_factor(b);
    if (factor == -1)
    {
        for (Index i = 0; i < size; ++i) { data[i] += b.data[i]; }
    }
    else
    {
        for (Index i = 0; i < size; ++i) { data[i] -= factor*b.data[i]; }
    }
    //Statistics::incr_num_reductions();
}

inline
IntegerType
Binomial::reduction_factor(
                const Binomial& b)
{
    Index i = 0;
    while (b.data[i] <= 0) { ++i; }
    assert(i < size);
    IntegerType factor = data[i] / b.data[i];
    if (factor == 1) { return factor; }

    IntegerType tmp;
    ++i;
    for (; i < rs_end; ++i)
    {
        if (b.data[i] > 0)
        {
            tmp = data[i] / b.data[i];
            if (tmp < factor)
            {
                factor = tmp;
                if (factor == 1) { return factor; }
            }
        }
    }

    return factor;
}

inline
IntegerType
Binomial::reduction_negative_factor(
                const Binomial& b)
{
    Index i = 0;
    while (b.data[i] <= 0) { ++i; }
    assert(i < size);
    IntegerType factor = data[i] / b.data[i];
    if (factor == -1) { return factor; }

    IntegerType tmp;
    ++i;
    for (; i < rs_end; ++i)
    {
        if (b.data[i] > 0)
        {
            tmp = data[i] / b.data[i];
            if (tmp > factor)
            {
                factor = tmp;
                if (factor == -1) { return factor; }
            }
        }
    }

    return factor;
}

inline
bool
Binomial::is_positive_disjoint(
                const Binomial& b1,
                const Binomial& b2)
{
    for (Index i = 0; i < rs_end; ++i)
    {
        if (b1.data[i] > 0 && b2.data[i] > 0) { return false; }
    }
    return true;
}

inline
bool
Binomial::is_negative_disjoint(
                const Binomial& b1,
                const Binomial& b2)
{
    for (Index i = 0; i < bnd_end; ++i)
    {
        if (b1.data[i] < 0 && b2.data[i] < 0) { return false; }
    }
    return true;
}

inline
bool
Binomial::is_negative_disjoint(
                    const Binomial& b1,
                    const Binomial& b2,
                    const BitSet& is)
{
    for (Index i = 0; i < urs_end; ++i)
    {
        if (b1.data[i] < 0 && b2.data[i] < 0 && is[i]) { return false; }
    }
    return true;
}

inline
bool
Binomial::is_neg_pos_disjoint(
                const Binomial& b1,
                const Binomial& b2)
{
    for (Index i = 0; i < bnd_end; ++i)
    {
        if (b1.data[i] < 0 && b2.data[i] > 0) { return false; }
    }
    return true;
}

inline
bool
Binomial::is_pos_neg_disjoint(
                const Binomial& b1,
                const Binomial& b2)
{
    for (Index i = 0; i < bnd_end; ++i)
    {
        if (b1.data[i] > 0 && b2.data[i] < 0) { return false; }
    }
    return true;
}

inline
void
Binomial::spair(
                const Binomial& b1,
                const Binomial& b2,
                Binomial& b)
{
    for (Index i = 0; i < size; ++i) { b.data[i] = b1.data[i]-b2.data[i]; }
}

inline
IntegerType
Binomial::l1_norm() const
{
    IntegerType norm = 0;
    for (Index i = 0; i < rs_end; ++i)
    {
        if (data[i] > 0) { norm += data[i]; }
    }
    return norm;
}

inline
IntegerType
Binomial::l1_norm_negative() const
{
    IntegerType norm = 0;
    for (Index i = 0; i < rs_end; ++i)
    {
        if (data[i] < 0) { norm -= data[i]; }
    }
    return norm;
}

inline
IntegerType
Binomial::max_l1_norm(
                const Binomial& b1,
                const Binomial& b2)
{
    IntegerType norm = 0;
    for (Index i = 0; i < rs_end; ++i)
    {
        if (b1.data[i] > 0 && b1.data[i]>=b2.data[i]) { norm += b1.data[i]; }
        else if (b2.data[i] > 0) { norm += b2.data[i]; }
    }
    return norm;
}

// Returns false if binomial is zero.
// Will orientate the binomial according to degrevlex
// x_1 < x_2 < ... < x_n
inline
bool
Binomial::orientate()
{
    Index i = cost_start;
    while (i < cost_end && data[i] == 0) { ++i; }
    if (i == cost_end)
    {
        i = 0;
        while (i < rs_end && data[i] == 0) { ++i; }
        if (i == rs_end) { return false; } // the binomial is zero.
        else if (data[i] > 0) { flip(); }
    }
    else if (data[i] < 0) { flip(); }
    return true;
}

inline
void
Binomial::flip()
{
    for (Index i = 0; i < size; ++i) { data[i] *= -1; }
}

inline
bool
Binomial::is_non_positive() const
{
    for (Index i = 0; i < rs_end; ++i)
    {
        if (data[i] > 0) { return false; }
    }
    return true;
}

inline
bool
Binomial::is_positive(Index i) const
{
    assert(i >= 0 && i < urs_end);
    return data[i] > 0;
}

inline
Size
Binomial::get_num_vars()
{
    return urs_end;
}

inline
Size
Binomial::get_num_svars()
{
    return rs_end;
}

inline
Size
Binomial::get_num_bvars()
{
    return bnd_end;
}

#if 0
inline
bool 
operator==(const Binomial& b1, const Binomial& b2)
{
    for (Index i = 0; i < Binomial::urs_end; ++i)
    {
        if (b1.data[i] != b2.data[i]) { return false; }
    }
    return true;
}
#endif

inline
bool
operator!=(const Binomial& b1, const Binomial& b2)
{
    return !(b1 == b2);
}

// Lexicographic ordering.
inline
bool
operator<(const Binomial& b1, const Binomial& b2)
{
    Index i = 0;
    while (i < Binomial::urs_end && b1.data[i] == b2.data[i])
    {
        ++i;
    }
    if (i < Binomial::urs_end && b1.data[i] < b2.data[i])
    {
        return true;
    }
    return false;
}

// Lexicographic ordering.
inline
bool
operator<=(const Binomial& b1, const Binomial& b2)
{
    return !(b2 < b1);
}

// Lexicographic ordering.
inline
bool
operator>(const Binomial& b1, const Binomial& b2)
{
    return b2 < b1;
}

// Lexicographic ordering.
inline
bool
operator>=(const Binomial& b1, const Binomial& b2)
{
    return b2 <= b1;
}

inline
Binomial::Binomial()
{
    data = new IntegerType[size];
}

inline
Binomial::Binomial(const Binomial& b)
{
    data = new IntegerType[size];
    *this = b;
}

inline
Binomial::~Binomial()
{
    delete [] data;
}

inline
Binomial&
Binomial::operator=(const Binomial& b)
{
    for (Index i = 0; i < size; ++i) { data[i] = b.data[i]; }
    return *this;
}

inline
void
Binomial::positive_support(BitSet& supp) const
{
    // Assuming supp is all zeros.
    assert(supp.get_size() == rs_end);
    for (Index i = 0; i < rs_end; ++i)
    {
        if (data[i] > 0) { supp.set(i); }
    }
}

inline
BitSet
Binomial::positive_support() const
{
    BitSet supp(rs_end);
    for (Index i = 0; i < rs_end; ++i)
    {
        if (data[i] > 0) { supp.set(i); }
    }
    return supp;
}

inline
void
Binomial::get_filter(Filter& filter) const
{
    for (Index i = 0; i < rs_end; ++i)
    {
        if (data[i] > 0) { filter.push_back(i); }
    }
}

inline
void
Binomial::negative_support(BitSet& supp) const
{
    // Assuming supp is all zeros.
    assert(supp.get_size() == bnd_end);
    for (Index i = 0; i < bnd_end; ++i)
    {
        if (data[i] < 0) { supp.set(i); }
    }
}

inline
BitSet
Binomial::negative_support() const
{
    BitSet supp(bnd_end);
    for (Index i = 0; i < bnd_end; ++i)
    {
        if (data[i] < 0) { supp.set(i); }
    }
    return supp;
}

inline
IntegerType
Binomial::degree(const Grading& g) const
{
    // TODO: What is the right length?
    assert(g.get_size() == urs_end);
    IntegerType d = 0;
    for (Index i = 0; i < rs_end; ++i)
    {
        if (data[i] > 0) { d += data[i]*g[i]; }
    }
    return d;
}

inline
IntegerType
Binomial::degree() const
{
    return degree(*grading);
}

// True if (b0+ wedge b1+) reduces (b0+ wedge b2+)
inline
bool
Binomial::reduces(
                const Binomial& b0,
                const Binomial& b1,
                const Binomial& b2)
{
    for (Index i = 0; i < rs_end; ++i)
    {
        if (b1.data[i]>0 && b1.data[i]>b2.data[i] && b1.data[i]>b0.data[i])
        {
            return false;
        }
    }
    return true;
}

inline
bool
Binomial::overweight(const Binomial& b)
{
    if (max_weights != 0)
    {
        assert(weights != 0);
        assert(max_weights->get_size() == weights->get_number());
        const Weight& max = *max_weights;
        const VectorArray& ws = *weights;
        for (Index i = 0; i < ws.get_number(); ++i)
        {
            if (b.degree(ws[i]) > max[i]) { return true; }
        }
    }
    return false;
}

} // namespace _4ti2_

#endif
