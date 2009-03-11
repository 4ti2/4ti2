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

#ifndef _4ti2_groebner__Vector_
#define _4ti2_groebner__Vector_

#include "groebner/Permutation.h"
#include "groebner/DataType.h"
#include <cassert>
#include "groebner/Index.h"
#include "groebner/Size.h"

#include <iostream>

namespace _4ti2_
{

// TODO: It is probably quicker when using gmp that we always pass IntegerType
// as a const reference instead of pass by value. It does not seem to have any
// effect when we are not using gmp.

class Vector
{
public:
    Vector();
    Vector(const Vector&);
    explicit Vector(Size size);
    Vector(Size size, IntegerType value);
    Vector& operator=(const Vector&);
    ~Vector();

    const Size& get_size() const;

    const IntegerType& operator[](Index) const;
    IntegerType& operator[](Index);

    void permute(const Permutation& p);

    void normalise();

    // TODO: sort out which operations are needed.
    static void add(const Vector&,
                    const Vector&,
                    Vector&);
    static void add(const Vector&,
                    IntegerType,
                    const Vector&,
                    IntegerType,
                    Vector&);
    static void sub(const Vector&,
                    const Vector&,
                    Vector&);
    static void sub(const Vector&,
                    IntegerType,
                    const Vector&,
                    IntegerType,
                    Vector&);
    static void mul(const Vector&,
                    IntegerType,
                    Vector&);
    static void div(const Vector&,
                    IntegerType,
                    Vector&);

    friend Vector operator+ (
                    const Vector&,
                    const Vector&);
    friend Vector operator- (
                    const Vector&,
                    const Vector&);
    friend Vector operator* (
                    const Vector&,
                    IntegerType);
    friend Vector operator/ (
                    const Vector&,
                    IntegerType);
    friend Vector operator+ (
                    const Vector&);
    friend Vector operator- (
                    const Vector&);
    void operator+=(const Vector&);
    void operator-=(const Vector&);
    void operator*=(IntegerType);
    void operator/=(IntegerType);

    friend bool operator==(const Vector&, const Vector&);
    friend bool operator!=(const Vector&, const Vector&);
    friend bool operator< (const Vector&, const Vector&);
    friend bool operator<=(const Vector&, const Vector&);
    friend bool operator> (const Vector&, const Vector&);
    friend bool operator>=(const Vector&, const Vector&);

    static IntegerType dot(
                    const Vector&,
                    const Vector&);
    static void dot(const Vector&,
                    const Vector&,
                    IntegerType&);
    static void concat(
                    const Vector& v1,
                    const Vector& v2,
                    Vector& v);
    template <class IndexSet>
    static void concat(
                    const Vector& v1,
                    const Vector& v2,
                    const IndexSet& s,
                    Vector& v);
    static void split(
                    const Vector& v,
                    Vector& v1,
                    Vector& v2);
    template <class IndexSet>
    static void split(
                    const Vector& v,
                    const IndexSet& s,
                    Vector& v1,
                    Vector& v2);
    static void project(
                    const Vector& v1,
                    Index start,
                    Index end,
                    Vector& v);
    template <class IndexSet>
    static void project(
                    const Vector& v1,
                    const IndexSet& proj,
                    Vector& v);
    template <class IndexSet>
    static void lift(
                    const Vector& v1,
                    const IndexSet& lift,
                    Vector& v);
    static void lift(
                    const Vector& v1,
                    Index start,
                    Index end,
                    Vector& v);
    static void minimum(
                    const Vector& v1,
                    const Vector& v2,
                    Vector& v);

    // TODO: We should not need this function.
    template <class IndexSet>
    void project(const IndexSet& proj);

    void add(const Vector&);
    void add(const Vector&, IntegerType);
    void sub(const Vector&);
    void sub(const Vector&, IntegerType);
    void mul(IntegerType);
    void div(IntegerType);

    bool is_zero() const;
    bool is_non_negative() const;
    bool is_non_positive() const;

    template <class IndexSet>
    bool is_zero(const IndexSet& proj) const;
    template <class IndexSet>
    bool is_non_negative(const IndexSet& proj) const;
    template <class IndexSet>
    bool is_non_positive(const IndexSet& proj) const;

protected:
    IntegerType *vector;
    Size size;
};

inline
const IntegerType&
Vector::operator[](Index index) const
{
    assert(index >= 0 && index < size);
    return vector[index];
}

inline
IntegerType&
Vector::operator[](Index index)
{
    assert(index >= 0 && index < size);
    return vector[index];
}

inline
const Size&
Vector::get_size() const
{
    return size;
}

inline
Vector 
operator+(const Vector& v1, const Vector& v2)
{
    Vector r(v1.size);
    Vector::add(v1,v2,r);
    return r;
}

inline
Vector 
operator-(const Vector& v1, const Vector& v2)
{
    Vector r(v1.size);
    Vector::sub(v1,v2,r);
    return r;
}

inline
Vector 
operator*(const Vector& v1, IntegerType m)
{
    Vector r(v1.size);
    Vector::mul(v1,m,r);
    return r;
}

inline
Vector 
operator/(const Vector& v1, IntegerType d)
{
    Vector r(v1.size);
    Vector::div(v1,d,r);
    return r;
}

inline
Vector 
operator+(const Vector& v)
{
    return v;
}

inline
Vector 
operator-(const Vector& v)
{
    Vector r(v.size);
    Vector::mul(v, -1, r);
    return r;
}

inline
void 
Vector::operator+=(const Vector& v)
{
    add(v);
}

inline
void 
Vector::operator-=(const Vector& v)
{
    sub(v);
}

inline
void 
Vector::operator*=(IntegerType m)
{
    mul(m);
}

inline
void 
Vector::operator/=(IntegerType d)
{
    mul(d);
}

inline
bool 
operator==(const Vector& v1, const Vector& v2)
{
    assert(v1.size == v2.size);
    for (Index i = 0; i < v1.size; ++i)
    {
        if (v1.vector[i] != v2.vector[i]) { return false; }
    }

    return true;
}

inline
bool 
operator!=(const Vector& v1, const Vector& v2)
{
    return !(v1 == v2);
}

// Lexicographic ordering.
inline
bool
operator<(const Vector& v1, const Vector& v2)
{
    assert(v1.size == v2.size);
    Index i = 0;
    while (i < v1.size && v1.vector[i] == v2.vector[i]) { ++i; }
    if (i < v1.size && v1.vector[i] < v2.vector[i]) { return true; }
    return false;
}

// Lexicographic ordering.
inline
bool
operator<=(const Vector& v1, const Vector& v2)
{
    return !(v2 < v1);
}

// Lexicographic ordering.
inline
bool
operator>(const Vector& v1, const Vector& v2)
{
    return v2 < v1;
}

// Lexicographic ordering.
inline
bool
operator>=(const Vector& v1, const Vector& v2)
{
    return v2 <= v1;
}

inline
void
Vector::sub(
                const Vector& v)
{
    assert(v.size == size);
    for (Index i = 0; i < size; ++i) { vector[i] -= v.vector[i]; }
}

inline
void
Vector::sub(
                const Vector& v,
                IntegerType mul)
{
    assert(v.size == size);
    for (Index i = 0; i < size; ++i) { vector[i] -= mul*v.vector[i]; }
}

inline
void
Vector::sub(
                const Vector& v1,
                const Vector& v2,
                Vector& r)
{
    assert(v1.size == v2.size);
    assert(v1.size == r.size);
    for (Index i = 0; i < v1.size; ++i)
    {
        r.vector[i] = v1.vector[i]-v2.vector[i];
    }
}

inline
void
Vector::mul(IntegerType m)
{
    for (Index i = 0; i < size; ++i) { vector[i] *= m; }
}

inline
IntegerType
Vector::dot(
                const Vector& v1,
                const Vector& v2)
{
    IntegerType r;
    Vector::dot(v1, v2, r);
    return r;
}

inline
void
Vector::dot(
                const Vector& v1,
                const Vector& v2,
                IntegerType& r)
{
    assert(v1.size == v2.size);
    r = 0;
    for (Index i = 0; i < v1.size; ++i) { r += v1[i]*v2[i]; }
}

inline
void
Vector::concat( const Vector& v1,
                const Vector& v2,
                Vector& v)
{
    assert(v1.get_size()+v2.get_size() == v.get_size());
    for (Index i = 0; i < v1.get_size(); ++i) { v[i] = v1[i]; }
    for (Index i = 0; i < v2.get_size(); ++i) { v[i+v1.get_size()] = v2[i]; }
}

template <class IndexSet>
inline
void
Vector::concat( const Vector& v1,
                const Vector& v2,
                const IndexSet& s,
                Vector& v)
{
    assert(v1.get_size() == s.count());
    assert(v1.get_size()+v2.get_size() == v.get_size());
    assert(v.get_size() == s.get_size());
    int index1 = 0;
    int index2 = 0;
    for (Index i = 0; i < v.get_size(); ++i)
    {
        if (s[i]) { v[i] = v1[index1]; ++index1; }
        else { v[i] = v2[index2]; ++index2; }
    } 
}

inline
void
Vector::split(  const Vector& v,
                Vector& v1,
                Vector& v2)
{
    assert(v1.get_size()+v2.get_size() == v.get_size());
    for (int i = 0; i < v1.get_size(); ++i) { v1[i] = v[i]; }
    for (int i = 0; i < v2.get_size(); ++i) { v2[i] = v[i+v1.get_size()]; }
}

template <class IndexSet>
inline
void
Vector::split(  const Vector& v,
                const IndexSet& s,
                Vector& v1,
                Vector& v2)
{
    assert(v1.get_size()+v2.get_size() == v.get_size());
    assert(v1.get_size() == s.count());
    int index1 = 0;
    int index2 = 0;
    for (Index i = 0; i < v.get_size(); ++i)
    {
        if (s[i]) { v1[index1] = v[i]; ++index1; }
        else { v2[index2] = v[i]; ++index2; }
    } 
}

inline
void
Vector::project(const Vector& v1,
                Index start,
                Index end,
                Vector& v)
{
    assert(start <= end && start < v1.get_size());
    assert(v.get_size() == end-start);
    for (Index i = 0; i < v.get_size(); ++i) { v[i] = v1[i+start]; } 
}

template <class IndexSet>
inline
void
Vector::project(const Vector& v1,
                const IndexSet& proj,
                Vector& v)
{
    assert(v1.get_size() == proj.get_size());
    assert(v.get_size() == proj.count());
    int index = 0;
    for (Index i = 0; i < v1.get_size(); ++i)
    {
        if (proj[i]) { v[index] = v1[i]; ++index; }
    } 
}

template <class IndexSet>
inline
void
Vector::project(const IndexSet& proj)
{
    int index = 0;
    for (Index i = 0; i < size; ++i)
    {
        if (proj[i]) { vector[index] = vector[i]; ++index; }
    }
    size = index; // I do not free up unused memory.
}

template <class IndexSet>
inline
void
Vector::lift(   const Vector& v1,
                const IndexSet& lift,
                Vector& v)
{
    assert(v1.get_size() == lift.count());
    assert(v.get_size() == lift.get_size());
    int index = 0;
    for (Index i = 0; i < v.get_size(); ++i)
    {
        if (lift[i]) { v[i] = v1[index]; ++index; }
    } 
}

inline
void
Vector::lift(   const Vector& v1,
                Index start,
                Index end,
                Vector& v)
{
    assert(v1.get_size() == end-start);
    assert(start <= end && end <= v.get_size());
    for (Index i = 0; i < v1.get_size(); ++i)
    {
        v[i+start] = v1[i];
    } 
}

inline
void
Vector::minimum(
                    const Vector& v1,
                    const Vector& v2,
                    Vector& v)
{
    assert(v2.get_size() == v1.get_size() && v.get_size() == v1.get_size());
    for (Index i = 0; i < v.get_size(); ++i)
    {
        if (v1[i] <= v2[i]) { v[i] = v1[i]; }
        else { v[i] = v2[i]; }
    } 
}

inline
bool
Vector::is_zero() const
{
    Index i = 0;
    while (i < size && vector[i] == 0) { ++i; }
    if (i < size) { return false; }
    return true;
}

template <class IndexSet>
inline
bool
Vector::is_zero(const IndexSet& proj) const
{
    Index i = 0;
    while (i < size && (!proj[i] || vector[i] == 0)) { ++i; }
    if (i < size) { return false; }
    return true;
}

inline
bool
Vector::is_non_negative() const
{
    for (Index i = 0; i < size; ++i)
    {
        if (vector[i] < 0) { return false; }
    }
    return true;
}

template <class IndexSet>
inline
bool
Vector::is_non_negative(const IndexSet& proj) const
{
    for (Index i = 0; i < size; ++i)
    {
        if (proj[i] && vector[i] < 0) { return false; }
    }
    return true;
}

inline
bool
Vector::is_non_positive() const
{
    for (Index i = 0; i < size; ++i)
    {
        if (vector[i] > 0) { return false; }
    }
    return true;
}

template <class IndexSet>
inline
bool
Vector::is_non_positive(const IndexSet& proj) const
{
    for (Index i = 0; i < size; ++i)
    {
        if (proj[i] && vector[i] > 0) { return false; }
    }
    return true;
}

inline
void
Vector::add(
                const Vector& v1,
                const Vector& v2,
                Vector& r)
{
    assert(v1.size == v2.size);
    assert(v1.size == r.size);
    for (Index i = 0; i < v1.size; ++i)
    {
        r.vector[i] = v1.vector[i]+v2.vector[i];
    }
}

inline
void
Vector::mul(
                const Vector& v,
                IntegerType m,
                Vector& r)
{
    assert(v.size == r.size);
    for (Index i = 0; i < v.size; ++i) { r.vector[i] = m*v.vector[i]; }
}

inline
void
Vector::div(
                const Vector& v,
                IntegerType d,
                Vector& r)
{
    assert(v.size == r.size);
    for (Index i = 0; i < v.size; ++i) { r.vector[i] = v.vector[i]/d; }
}

// r = m1*v1 + m2*v2
inline
void
Vector::add(
                const Vector& v1,
                IntegerType m1, 
                const Vector& v2,
                IntegerType m2, 
                Vector& r)
{
    assert(v1.size == v2.size);
    assert(v1.size == r.size);
    for (Index i = 0; i < v1.size; ++i)
    {
        r.vector[i] = m1*v1.vector[i]+m2*v2.vector[i];
    }
}

// r = m1*v1 - m2*v2
inline
void
Vector::sub(
                const Vector& v1,
                IntegerType m1, 
                const Vector& v2,
                IntegerType m2, 
                Vector& r)
{
    assert(v1.size == v2.size);
    assert(v1.size == r.size);
    for (Index i = 0; i < v1.size; ++i)
    {
        r.vector[i] = m1*v1.vector[i]-m2*v2.vector[i];
    }
}

inline
void
Vector::add(const Vector& v)
{
    assert(v.size == size);
    for (Index i = 0; i < size; ++i) { vector[i] += v.vector[i]; }
}

inline
void
Vector::add(const Vector& v, IntegerType mul)
{
    assert(v.size == size);
    for (Index i = 0; i < size; ++i) { vector[i] += mul*v.vector[i]; }
}

inline
void
Vector::div(IntegerType d)
{
    for (Index i = 0; i < size; ++i) { vector[i] /= d; }
}

inline
Vector&
Vector::operator=(const Vector& v)
{
    assert(size == v.size);
    for (Index i = 0; i < size; ++i) { vector[i] = v.vector[i]; }
    return *this;
}

} // namespace 4ti2

#endif
