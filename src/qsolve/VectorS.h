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

#ifndef _4ti2_qsolve__VectorS_
#define _4ti2_qsolve__VectorS_

#include "qsolve/Index.h"
#include "qsolve/Size.h"
#include "qsolve/Euclidean.h"

#include <vector>
#include <cassert>
#include <algorithm>

// TODO: Remove the following dependency.
#include <iostream>

namespace _4ti2_
{

// TODO: It is probably quicker when using gmp that we always pass T
// as a const reference instead of pass by value. It does not seem to have any
// effect when we are not using gmp.

template <class T>
class VectorS
{
public:
    VectorS();
    VectorS(const VectorS<T>&);
    explicit VectorS(Size size);
    VectorS(Size size, T value);
    VectorS<T>& operator=(const VectorS<T>&);
    ~VectorS();

    Size get_size() const;

    const T operator[](Index) const;
    void set(Index, T);

    static void add(const VectorS<T>&, const VectorS<T>&, VectorS<T>&);
    static void add(const VectorS<T>&, const VectorS<T>&, T, VectorS<T>&);
    static void add(const VectorS<T>&, T, const VectorS<T>&, T, VectorS<T>&);
    static void sub(const VectorS<T>&, const VectorS<T>&, VectorS<T>&);
    static void mul(const VectorS<T>&, T, VectorS<T>&);
    static void div(const VectorS<T>&, T, VectorS<T>&);

    void add(const VectorS<T>&);
    void add(const VectorS<T>&, T);
    void add(T, const VectorS<T>&, T);
    void sub(const VectorS<T>&);
    void mul(T);
    void div(T);

    bool operator==(const VectorS<T>&) const;
    bool operator!=(const VectorS<T>&) const;
    bool operator< (const VectorS<T>&) const;

    static void swap(VectorS<T>&, VectorS<T>&);

    static T dot( const VectorS<T>&, const VectorS<T>&);
    static void dot(const VectorS<T>&, const VectorS<T>&, T&);

    //template <class IndexSet>
    //static void project(const VectorS<T>& v1, const IndexSet& proj, VectorS<T>& v);
    template <class IndexSet1, class IndexSet2>
    void assign(const VectorS<T>& v1, const IndexSet1& is1, const IndexSet2& is2);

    void normalise();
    void print();

protected:
    // An exponent.
    struct Exp {
        public:
            Exp() : i(0), v(0) {}
            Exp(Index index, T value) : i(index), v(value) {}
            Index i;
            T v;
            bool operator<(Exp e) const;
    };
    typedef typename std::vector<Exp> Data;
    typedef typename Data::iterator Iter;
    typedef typename Data::const_iterator CIter;
    Data data;
    Size size;
};

template <class T>
inline
bool
VectorS<T>::Exp::operator<(Exp e) const
{
    if (i < e.i) { return true; }
    return false;
}

template <class T>
inline
const T
VectorS<T>::operator[](Index index) const
{
    assert(index >= 0 && index < size);
    CIter it = std::lower_bound(data.begin(), data.end(), Exp(index,0));
#if 0
    CIter it = data.begin();
    while (it != data.end() && it->i < index) { ++it; }
#endif
    if (it != data.end() && it->i == index) { return it->v; }
    return 0;
}

template <class T>
inline void
VectorS<T>::set(Index index, T v)
{
    assert(index >= 0 && index < size);
    Iter it = std::lower_bound(data.begin(), data.end(), Exp(index,0));
#if 0
    Iter it = data.begin();
    while (it != data.end() && it->i < index) { ++it; }
#endif
    if (it != data.end() && it->i == index) {
        if (v == 0) { data.erase(it); }
        else { it->v = v; }
    }
    if (v != 0) { data.insert(it, Exp(index, v)); }
}

template <class T>
inline
Size
VectorS<T>::get_size() const
{
    return size;
}

template <class T>
inline
bool 
VectorS<T>::operator==(const VectorS<T>& v2) const
{
    assert(size == v2.size);
    if (data.size() != v2.data.size()) { return false; }
    CIter it = data.begin();
    CIter it2 = v2.data.begin();
    while (it != data.end()) {
        if (it->i != it2->i || it->v != it2->v) { return false; }
        ++it; ++it2;
    }
    return true;
}

template <class T>
inline
bool 
VectorS<T>::operator!=(const VectorS<T>& v2) const
{
    assert(size == v2.size);
    return !(*this == v2);
}

// Lexicographic ordering.
template <class T>
inline
bool
VectorS<T>::operator<(const VectorS<T>& v2) const
{
    assert(size == v2.size);
    CIter it = data.begin();
    CIter it2 = v2.data.begin();
    while (it != data.end() && it2 != v2.data.end()) {
        if (it->i != it2->i || it->v != it2->v) { break; }
        ++it; ++it2;
    }
    if (it2 == v2.data.end()) { return false; }
    if (it == data.end()) { return true; }
    if (it->i < it2->i) { return true; }
    if (it->i == it2->i && it->v < it2->v) { return true; }
    return false;
}

template <class T>
inline
void
VectorS<T>::swap(VectorS<T>& v1, VectorS<T>& v2)
{
    assert(v1.size == v2.size);
    v1.data.swap(v2.data);
}

template <class T>
inline
T
VectorS<T>::dot(const VectorS<T>& v1, const VectorS<T>& v2)
{
    assert(v1.size == v2.size);
    T r;
    VectorS<T>::dot(v1, v2, r);
    return r;
}

template <class T>
inline
void
VectorS<T>::dot(const VectorS<T>& v1, const VectorS<T>& v2, T& r)
{
    assert(v1.size == v2.size);
    CIter it1 = v1.data.begin();
    CIter it2 = v2.data.begin();
    r = 0;
    while (it1 != v1.data.end() && it2 != v2.data.end()) {
        if (it1->i == it2->i) { r += it1->v * it2->v; ++it1; ++it2; }
        else if (it1->i < it2->i) { ++it1; }
        else { ++it2; }
    }
}

#if 0
template <class T> template <class IndexSet>
inline
void
VectorS<T>::project(const VectorS<T>& v1,
                const IndexSet& proj,
                VectorS<T>& v)
{
    assert(v.size == proj.count() && v.size <= v1.size);
    v.data.clear();
    typename IndexSet::Iter p  = proj.begin();
    Index i = 0;
    CIter it1 = v1.data.begin();
    while (p != proj.end()) {
        while (it1 != v1.data.end() && it1->i < *p) { ++it1; }
        if (it1 == v1.data.end()) { return; }
        if (it1->i == *p) { v.data.push_back(Exp(i, it1->v)); }
        ++p; ++i;
    }
}
#endif

template <class T>
inline
void
VectorS<T>::add(const VectorS<T>& v1, const VectorS<T>& v2, VectorS<T>& r)
{
    assert(v1.size == v2.size && v1.size == r.size);
    r.data.clear();

    CIter it1 = v1.data.begin();
    CIter it2 = v2.data.begin();
    while (it1 != v1.data.end() && it2 != v2.data.end()) {
        if (it1->i < it2->i) {
            r.data.push_back(*it1); ++it1;
        } else if (it2->i < it1->i) {
            r.data.push_back(*it2); ++it2;
        } else {
            T v = it1->v+it2->v;
            if (v!=0) { r.data.push_back(Exp(it1->i, v)); }
            ++it1; ++it2;
        }
    }

    r.data.insert(r.data.end(), it1, v1.data.end());
    r.data.insert(r.data.end(), it2, v2.data.end());
}

template <class T>
inline
void
VectorS<T>::sub(const VectorS<T>& v1, const VectorS<T>& v2, VectorS<T>& r)
{
    assert(v1.size == v2.size && v1.size == r.size);
    r.data.clear();

    CIter it1 = v1.data.begin();
    CIter it2 = v2.data.begin();
    while (it1 != v1.data.end() && it2 != v2.data.end()) {
        if (it1->i < it2->i) {
            r.data.push_back(*it1); ++it1;
        } else if (it2->i < it1->i) {
            r.data.push_back(Exp(it2->i,-it2->v)); ++it2;
        } else {
            T v = it1->v-it2->v;
            if (v!=0) { r.data.push_back(Exp(it1->i, v)); }
            ++it1; ++it2;
        }
    }

    r.data.insert(r.data.end(), it1, v1.data.end());
    while (it2 != v2.data.end()) {
        r.data.push_back(Exp(it2->i,-it2->v));
        ++it2;
    }
}

// This function assumes that v*m==0 iff m==0 or v==0.
template <class T>
inline
void
VectorS<T>::mul(const VectorS<T>& v1, T m, VectorS<T>& r)
{
    assert(v1.size == r.size);
    r.data.clear();
    if (m == 0) { return; }
    for (CIter it1 = v1.data.begin(); it1 != v1.data.end(); ++it1) {
        r.data.push_back(Exp(it1->i, m*it1->v));
    }
}

template <class T>
inline
void
VectorS<T>::div(const VectorS<T>& v1, T d, VectorS<T>& r)
{
    assert(v1.size == r.size);
    assert(d!=0);
    r.data.clear();
    for (CIter it1 = v1.data.begin(); it1 != v1.data.end(); ++it1) {
        T v = it1->v / d;
        if (v!=0) { r.data.push_back(Exp(it1->i, v)); }
    }
}

// r = v1 + m2*v2
template <class T>
inline
void
VectorS<T>::add(
                const VectorS<T>& v1,
                const VectorS<T>& v2,
                T m2, 
                VectorS<T>& r)
{
    r.data.clear();
    CIter it1 = v1.data.begin();
    CIter it2 = v2.data.begin();
    while (it1 != v1.data.end() && it2 != v2.data.end())
    {
        if (it1->i < it2->i) {
            r.data.push_back(*it1); ++it1;
        } else if (it2->i < it1->i) {
            r.data.push_back(Exp(it2->i,m2*it2->v)); ++it2;
        } else {
            T v = it1->v+m2*it2->v;
            if (v!=0) { r.data.push_back(Exp(it1->i, v)); }
            ++it1; ++it2;
        }
    }

    r.data.insert(r.data.end(), it1, v1.data.end());
    while (it2 != v2.data.end()) {
        r.data.push_back(Exp(it2->i,m2*it2->v));
        ++it2;
    }
}

// r = m1*v1 + m2*v2
template <class T>
inline
void
VectorS<T>::add(
                const VectorS<T>& v1,
                T m1, 
                const VectorS<T>& v2,
                T m2, 
                VectorS<T>& r0)
{
    //r.data.clear();
    VectorS<T> r(r0.size);
    CIter it1 = v1.data.begin();
    CIter it2 = v2.data.begin();
    while (it1 != v1.data.end() && it2 != v2.data.end())
    {
        if (it1->i < it2->i) {
            r.data.push_back(Exp(it1->i,m1*it1->v)); ++it1;
        } else if (it2->i < it1->i) {
            r.data.push_back(Exp(it2->i,m2*it2->v)); ++it2;
        } else {
            T v = m1*it1->v+m2*it2->v;
            if (v!=0) { r.data.push_back(Exp(it1->i, v)); }
            ++it1; ++it2;
        }
    }

    while (it1 != v1.data.end()) {
        r.data.push_back(Exp(it1->i,m1*it1->v));
        ++it1;
    }
    while (it2 != v2.data.end()) {
        r.data.push_back(Exp(it2->i,m2*it2->v));
        ++it2;
    }
    r0 = r;
}

template <class T>
inline
void
VectorS<T>::add(const VectorS<T>& v1)
{
    assert(size == v1.size);
    VectorS<T> r(size);
    VectorS<T>::add(*this, v1, r);
    VectorS<T>::swap(*this, r);
}

template <class T>
inline
void
VectorS<T>::add(T m, const VectorS<T>& v1, T m1)
{
    VectorS<T> r(size);
    VectorS<T>::add(*this, m, v1, m1, r);
    VectorS<T>::swap(*this, r);
}

template <class T>
inline
void
VectorS<T>::add(const VectorS<T>& v1, T mul)
{
    VectorS<T> r(size);
    VectorS<T>::add(*this, v1, mul, r);
    VectorS<T>::swap(*this, r);
}

template <class T>
inline
void
VectorS<T>::sub(const VectorS<T>& v1)
{
    VectorS<T> r(size);
    VectorS<T>::sub(*this, v1, r);
    VectorS<T>::swap(*this, r);
}

// This function assumes that v*m==0 iff m==0 or v==0.
template <class T>
inline
void
VectorS<T>::mul(T m)
{
    if (m==0) { data.clear(); return; }
    for (Iter it = data.begin(); it != data.end(); ++it) { it->v *= m; }
}

template <class T>
inline
void
VectorS<T>::div(T d)
{
    Iter it0 = data.begin();
    for (Iter it = data.begin(); it != data.end(); ++it) {
        it0->v = it->v / d;
        if (it0->v != 0) { ++it0; }
    }
    data.erase(it0, data.end());
}

template <class T>
inline
VectorS<T>&
VectorS<T>::operator=(const VectorS<T>& v)
{
    data = v.data;
    return *this;
}

template <class T>
inline
VectorS<T>::VectorS()
    : size(0)
{
}

template <class T>
inline
VectorS<T>::VectorS(const VectorS<T>& v)
    : size(v.size)
{
    data = v.data;
}

template <class T>
inline
VectorS<T>::VectorS(Size s)
    : size(s)
{
}

template <class T>
inline
VectorS<T>::VectorS(Size s, T v)
    : size(s)
{
    if (v==0) { return; }
    for (Index i = 0; i < s; ++i) { data.push_back(Exp(i,v)); }
}

template <class T>
inline
VectorS<T>::~VectorS()
{
}

template <class T>
inline
void
VectorS<T>::normalise()
{
    if (data.empty()) { return; }
    Iter it = data.begin();
    T gcd = it->v;
    if (gcd == 1) return;
    ++it;
    while (it != data.end()) {
        euclidean(gcd, it->v, gcd);
        if (gcd == 1) return;
        ++it;
    }
    if (gcd != 1) { div(gcd); }
}

template <class T>
inline
void
VectorS<T>::print()
{
    for (Iter it = data.begin(); it != data.end(); ++it) {
        std::cout << " (" << it->i << "," << it->v << ")";
    }
}

} // namespace 4ti2

#endif
