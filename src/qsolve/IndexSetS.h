/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Peter Malkin.

This program indices free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program indices distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 
*/

#ifndef _4ti2_qsolve__IndexSetS_
#define _4ti2_qsolve__IndexSetS_

#include "qsolve/Size.h"
#include "qsolve/Index.h"

#include <cassert>
#include <algorithm>

// TODO: Remove this dependency on iostream.
#include <iostream>

#include "IndexSetR.h"

namespace _4ti2_
{

class IndexSetS
{
public:
    explicit IndexSetS(Size _size, bool v = false);
    IndexSetS(const IndexSetS&);
    IndexSetS& operator=(const IndexSetS&);
    ~IndexSetS();

    static bool set_subset(const IndexSetS& b1, const IndexSetS& b2);
    static bool set_subset(const IndexSetS& b1, const IndexSetS& b2, const IndexSetS& mask);
    static bool set_disjoint(const IndexSetS& b1, const IndexSetS& b2);
    static bool set_disjoint(const IndexSetS& b1, const IndexSetS& b2, const IndexSetS& mask);

    static void set_union(const IndexSetS& b1, const IndexSetS& b2, IndexSetS& b3);
    static void set_intersection(const IndexSetS& b1, const IndexSetS& b2, IndexSetS& b3);
    static void set_difference(const IndexSetS& b1, const IndexSetS& b2, IndexSetS& b3);
    static void set_complement(const IndexSetS& b1, IndexSetS& b2);

    void set_union(const IndexSetS& b);
    void set_intersection(const IndexSetS& b);
    void set_difference(const IndexSetS& b);
    void set_complement();

    static void extend(const IndexSetS& b1, IndexSetS& b2);
    static void shrink(const IndexSetS& b1, IndexSetS& b2);

    bool singleton() const;
    bool less_than_equal(Size s) const;
    bool empty() const;

    Size count() const;

    bool operator[](Index index) const;

    Size get_size() const;
    void resize(Size s);

    void set(Index index);
    void unset(Index index);
    void flip(Index index);
    void zero();
    void one();

    static void swap(IndexSetS& b1, IndexSetS& b2);

    friend bool operator==(const IndexSetS& b1, const IndexSetS& b2);

    class Iter
    {
    public:
        operator Index() const { return *i; }
        Index operator*() const { return *i; }
        bool operator!=(Iter it) const { return i != it.i; }
        Iter& operator=(Iter it) { i = it.i; return *this; }
        void operator++() { ++i; }
    private:
        Iter(Index* _i) : i(_i) {}
        Index* i;
        friend class IndexSetS;
    };

    Iter begin() const { return indices; }
    Iter end() const { return finish; }

protected:
    IndexSetS();

    Index *indices;
    Index *finish;
    Size size;
};

// TODO: Should this return true if and only if there is exactly zero or one elements?
inline
bool
IndexSetS::singleton() const
{
    return (finish-indices <= 1);
}

inline
bool
IndexSetS::less_than_equal(Size s) const
{
    return (finish-indices <= s);
}

inline
bool
IndexSetS::empty() const
{
    return (finish == indices);
}

#if 0
inline
void
IndexSetS::set(Index index)
{
    assert(index >= 0 && index < size);
    Index* i = indices;
    while (i != finish && *i < index)  { ++i; }
    if (i != finish && *i == index) { return; }
    for (Index* j = finish; j != i; --j) { *j = *(j-1); }
    *i = index;
    ++finish;
}
#endif

inline
void
IndexSetS::set(Index index)
{
    assert(index >= 0 && index < size);
    Index* i = finish-1;
    while (i >= indices && *i > index)  { --i; }
    if (i >= indices && *i == index) { return; }
    ++i;
    for (Index* j = finish; j != i; --j) { *j = *(j-1); }
    *i = index;
    ++finish;
}

inline
void
IndexSetS::unset(Index index)
{
    assert(index >= 0 && index < size);
    Index* i = indices;
    while (i != finish && *i < index)  { ++i; }
    if (i == finish || *i != index) { return; }
    --finish;
    while (i != finish) { *i = *(i+1); ++i; }
}

inline
void
IndexSetS::flip(Index index)
{
    assert(index >= 0 && index < size);
    Index* i = indices;
    while (i != finish && *i < index)  { ++i; }
    if (i == finish || *i != index) {
        for (Index* j = finish; j != i; --j) { *j = *(j-1); }
        *i = index;
        ++finish;
    }
    else {
        --finish;
        while (i != finish) { *i = *(i+1); ++i; }
    }
}

inline
void
IndexSetS::zero()
{
    finish = indices;
}

inline
void
IndexSetS::one()
{
    finish = indices+size;
    for (Index* i = indices; i != finish; ++i) { *i = (i-indices); }
}

inline
void
IndexSetS::swap(IndexSetS& b1, IndexSetS& b2)
{
    std::swap(b1.indices, b2.indices);
    std::swap(b1.finish, b2.finish);
}

inline
bool
IndexSetS::operator[](Index index) const
{
    assert(index >= 0 && index <= size);
    Index* i = indices;
    while (i != finish && *i < index) { ++i; }
    if (i != finish && *i == index) { return true; }
    return false;
}

inline
Size
IndexSetS::get_size() const
{
    return size;
}

inline
bool
operator==(const IndexSetS& b1, const IndexSetS& b2)
{
    if (b1.count() != b2.count()) { return false; }
    Index* i2 = b2.indices;
    for (Index* i1 = b1.indices; i1 != b1.finish; ++i1) {
        if (*i1 != *i2) { return false; }
        ++i2;
    }
    return true;
}

inline
bool
IndexSetS::set_subset(const IndexSetS& b1, const IndexSetS& b2)
{
    Index* i1 = b1.indices;
    Index* i2 = b2.indices;
    while (i1 != b1.finish) {
        while (i2 != b2.finish && *i2 < *i1) { ++i2; }
        if (i2 == b2.finish || *i2 != *i1) { return false; }
        ++i1;
    }
    return true;
}

inline
bool
IndexSetS::set_subset(const IndexSetS& b1, const IndexSetS& b2, const IndexSetS& mask)
{
    IndexSetS b(b1.size);
    IndexSetS::set_difference(b1, mask, b);
    return IndexSetS::set_subset(b, b2);
}

inline
bool
IndexSetS::set_disjoint(const IndexSetS& b1, const IndexSetS& b2)
{
    Index* i1 = b1.indices;
    Index* i2 = b2.indices;
    while (i1 != b1.finish && i2 != b2.finish) {
        if (*i1 < *i2) { ++i1; }
        else if (*i1 > *i2) { ++i2; }
        else { return false; }
    }
    return true;
}

inline
bool
IndexSetS::set_disjoint(const IndexSetS& b1, const IndexSetS& b2, const IndexSetS& mask)
{
    IndexSetS b(b1.size);
    IndexSetS::set_intersection(b1, mask, b);
    return IndexSetS::set_disjoint(b, b2);
}

inline
void
IndexSetS::set_union(
                const IndexSetS& b1,
                const IndexSetS& b2,
                IndexSetS& b)
{
    Index* i1 = b1.indices;
    Index* i2 = b2.indices;
    Index* i = b.indices;
    while (i1 != b1.finish && i2 != b2.finish) {
        if (*i1 < *i2) { *i = *i1; ++i; ++i1; }
        else if (*i1 > *i2) { *i = *i2; ++i; ++i2; }
        else { *i = *i1; ++i; ++i1; ++i2; }
    }
    while (i1 != b1.finish) { *i = *i1; ++i; ++i1; }
    while (i2 != b2.finish) { *i = *i2; ++i; ++i2; }
    b.finish = i;
}

inline
void
IndexSetS::set_intersection(
                const IndexSetS& b1,
                const IndexSetS& b2,
                IndexSetS& b)
{
    Index* i1 = b1.indices;
    Index* i2 = b2.indices;
    Index* i = b.indices;
    while (i1 != b1.finish && i2 != b2.finish) {
        if (*i1 < *i2) { ++i1; }
        else if (*i1 > *i2) { ++i2; }
        else { *i = *i1; ++i; ++i1; ++i2; }
    }
    b.finish = i;
}

inline
void
IndexSetS::set_difference(
                const IndexSetS& b1,
                const IndexSetS& b2,
                IndexSetS& b)
{
    Index* i1 = b1.indices;
    Index* i2 = b2.indices;
    Index* i = b.indices;
    while (i1 != b1.finish && i2 != b2.finish) {
        if (*i1 < *i2) { *i = *i1; ++i; ++i1; }
        else if (*i1 > *i2) { ++i2; }
        else { ++i1; ++i2; }
    }
    while (i1 != b1.finish) { *i = *i1; ++i; ++i1; }
    b.finish = i;
}

inline
void
IndexSetS::set_complement(
                const IndexSetS& b1,
                IndexSetS& b)
{
    Index* i = b.indices;
    Index j = 0;
    for (Index* i1 = b1.indices; i1 != b1.finish; ++i1) {
        while (j < *i1) { *i = j; ++i; ++j; }
        ++j;
    }
    while (j < b.size) { *i = j; ++j; ++i; }
    b.finish = i;
}

inline
void
IndexSetS::set_union(const IndexSetS& b1)
{
    IndexSetS b(size);
    IndexSetS::set_union(*this, b1, b);
    IndexSetS::swap(*this, b);
}

inline
void
IndexSetS::set_intersection(const IndexSetS& b2)
{
    Index* i1 = indices;
    Index* i2 = b2.indices;
    Index* i = indices;
    while (i1 != finish && i2 != b2.finish) {
        if (*i1 < *i2) { ++i1; }
        else if (*i1 > *i2) { ++i2; }
        else { *i = *i1; ++i; ++i1; ++i2; }
    }
    finish = i;
}

inline
void
IndexSetS::set_difference(const IndexSetS& b2)
{
    Index* i1 = indices;
    Index* i2 = b2.indices;
    Index* i = indices;
    while (i1 != finish && i2 != b2.finish) {
        if (*i1 < *i2) { *i = *i1; ++i; ++i1; }
        else if (*i1 > *i2) { ++i2; }
        else { ++i1; ++i2; }
    }
    while (i1 != finish) { *i = *i1; ++i; ++i1; }
    finish = i;
}

inline
void
IndexSetS::set_complement()
{
    IndexSetS b(size);
    IndexSetS::set_complement(*this, b);
    IndexSetS::swap(*this, b);
}

inline
void
IndexSetS::extend(const IndexSetS& b1, IndexSetS& b2)
{
    b2 = b1;
}

inline
void
IndexSetS::shrink(const IndexSetS& b1, IndexSetS& b2)
{
    Index* i1 = b1.indices;
    Index* i2 = b2.indices;
    while (i1 != b1.finish && *i1 < b2.size) {
        *i2 = *i1;
        ++i1; ++i2;
    }
    b2.finish = i2;
}

inline
Size
IndexSetS::count() const
{
    return finish-indices;
}

inline
IndexSetS::IndexSetS(Size _size, bool v)
        : size(_size) 
{
    indices = new Index[size];
    finish = indices;
    if (v) { one(); }
}

inline
IndexSetS::IndexSetS(const IndexSetS& b1)
        : size(b1.size) 
{
    indices = new Index[size];
    Index* i = indices;
    for (Index* i1 = b1.indices; i1 != b1.finish; ++i1) { *i = *i1; ++i; }
    finish = i;
}

inline
IndexSetS&
IndexSetS::operator=(const IndexSetS& b1)
{
    assert(b1.count() <= size);
    Index* i = indices;
    for (Index* i1 = b1.indices; i1 != b1.finish; ++i1) { *i = *i1; ++i; }
    finish = i;
    return *this;
}

inline
IndexSetS::~IndexSetS()
{
    delete [] indices;
}

inline
void
print(const IndexSetS& is)
{
    for (IndexSetS::Iter i = is.begin(); i != is.end(); ++i) {
        std::cout << *i << " ";
    }
    std::cout << "\n";
}

} // namespace _4ti2_

#endif
