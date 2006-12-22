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

#ifndef _4ti2__VectorArray_
#define _4ti2__VectorArray_

#include "Vector.h"
#include <vector>
#include "Permutation.h"
#include "Index.h"
#include "Size.h"

namespace _4ti2_
{

class VectorArray
{
public:
    VectorArray();
    VectorArray(Size number, Size size);
    VectorArray(Size number, Size size, IntegerType v);
    VectorArray(const VectorArray& vs);
    VectorArray& operator=(const VectorArray& vs);
    ~VectorArray();

    void renumber(Size m);
    void renumber(Size m, const Vector& v);

    const Vector& operator[](Index) const;
    Vector& operator[](Index);

    static void dot(const VectorArray& vs1,
                    const Vector& v2,
                    Vector& v);
    static void dot(const VectorArray& vs1,
                    const VectorArray& vs2,
                    VectorArray& vs);
    static void transpose(
                    const VectorArray& vs1,
                    VectorArray& vs);
    static void concat(
                    const VectorArray& vs1,
                    const VectorArray& vs2,
                    VectorArray& vs);
    template <class IndexSet>
    static void concat(
                    const VectorArray& vs1,
                    const VectorArray& vs2,
                    const IndexSet& s,
                    VectorArray& vs);
    static void split(
                    const VectorArray& vs,
                    VectorArray& vs1,
                    VectorArray& vs2);
    template <class IndexSet>
    static void split(
                    const VectorArray& vs,
                    const IndexSet& s,
                    VectorArray& vs1,
                    VectorArray& vs2);
    static void project(
                    const VectorArray& vs1,
                    Index start,
                    Index end,
                    VectorArray& vs);
    template <class IndexSet>
    static void project(
                    const VectorArray& vs1,
                    const IndexSet& proj,
                    VectorArray& vs);
    template <class IndexSet>
    static void lift(
                    const VectorArray& vs1,
                    const IndexSet& lift,
                    VectorArray& vs);
    static void lift(
                    const VectorArray& vs1,
                    Index start,
                    Index end,
                    VectorArray& vs);

    void project(Index start, Index end);
    template <class IndexSet>
    void project(const IndexSet& proj);
    void lift(Index start, Index end);
    template <class IndexSet>
    void lift(const IndexSet& lift);

    static void transfer(
                    VectorArray& vs1,
                    VectorArray& vs2);
    static void transfer(
                    VectorArray& vs1,
                    Index start,
                    Index end,
                    VectorArray& vs2,
                    Index pos);
    
    void mul(IntegerType m);
    
    bool is_index_zero(Index index) const;

    void insert(const Vector& v);
    void insert(const Vector& v, Index i);
    void insert(const VectorArray& vs);
    void insert(const VectorArray& vs, Index i);
    void remove(Index i);
    void remove(Index start, Index end);
    void clear();

    void permute(const Permutation& p);

    void normalise();

    void swap_vectors(Index i1, Index i2);
    void swap_indices(Index i1, Index i2);

    void sort();

    Size get_number() const;
    Size get_size() const;

    friend bool operator==(const VectorArray&, const VectorArray&);

protected:
    void insert(Vector* v);

    std::vector<Vector*> vectors;
    Size number;
    Size size;
};

inline
const Vector&
VectorArray::operator[](Index index) const
{
    assert(index >= 0 && index < number);
    return *vectors[index];
}

inline
Vector&
VectorArray::operator[](Index index)
{
    assert(index >= 0 && index < number);
    return *vectors[index];
}

inline
Size
VectorArray::get_number() const
{
    return number;
}

inline
Size
VectorArray::get_size() const
{
    return size;
}

inline
bool
operator==(const VectorArray& vs1, const VectorArray& vs2)
{
    if (vs1.get_size() != vs2.get_size()) { return false; }
    if (vs1.get_number() != vs2.get_number()) {return false;}
    for (int i = 0; i < vs1.get_number(); ++i)
    {
        if (vs1[i] != vs2[i]) { return false; }
    }
    return true;
}

template <class IndexSet>
void
VectorArray::concat(
                const VectorArray& vs1,
                const VectorArray& vs2,
                const IndexSet& s,
                VectorArray& vs)
{
    assert(vs1.get_number() == vs2.get_number());
    assert(vs1.get_number() == vs.get_number());
    for (Index i = 0; i < vs1.get_number(); ++i)
    {
        Vector::concat(vs1[i], vs2[i], s, vs[i]);
    }
}

template <class IndexSet>
void
VectorArray::split(
                const VectorArray& vs,
                const IndexSet& s,
                VectorArray& vs1,
                VectorArray& vs2)
{
    assert(vs1.get_number() == vs2.get_number());
    assert(vs1.get_number() == vs.get_number());
    for (Index i = 0; i < vs1.get_number(); ++i)
    {
        Vector::split(vs[i], s, vs1[i], vs2[i]);
    }
}

template <class IndexSet>
void
VectorArray::project(
                const VectorArray& vs1,
                const IndexSet& proj,
                VectorArray& vs)
{
    assert(vs1.get_number() == vs.get_number());
    assert(proj.get_size() == vs1.get_size());
    assert(vs.get_size() == proj.count());
    for (Index i = 0; i < vs1.get_number(); ++i)
    {
        Vector::project(vs1[i], proj, vs[i]);
    }
}

template <class IndexSet>
void
VectorArray::lift(
                const VectorArray& vs1,
                const IndexSet& proj,
                VectorArray& vs)
{
    assert(vs1.get_number() == vs.get_number());
    assert(proj.get_size() == vs.get_size());
    assert(vs1.get_size() == proj.count());
    for (Index i = 0; i < vs1.get_number(); ++i)
    {
        Vector::lift(vs1[i], proj, vs[i]);
    }
}

inline
void
VectorArray::transfer(
                    VectorArray& vs1,
                    VectorArray& vs2)
{
    transfer(vs1, 0, vs1.get_number(), vs2, vs2.get_number());
}

} // namespace _4ti2_

#endif
