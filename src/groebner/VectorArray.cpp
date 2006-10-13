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

#include "VectorArray.h"
#include <algorithm>

using namespace _4ti2_;

VectorArray::VectorArray()
    : vectors(0), number(0), size(0)
{
}

VectorArray::VectorArray(Size _number, Size _size)
{
    number = _number;
    size = _size;
    for (Index i = 0; i < number; ++i)
    {
        vectors.push_back(new Vector(size));
    }
}

VectorArray::VectorArray(Size _number, Size _size, IntegerType v)
{
    number = _number;
    size = _size;
    for (Index i = 0; i < number; ++i)
    {
        vectors.push_back(new Vector(size, v));
    }
}

VectorArray::VectorArray(const VectorArray& vs)
{
    number = vs.number;
    size = vs.size;
    for (Index i = 0; i < number; i++)
    {
        vectors.push_back(new Vector(vs[i]));
    }
}

VectorArray&
VectorArray::operator=(const VectorArray& vs)
{
    for (Index i = 0; i < number; i++) { delete vectors[i]; }
    vectors.clear();
    number = vs.number;
    size = vs.size;
    for (Index i = 0; i < number; i++)
    {
        vectors.push_back(new Vector(vs[i]));
    }
    return *this;
}

VectorArray::~VectorArray()
{
    clear();
}

void
VectorArray::insert(const Vector& v)
{
    assert(v.get_size() == size);
    ++number;
    vectors.push_back(new Vector(v));
}

void
VectorArray::insert(Vector* v)
{
    assert(v->get_size() == size);
    ++number;
    vectors.push_back(v);
}

void
VectorArray::insert(const Vector& v, Index i)
{
    assert(v.get_size() == size);
    ++number;
    vectors.insert(vectors.begin()+i, new Vector(v));
}

void
VectorArray::insert(const VectorArray& vs)
{
    assert(vs.get_size() == size);
    vectors.reserve(size + vs.get_number());
    for (Index i = 0; i < vs.get_number(); ++i) { insert(vs[i]); }
}

void
VectorArray::insert(const VectorArray& vs, Index i)
{
    assert(vs.get_size() == size);
    vectors.reserve(size + vs.get_number());
    for (Index i = 0; i < vs.get_number(); ++i) { insert(vs[i], i); }
}

void
VectorArray::clear()
{
    for (Index i = 0; i < number; ++i) { delete vectors[i]; }
    vectors.clear();
    number = 0;
}

void
VectorArray::remove(Index index)
{
    assert(index < number);
    delete vectors[index];
    --number;
    vectors.erase(vectors.begin()+index);
}

void
VectorArray::remove(Index start, Index end)
{
    assert(start <= end && end <= number);
    for (Index i = start; i < end; ++i) { delete vectors[i]; }
    number -= end-start;
    vectors.erase(vectors.begin()+start, vectors.begin()+end);
    assert((Size) vectors.size() == number);
}

// Matrix multiplication.  Computes AxB^T.
void
VectorArray::dot(
                const VectorArray& vs1,
                const Vector& v,
                Vector& r)
{
    assert(vs1.size == v.get_size());
    assert(r.get_size() == vs1.number);

    for (Index i = 0; i < vs1.number; ++i)
    {
        Vector::dot(vs1[i], v, r[i]);
    }
}

// Matrix multiplication.  Computes BxA^T.
void
VectorArray::dot(
                const VectorArray& vs1,
                const VectorArray& vs2,
                VectorArray& rs)
{
    assert(vs1.size == vs2.size);
    assert(rs.size == vs1.number);
    assert(rs.number == vs2.number);

    for (Index i = 0; i < vs2.number; ++i)
    {
        VectorArray::dot(vs1, vs2[i], rs[i]);
    }
}

void
VectorArray::transpose(
                const VectorArray& vs,
                VectorArray& t)
{
    assert(t.get_number() == vs.get_size());
    assert(t.get_size() == vs.get_number());

    for (Index i = 0; i < vs.get_number(); ++i)
    {
        for (Index j = 0; j < vs.get_size(); ++j)
        {
            t[j][i] = vs[i][j];
        }
    }
}

void
VectorArray::concat(
                const VectorArray& vs1,
                const VectorArray& vs2,
                VectorArray& vs)
{
    assert(vs1.get_number() == vs2.get_number());
    assert(vs1.get_number() == vs.get_number());
    assert(vs.get_size() == vs1.get_size()+vs2.get_size());
    for (Index i = 0; i < vs1.get_number(); ++i)
    {
        Vector::concat(vs1[i], vs2[i], vs[i]);
    }
}

void
VectorArray::split(
                const VectorArray& vs,
                VectorArray& vs1,
                VectorArray& vs2)
{
    assert(vs1.get_number() == vs2.get_number());
    assert(vs1.get_number() == vs.get_number());
    for (Index i = 0; i < vs1.get_number(); ++i)
    {
        Vector::split(vs[i], vs1[i], vs2[i]); 
    }
}

void
VectorArray::project(
                const VectorArray& vs1,
                Index start,
                Index end,
                VectorArray& vs)
{
    assert(start <= end && end <= vs1.get_size());
    assert(vs1.get_number() == vs.get_number());
    assert(vs.get_size() == end-start);
    for (Index i = 0; i < vs1.get_number(); ++i)
    {
        Vector::project(vs1[i], start, end , vs[i]);
    }
}

void
VectorArray::lift(
                const VectorArray& vs1,
                Index start,
                Index end,
                VectorArray& vs)
{
    assert(start <= end && end <= vs.get_size());
    assert(vs1.get_number() == vs.get_number());
    assert(vs1.get_size() == end-start);
    for (Index i = 0; i < vs1.get_number(); ++i)
    {
        Vector::lift(vs1[i], start, end , vs[i]);
    }
}

void
VectorArray::transfer(
                VectorArray& vs1,
                Index start,
                Index end,
                VectorArray& vs2,
                Index pos)
{
    assert(start >= 0 && start <= end);
    assert(end <= vs1.get_number());
    assert(pos >= 0 && pos <= vs2.get_number());
    vs2.vectors.insert(vs2.vectors.begin()+pos, vs1.vectors.begin()+start, vs1.vectors.begin()+end);
    vs1.vectors.erase(vs1.vectors.begin()+start, vs1.vectors.begin()+end);
    vs1.number -= end-start;
    vs2.number += end-start;
}

void
VectorArray::mul(IntegerType m)
{
    for (Index i = 0; i < number; ++i)
    {
        vectors[i]->mul(m);
    }
}

bool
VectorArray::is_index_zero(Index index) const
{
    for (Index i = 0; i < number; ++i)
    {
        if ((*vectors[i])[index] != 0) { return false; }
    }
    return true;
}

bool compare(const Vector* p1, const Vector* p2) { return *p1 < *p2; }

void
VectorArray::sort()
{
    std::sort(vectors.begin(), vectors.end(), &compare);
}

void
VectorArray::permute(const Permutation& p)
{
    for (Index i = 0; i < number; ++i) { vectors[i]->permute(p); }
}

void
VectorArray::normalise()
{
    for (Index i = 0; i < number; ++i) { vectors[i]->normalise(); }
}

void
VectorArray::swap_indices(Index i1, Index i2)
{
    assert(i1 >= 0 && i2 >= 0);
    assert(i1 < size && i2 < size);
    if (i1 == i2) return;
    for (Index i = 0; i < number; ++i)
    {
        IntegerType temp = (*vectors[i])[i1];
        (*vectors[i])[i1] = (*vectors[i])[i2];
        (*vectors[i])[i2] = temp;
    }
}

void
VectorArray::swap_vectors(Index i1, Index i2)
{
    assert(i1 >= 0 && i2 >= 0);
    assert(i1 < number && i2 < number);
    if (i1 == i2) return;
    Vector* temp = vectors[i1];
    vectors[i1] = vectors[i2];
    vectors[i2] = temp;
}

void
VectorArray::renumber(Size m)
{
    renumber(m, Vector(size));
}

void
VectorArray::renumber(Size m, const Vector& v)
{
    assert(v.get_size() == size);
    if (m == number)
    {
        return;
    }
    else if (m < number)
    {
        for (Index i = m; i < number; ++i) { delete vectors[i]; }
        vectors.resize(m);
        number = m;
    }
    else // m > number
    {
        for (Index i = number; i < m; ++i)
        {
            vectors.push_back(new Vector(v));
        }
        number = m;
    }
}
