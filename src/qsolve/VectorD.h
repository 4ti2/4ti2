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

#ifndef _4ti2_qsolve__VectorD_
#define _4ti2_qsolve__VectorD_

#include "qsolve/Index.h"
#include "qsolve/Size.h"
#include "qsolve/Euclidean.h"

#include <cassert>
#include <algorithm>

namespace _4ti2_
{

// TODO: Should we pass gmp by reference?
template <class T>
class VectorDBase
{
public:
    const Size get_size() const;

    const T& operator[](Index) const;
    T& operator[](Index);
    void set(Index, T);

    void add(const VectorDBase<T>&, T, const VectorDBase<T>&, T);
    void add(const VectorDBase<T>&, const VectorDBase<T>&);
    void sub(const VectorDBase<T>&, const VectorDBase<T>&);
    void mul(const VectorDBase<T>&, T);
    void div(const VectorDBase<T>&, T);

    void addeq(const VectorDBase<T>&);
    void addeq(const VectorDBase<T>&, T);
    void subeq(const VectorDBase<T>&);
    void muleq(T);
    void diveq(T);

    bool same(const VectorDBase<T>&) const;
    bool operator==(const VectorDBase<T>&) const;
    bool operator!=(const VectorDBase<T>&) const;
    bool operator< (const VectorDBase<T>&) const;

    void swap(VectorDBase<T>&);

    T dot(const VectorDBase<T>&) const;
    void dot(const VectorDBase<T>&, T&) const;

    void assignT(T v);
    template <class IndexSet>
    void assignT(T v, const IndexSet& is);

    void assign(const VectorDBase<T>& v1);
    template <class IndexSet>
    void assign(const VectorDBase<T>& v1, const IndexSet& is);
    template <class IndexSet1, class IndexSet2>
    void assign(const VectorDBase<T>& v1, const IndexSet1& is1,
                const IndexSet2& is2);

    void normalise();

protected:
    T* start;
    T* end;

protected:
    VectorDBase();
    VectorDBase(Size s);
    VectorDBase(Size s, T v);
    VectorDBase(T*, T*);
    VectorDBase(const VectorDBase<T>&); // Not implemented
    VectorDBase<T>& operator=(const VectorDBase<T>&); // Not implemented.

    void allocate_memory(Size s);
    void free_memory();
};

template <class T>
class VectorD : public VectorDBase<T>
{
public:
    VectorD();
    explicit VectorD(Size s);
    VectorD(const VectorD<T>&);
    VectorD(Size s, T value);
    VectorD<T>& operator=(const VectorD<T>&);
    ~VectorD();
};

template <class T>
class VectorDLW : public VectorDBase<T>
{
public:
    VectorDLW();
    VectorDLW(const VectorDBase<T>&);
    VectorDLW(T* start, T* end);
    VectorDLW<T>& operator=(const VectorDBase<T>&);
};

} // namespace 4ti2

#include "VectorD.hpp"

#endif
