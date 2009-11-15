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

#ifndef _4ti2_qsolve__VectorArrayT_
#define _4ti2_qsolve__VectorArrayT_

#include "qsolve/Index.h"
#include "qsolve/Size.h"
#include "qsolve/Vector.h"
#include "qsolve/IndexSetR.h"
#include <vector>

namespace _4ti2_
{

template <class T>
class VectorArrayT
{
public:
    typedef T DataType;

    VectorArrayT();
    VectorArrayT(Size m, Size n);
    VectorArrayT(const VectorArrayT<T>& vs);
    VectorArrayT<T>& operator=(const VectorArrayT<T>& vs);
    void init(Size m, Size n);
    ~VectorArrayT();


    const VectorR<T>& operator[](Index) const;
    VectorR<T>& operator[](Index);

    static void dot(const VectorArrayT<T>& vs1,
                    const VectorR<T>& v2,
                    VectorR<T>& v);
    static void dot(const VectorArrayT<T>& vs1,
                    const VectorArrayT<T>& vs2,
                    VectorArrayT<T>& vs);
    void transpose(const VectorArrayT<T>& vs1);

    void assignT(T v);
    template <class RowSet, class ColSet>
    void assignT(T v, const RowSet& rows, const ColSet& cols);

    template <class RowSet1, class ColSet1>
    void assign(const VectorArrayT<T>& vs1, const RowSet1& rows1, const ColSet1& cols1);
    template <class RowSet1, class ColSet1, class RowSet0, class ColSet0>
    void assign(const VectorArrayT<T>& vs1, const RowSet1& rows1, const ColSet1& cols1,
                                                const RowSet0& rows0, const ColSet0& cols0);

    template <class RowSet1, class ColSet1>
    void assign_trans(const VectorArrayT<T>& vs1, const RowSet1& rows1, const ColSet1& cols1);
    template <class RowSet1, class ColSet1, class RowSet0, class ColSet0>
    void assign_trans(const VectorArrayT<T>& vs1, const RowSet1& rows1, const ColSet1& cols1,
                                                const RowSet0& rows0, const ColSet0& cols0);

    void transfer(VectorArrayT<T>& vs1, Index start1, Index end1, Index pos);

    void mul(T m);
 
    bool is_index_zero(Index index) const;

    void insert(const VectorR<T>& v);
    void insert(const VectorR<T>& v, Index i);
    void insert(const VectorArrayT<T>& vs);
    void insert(const VectorArrayT<T>& vs, Index i);
    void remove(Index i);
    void remove(Index start, Index end);
    void clear();

    void normalise();

    void swap_rows(Index i1, Index i2);
    void swap_vectors(Index i1, Index i2);

    void sort();

    Size get_number() const;
    Size get_size() const;

    bool operator==(const VectorArrayT<T>&);
    bool operator!=(const VectorArrayT<T>&);

protected:
    class Row : public VectorR<T>
    {
    public:
        Row() : VectorR<T>() {}
        Row(Size s) : VectorR<T>(s) {}
        Row(Size s, T v) : VectorR<T>(s,v) {}
        Row(const Row& r) : VectorR<T>(r.start, r.end) {}
        Row& operator=(const Row& r) 
        { VectorR<T>::start = r.start; VectorR<T>::end = r.end; return *this; }
        void free_memory() { VectorR<T>::free_memory(); }
    };

    std::vector<Row> vectors;
    Size number;
    Size size;
};

template <class T> inline
const VectorR<T>&
VectorArrayT<T>::operator[](Index index) const
{
    assert(index >= 0 && index < number);
    return vectors[index];
}

template <class T> inline
VectorR<T>&
VectorArrayT<T>::operator[](Index index)
{
    assert(index >= 0 && index < number);
    return vectors[index];
}

template <class T> inline
Size
VectorArrayT<T>::get_number() const
{
    return number;
}

template <class T> inline
Size
VectorArrayT<T>::get_size() const
{
    return size;
}

template <class T> inline
bool
VectorArrayT<T>::operator==(const VectorArrayT<T>& vs2)
{
    if (get_size() != vs2.get_size()) { return false; }
    if (get_number() != vs2.get_number()) {return false;}
    for (int i = 0; i < get_number(); ++i)
    {
        if (vectors[i] != vs2[i]) { return false; }
    }
    return true;
}

template <class T> inline
bool
VectorArrayT<T>::operator!=(const VectorArrayT<T>& vs2)
{
    return !operator==(vs2);
}

template <class T>
void
VectorArrayT<T>::assignT(T v)
{
    assignT(v, IndexSetR(0,number), IndexSetR(0,size));
}

template <class T>
template <class RowSet, class ColSet>
inline
void
VectorArrayT<T>::assignT(T v, const RowSet& rows, const ColSet& cols)
{
    for (typename RowSet::Iter i = rows.begin(); i != rows.end(); ++i) {
        vectors[i].assignT(v, cols);
    }
}

template <class T>
template <class RowSet1, class ColSet1>
inline
void
VectorArrayT<T>::assign(const VectorArrayT<T>& vs1, const RowSet1& rows1, const ColSet1& cols1)
{
    assign(vs1, rows1, cols1, IndexSetR(0,number), IndexSetR(0,size));
}

template <class T>
template <class RowSet1, class ColSet1, class RowSet0, class ColSet0>
void
VectorArrayT<T>::assign(  const VectorArrayT<T>& vs1,
                            const RowSet1& rows1, const ColSet1& cols1,
                            const RowSet0& rows0, const ColSet0& cols0)
{
    assert(rows1.count() == rows0.count() && cols1.count() == cols0.count());
    typename RowSet0::Iter i = rows0.begin();
    for (typename RowSet1::Iter j = rows1.begin(); j != rows1.end(); ++j) {
        vectors[i].assign(vs1[j], cols1, cols0);
        ++i;
    }
}

template <class T>
template <class RowSet1, class ColSet1>
inline
void
VectorArrayT<T>::assign_trans(const VectorArrayT<T>& vs1, const RowSet1& rows1, const ColSet1& cols1)
{
    assign_trans(vs1, rows1, cols1, IndexSetR(0,number), IndexSetR(0,size));
}

template <class T>
template <class RowSet1, class ColSet1, class RowSet0, class ColSet0>
void
VectorArrayT<T>::assign_trans(const VectorArrayT<T>& vs1,
                            const RowSet1& rows1, const ColSet1& cols1,
                            const RowSet0& rows0, const ColSet0& cols0)
{
    assert(rows1.count() == cols0.count() && cols1.count() == rows0.count());
    typename ColSet0::Iter j0 = cols0.begin();
    for (typename RowSet1::Iter i1 = rows1.begin(); i1 != rows1.end(); ++i1) {
        typename RowSet0::Iter i0 = rows0.begin();
        for (typename ColSet1::Iter j1 = cols1.begin(); j1 != cols1.end(); ++j1) {
            vectors[i0].set(j0, vs1[i1][j1]);
            ++i0;
        }
        ++j0;
    }
}

template <class T>
void
VectorArrayT<T>::transfer(VectorArrayT<T>& vs1, Index start1, Index end1, Index pos)
{
    vectors.insert(vectors.begin()+pos, vs1.vectors.begin()+start1, vs1.vectors.begin()+end1);
    vs1.vectors.erase(vs1.vectors.begin()+start1, vs1.vectors.begin()+end1);
    vs1.number -= end1-start1;
    number += end1-start1;
}

template <class T>
VectorArrayT<T>::VectorArrayT()
    : vectors(0), number(0), size(0)
{
}

template <class T>
VectorArrayT<T>::VectorArrayT(Size _number, Size _size)
{
    number = _number;
    size = _size;
    for (Index i = 0; i < number; ++i) {
        vectors.push_back(Row(size));
    }
}

template <class T>
VectorArrayT<T>::VectorArrayT(const VectorArrayT<T>& vs)
{
    number = vs.number;
    size = vs.size;
    for (Index i = 0; i < number; i++) {
        vectors.push_back(Row(size));
        vectors.back().assign(vs[i]);
    }
}

template <class T>
VectorArrayT<T>&
VectorArrayT<T>::operator=(const VectorArrayT<T>& vs)
{
    for (Index i = 0; i < number; i++) { vectors[i].free_memory(); }
    vectors.clear();
    number = vs.number;
    size = vs.size;
    for (Index i = 0; i < number; i++) {
        vectors.push_back(Row(size));
        vectors.back.assign(vs[i]);
    }
    return *this;
}

template <class T>
VectorArrayT<T>::~VectorArrayT()
{
    clear();
}

template <class T>
void
VectorArrayT<T>::insert(const VectorR<T>& v)
{
    assert(v.get_size() == size);
    ++number;
    vectors.push_back(Row(size));
    vectors.back().assign(v);
}

template <class T>
void
VectorArrayT<T>::insert(const VectorR<T>& v, Index i)
{
    assert(v.get_size() == size);
    ++number;
    vectors.insert(vectors.begin()+i, Row(v))->assign(v);
}

template <class T>
void
VectorArrayT<T>::insert(const VectorArrayT<T>& vs)
{
    assert(vs.get_size() == size);
    vectors.reserve(size + vs.get_number());
    for (Index i = 0; i < vs.get_number(); ++i) { insert(vs[i]); }
}

template <class T>
void
VectorArrayT<T>::insert(const VectorArrayT<T>& vs, Index i)
{
    assert(vs.get_size() == size);
    vectors.reserve(size + vs.get_number());
    for (Index i = 0; i < vs.get_number(); ++i) { insert(vs[i], i); }
}

template <class T>
void
VectorArrayT<T>::clear()
{
    for (Index i = 0; i < number; ++i) { vectors[i].free_memory(); }
    vectors.clear();
    number = 0;
}

template <class T>
void
VectorArrayT<T>::remove(Index index)
{
    assert(index < number);
    vectors[index].free_memory();
    --number;
    vectors.erase(vectors.begin()+index);
}

template <class T>
void
VectorArrayT<T>::remove(Index start, Index end)
{
    assert(start <= end && end <= number);
    for (Index i = start; i < end; ++i) { vectors[i].free_memory(); }
    number -= end-start;
    vectors.erase(vectors.begin()+start, vectors.begin()+end);
    assert((Size) vectors.size() == number);
}

// Matrix multiplication.  Computes AxB^T.
template <class T>
void
VectorArrayT<T>::dot(
                const VectorArrayT<T>& vs1,
                const VectorR<T>& v,
                VectorR<T>& r)
{
    assert(vs1.size == v.get_size());
    assert(r.get_size() == vs1.number);

    T value;
    for (Index i = 0; i < vs1.number; ++i) {
        v.dot(vs1[i], value);
        r.set(i, value);
    }
}

// Matrix multiplication.  Computes BxA^T.
template <class T>
void
VectorArrayT<T>::dot(
                const VectorArrayT<T>& vs1,
                const VectorArrayT<T>& vs2,
                VectorArrayT<T>& rs)
{
    assert(vs1.size == vs2.size);
    assert(rs.size == vs1.number);
    assert(rs.number == vs2.number);

    for (Index i = 0; i < vs2.number; ++i) {
        dot(vs1, vs2[i], rs[i]);
    }
}

template <class T>
void
VectorArrayT<T>::transpose(const VectorArrayT<T>& vs)
{
    assert(get_number() == vs.get_size());
    // TODO: assert(get_size() == vs.get_number());

    for (Index i = 0; i < vs.get_number(); ++i) {
        for (Index j = 0; j < vs.get_size(); ++j) {
            vectors[j].set(i,vs[i][j]);
        }
    }
}

template <class T>
void
VectorArrayT<T>::mul(T m)
{
    for (Index i = 0; i < number; ++i) { vectors[i].muleq(m); }
}

template <class T>
bool
VectorArrayT<T>::is_index_zero(Index index) const
{
    for (Index i = 0; i < number; ++i) {
        if (vectors[i][index] != 0) { return false; }
    }
    return true;
}

template <class T>
void
VectorArrayT<T>::sort()
{
    std::sort(vectors.begin(), vectors.end());
}

template <class T>
void
VectorArrayT<T>::normalise()
{
    for (Index i = 0; i < number; ++i) { vectors[i].normalise(); }
}

template <class T>
inline void
VectorArrayT<T>::swap_rows(Index i1, Index i2)
{
    swap_vectors(i1, i2);
}

template <class T>
void
VectorArrayT<T>::swap_vectors(Index i1, Index i2)
{
    assert(i1 >= 0 && i2 >= 0);
    assert(i1 < number && i2 < number);
    if (i1 == i2) return;
    vectors[i1].swap(vectors[i2]);
}

template <class T>
void
VectorArrayT<T>::init(Size m, Size n)
{
    if (n != size) { clear(); size = n; }

    if (m == number) { return; }
    else if (m < number) {
        for (Index i = m; i < number; ++i) { vectors[i].free_memory(); }
        vectors.resize(m);
        number = m;
    }
    else { // m > number
        for (Index i = number; i < m; ++i) {
            vectors.push_back(Row(size));
        }
        number = m;
    }
}

} // namespace _4ti2_

#endif
