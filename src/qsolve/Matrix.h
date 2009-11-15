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

#ifndef _4ti2_qsolve__MatrixT_
#define _4ti2_qsolve__MatrixT_

#include "qsolve/Index.h"
#include "qsolve/Size.h"
#include "qsolve/Debug.h"
#include <vector>

namespace _4ti2_
{

template <class T>
class MatrixT
{
public:
    MatrixT();
    MatrixT(Size m, Size n);
    ~MatrixT();

    void init(Size m, Size n);

    const T operator()(Index i, Index j) const;

    template <class RowSet, class ColSet>
    void assign(const VectorArrayT<T>& vs1, const RowSet& rows, const ColSet& cols);
    template <class RowSet, class ColSet>
    void assign_trans(const VectorArrayT<T>& vs1, const RowSet& rows, const ColSet& cols);

    template <class RowSet, class ColSet>
    void assign(const MatrixT<T>& mat, const RowSet& rows, const ColSet& cols);
    template <class RowSet, class ColSet>
    void assign_trans(const MatrixT<T>& mat, const RowSet& rows, const ColSet& cols);

    typedef std::vector<std::pair<Index,Index> > Pivots;
    Size row_triangulate();
    template <class RowSet, class ColSet>
    Size row_triangulate(const RowSet& rows, const ColSet& cols, Pivots* pivots = 0);

    Size row_diagonalise();
    template <class RowSet, class ColSet>
    Size row_diagonalise(const RowSet& rows, const ColSet& cols, Pivots* pivots = 0);

    void row_diagonalise(const Pivots& pivots);
    template <class RowSet>
    void row_diagonalise(const Pivots& pivots, const RowSet& rows);
    void row_normalise();

    Size get_m() const;
    Size get_n() const;

    void print() const;

protected:
    T* start;
    T* end;

    Size m;
    Size n;
};

template <class T>
MatrixT<T>::MatrixT()
    : start(0), end(0), m(0), n(0)
{
}

template <class T>
MatrixT<T>::MatrixT(Size _m, Size _n)
    : m(_m), n(_n)
{
    start = new T[n*m];
    end = start + n*m;
}

template <class T>
MatrixT<T>::~MatrixT()
{
    delete [] start;
}

template <class T> inline
void
MatrixT<T>::init(Size _m, Size _n)
{
    m=_m; n=_n;
    if (m*n > end-start) { 
        delete [] start; 
        start = new T[n*m];
        end = start + n*m;
    }
}

template <class T> inline
const T
MatrixT<T>::operator()(Index i, Index j) const
{
    assert(i >= 0 && i < m);
    assert(j >= 0 && j < n);
    return *(start+i*n+j);
}

template <class T> inline
Size
MatrixT<T>::get_m() const
{
    return m;
}

template <class T> inline
Size
MatrixT<T>::get_n() const
{
    return n;
}

template <class T> template <class RowSet, class ColSet> inline
void
MatrixT<T>::assign(const VectorArrayT<T>& vs, const RowSet& rows, const ColSet& cols)
{
    T* s = start;
    for (typename RowSet::Iter i = rows.begin(); i != rows.end(); ++i) {
        const VectorR<T>& vi = vs[i];
        for (typename ColSet::Iter j = cols.begin(); j != cols.end(); ++j) {
            *s = vi[j];
            ++s;
        }
    }
}

template <class T> template <class RowSet, class ColSet> inline
void
MatrixT<T>::assign_trans(const VectorArrayT<T>& vs, const RowSet& rows, const ColSet& cols)
{
    T* s = start;
    for (typename ColSet::Iter j = cols.begin(); j != cols.end(); ++j) {
        for (typename RowSet::Iter i = rows.begin(); i != rows.end(); ++i) {
            *s = vs[i][j];
            ++s;
        }
    }
}

template <class T> template <class RowSet, class ColSet> inline
void
MatrixT<T>::assign(const MatrixT<T>& mat, const RowSet& rows, const ColSet& cols)
{
    T* s = start;
    for (typename RowSet::Iter i = rows.begin(); i != rows.end(); ++i) {
        T* mi = mat.start + i*mat.n;
        for (typename ColSet::Iter j = cols.begin(); j != cols.end(); ++j) {
            *s = *(mi+j);
            ++s;
        }
    }
}

template <class T> template <class RowSet, class ColSet> inline
void
MatrixT<T>::assign_trans(const MatrixT<T>& mat, const RowSet& rows, const ColSet& cols)
{
    T* s = start;
    for (typename ColSet::Iter j = cols.begin(); j != cols.end(); ++j) {
        T* mj = mat.start + j;
        for (typename RowSet::Iter i = rows.begin(); i != rows.end(); ++i) {
            *s = mj + i*mat.n;
            ++s;
        }
    }
}

template <class T> inline
Size
MatrixT<T>::row_triangulate()
{
    return row_triangulate(IndexSetR(0,m), IndexSetR(0,n));
}

template <class T> template <class RowSet, class ColSet>
Size
MatrixT<T>::row_triangulate(const RowSet& rows, const ColSet& cols, Pivots* pivots)
{
    std::vector<std::vector<T*> > col_entries(n);
    for (typename RowSet::Iter r = rows.begin(); r != rows.end(); ++r) {
        T* s = start + r*n;
        for (typename ColSet::Iter c = cols.begin(); c != cols.end(); ++c) {
            if (*(s+c) != 0) { col_entries[c].push_back(s); break; }
        }
    }

    Size num_pivots = 0;
    for (typename ColSet::Iter c = cols.begin(); c != cols.end(); ++c) {
        std::vector<T*>& col = col_entries[c];
        while (col.size() > 1) {
            T* next = col[0];
            T nv = *(next+c);
            if (nv < 0) { nv *= -1; }
            if (nv != 1) {
                for (Index i = 1; i < (Index) col.size(); ++i) {
                    T iv = *(col[i]+c);
                    if (iv < 0) { iv *= -1; }
                    if (iv < nv) {
                        next = col[i];
                        nv = iv;
                        if (nv == 1) { break; }
                    }
                }
            }
            for (Index i = col.size()-1; i != -1; --i) {
                T* r = col[i];
                if (r != next) {
                    T m = *(r+c) / *(next+c);
                    for (Index j = c; j < n; ++j) { *(r+j) -= m**(next+j); }
                    if (*(r+c) == 0) {
                        col.erase(col.begin()+i);
                        typename ColSet::Iter k = c; ++k;
                        while (k != cols.end() && *(r+k) == 0) { ++k; }
                        if (k != cols.end()) { col_entries[k].push_back(r); }
                    }
                }
            }
        }
        if (!col.empty()) {
            ++num_pivots;
            if (pivots) { pivots->push_back(Pivots::value_type((col[0]-start)/n,c)); }
        }
    }
    return num_pivots;
}

template <class T> inline
Size
MatrixT<T>::row_diagonalise()
{
    return row_diagonalise(IndexSetR(0,m), IndexSetR(0,n));
}

template <class T> template <class RowSet, class ColSet>
Size
MatrixT<T>::row_diagonalise(const RowSet& rows, const ColSet& cols, Pivots* _pivots)
{
    Pivots pivots;
    Size s = row_triangulate(rows, cols, &pivots);
    T g0,p0,q0,p1,q1;
    for (Index i = 0; i < (Index) pivots.size(); ++i) {
        T* ivec = start + pivots[i].first*n;
        for (Index j = i+1; j < (Index) pivots.size(); ++j) {
            Index col = pivots[j].second;
            if (*(ivec+col) != 0) {
                T* jvec = start + pivots[j].first*n;
                euclidean(*(ivec+col),*(jvec+col),g0,p0,q0,p1,q1);
                if (p1 != 1) {
                    for (Index c = 0; c < col; ++c) {
                        *(ivec+c) = *(ivec+c)*p1;
                    }
                }
                for (Index c = col; c < n; ++c) {
                     *(ivec+c) = *(ivec+c)*p1 + *(jvec+c)*q1;
                }
            }
        }
    }
    if (_pivots != 0) { _pivots->swap(pivots); }
    return s;
}

template <class T> inline
void
MatrixT<T>::row_diagonalise(const Pivots& pivots)
{
    return row_diagonalise(pivots, IndexSetR(0,m));
}

template <class T> template <class RowSet>
void
MatrixT<T>::row_diagonalise(const Pivots& pivots, const RowSet& rows)
{
    T g0,p0,q0,p1,q1;
    for (typename RowSet::Iter r = rows.begin(); r != rows.end(); ++r) {
        T* rvec = start + r*n;
        for (Index i = 0; i < (Index) pivots.size(); ++i) {
            Index col = pivots[i].second;
            if (*(rvec+col) != 0) {
                T* pvec = start + pivots[i].first*n;
                if (rvec != pvec) {
                    euclidean(*(rvec+col),*(pvec+col),g0,p0,q0,p1,q1);
                    if (p1 != 1) {
                        for (Index c = 0; c < col; ++c) {
                            *(rvec+c) = *(rvec+c)*p1;
                        }                
                    }
                    for (Index c = col; c < n; ++c) {
                        *(rvec+c) = *(rvec+c)*p1 + *(pvec+c)*q1;
                    }                
                }
            }
        }
    }
}

template <class T>
void
MatrixT<T>::row_normalise()
{
    for (Index r = 0; r < m; ++r) {
        VectorDLW<T>(start+r*n, start+(r+1)*n).normalise();
    }
}

template <class T>
void
MatrixT<T>::print() const
{
    for (Index r = 0; r < m; ++r) {
        T* s = start + r*n;
        for (Index c = 0; c < n; ++c) { out->width(2); *out << *(s+c) << " "; }
        *out << "\n";
    }
    *out << std::endl;
}

} // namespace _4ti2_

#endif
