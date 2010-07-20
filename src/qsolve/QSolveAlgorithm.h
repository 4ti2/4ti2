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

#ifndef _4ti2_qsolve__QSolveAlgorithm_
#define _4ti2_qsolve__QSolveAlgorithm_

#include "qsolve/Vector.h"
#include "qsolve/VectorArray.h"
#include "qsolve/Cone.h"
#include "qsolve/QSolveConsOrder.h"
#include "qsolve/QSolveVariant.h"

namespace _4ti2_
{

template <class T>
class QSolveAlgorithm
{
public:
    QSolveAlgorithm();
    QSolveAlgorithm(QSolveVariant v, QSolveConsOrder o);
    virtual ~QSolveAlgorithm();

    void compute(   const ConeT<T>& cone,
                    VectorArrayT<T>& rays,
                    VectorArrayT<T>& circuits,
                    VectorArrayT<T>& subspace);

protected:
    void compute(const ConeT<T>& cone, 
                VectorArrayT<T>& rays, std::vector<Index>& ray_ineqs, 
                VectorArrayT<T>& cirs, std::vector<Index>& cir_ineqs);

    template <class IndexSet>
    void split_rays(
                const ConeT<T>& cone,
                const std::vector<IndexSet>& supps,
                VectorArrayT<T>& rays,
                VectorArrayT<T>& cirs);

    ConsOrder order;
    QSolveVariant variant;
};

class IndexRanges {
public:
    IndexRanges();
    void init(Index _r1_start, Index _r1_end, Index _r2_start, Index _r2_index, Index _r2_end);
    void next(Index& _r1_start, Index& _r2_end, Index& _r2_start, Index& _r2_index, Index& _r2_end);
private:
    Index r1_start;
    Index r1_index;
    Index r1_end;
    Index r2_start;
    Index r2_index;
    Index r2_end;
    Size increment;
    pthread_mutex_t mutex1;
};

class IndexCirRanges {
public:
    IndexCirRanges();
    void init(Index pos_ray_start, Index pos_ray_start, Index neg_ray_start, Index neg_ray_end, Index cir_start, Index cir_end);
    void next(Index& _r1_start, Index& _r2_end, Index& _r2_start, Index& _r2_end);
private:
    Index r1_start;
    Index r1_index;
    Index r1_end;
    Index r2_start;
    Index r2_index;
    Index r2_end;
    Size increment;
    pthread_mutex_t mutex1;
};

} // namespace _4ti2_

// Definitions of template class functions.
#include "qsolve/QSolveAlgorithm.hpp"

#endif
