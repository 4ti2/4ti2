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

#ifndef _4ti2_qsolve__MatrixAlgorithm_
#define _4ti2_qsolve__MatrixAlgorithm_

#include "qsolve/Vector.h"
#include "qsolve/VectorArray.h"
#include "qsolve/IndexSet.h"
#include "qsolve/QSolveConsOrder.h"
#include "qsolve/QSolveVariant.h"
#include "qsolve/Cone.h"
#include "qsolve/QSolveAlgorithm.h"
#include "qsolve/ThreadedAlgorithm.h"
#include "qsolve/Matrix.h"
#include "qsolve/ConeC.h"
#include "qsolve/RayState.h"

namespace _4ti2_
{

template <class T>
class MatrixAlgorithm : public QSolveAlgorithm<T>
{
public:
    MatrixAlgorithm();
    MatrixAlgorithm(QSolveConsOrder o);
    ~MatrixAlgorithm();

protected:
    virtual void compute(const ConeT<T>& cone, 
                VectorArrayT<T>& rays, std::vector<Index>& ray_ineqs, 
                VectorArrayT<T>& cirs, std::vector<Index>& cir_ineqs);

    template <class IndexSet>
    void compute(const ConeT<T>& cone, VectorArrayT<T>& rays, std::vector<IndexSet>& supps, 
                Index& cons_added, std::vector<int>& ineqs);
    template <class IndexSet>
    void compute(const ConeT<T>& cone, VectorArrayT<T>& rays, std::vector<IndexSet>& supps,
                Index& cons_added, VectorArrayT<T>& cirs, std::vector<int>& ineqs);

    template <class IndexSet>
    void check(const ConeT<T>& cone, const IndexSet& rem, const VectorArrayT<T>& rays,
                const std::vector<IndexSet>& supps,
                const std::vector<IndexSet>& cir_supps);
};


template <class IndexSet>
class MatrixSubAlgorithmBase
{
public:
    MatrixSubAlgorithmBase(RayStateAPI<IndexSet>& helper, std::vector<IndexSet>& supps, std::vector<IndexSet>& cir_supps, 
                        const IndexSet& rel, const Index& cons_added, const Index& next);
    virtual ~MatrixSubAlgorithmBase();

    void compute_rays(Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end);
    void compute_cirs(Index r1_start, Index r1_end, Index r2_start, Index r2_end);

    void transfer();

protected:
    RayStateAPI<IndexSet>& helper;
    std::vector<IndexSet>& supps;
    //std::vector<IndexSet> new_supps;
    std::vector<IndexSet>& cir_supps;
    //std::vector<IndexSet> new_cir_supps;

    const IndexSet& rel;
    const Index& cons_added;
    const Index& next;

    Index _r1;

#if 0
    void set_r1_index(Index r1);
    // The next function is virtual just so the compiler doesn't inline it!!!
    virtual void create_ray(Index r2);
    void create_rays(Index r1, std::vector<Index>& r2s);
    void create_circuit(Index r2);
#endif
};


template <class IndexSet>
class MatrixRayAlgorithm : public MatrixSubAlgorithmBase<IndexSet>, public ThreadedAlgorithm
{
public:
    MatrixRayAlgorithm(RayStateAPI<IndexSet>& helper,
                std::vector<IndexSet>& supps, const IndexSet& rel, const
                Index& cons_added, const Index& next, IndexRanges& indices);

    virtual void compute();
    MatrixRayAlgorithm* clone();

protected:
    IndexRanges& indices;
};

template <class IndexSet>
class MatrixCirAlgorithm : public MatrixSubAlgorithmBase<IndexSet>, public ThreadedAlgorithm
{
public:
    MatrixCirAlgorithm(RayStateAPI<IndexSet>& helper,
                std::vector<IndexSet>& supps, std::vector<IndexSet>& cir_supps, 
                const IndexSet& rel, const Index& cons_added, const Index& next, IndexRanges& indices);

    virtual void compute();
    MatrixCirAlgorithm* clone();

protected:
    IndexRanges& indices;
};

} // namespace _4ti2_

#include "qsolve/MatrixAlgorithm.hpp"

#endif
