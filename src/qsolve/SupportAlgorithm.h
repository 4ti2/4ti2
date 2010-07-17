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

#ifndef _4ti2_qsolve__SupportAlgorithm_
#define _4ti2_qsolve__SupportAlgorithm_

#include "qsolve/QSolveAlgorithm.h"
#include "qsolve/Vector.h"
#include "qsolve/VectorArray.h"
#include "qsolve/QSolveConsOrder.h"
#include "qsolve/RayState.h"

#include "qsolve/BinaryTree.h"
#include "qsolve/FullTree.h"
#include "qsolve/MultiTree.h"

//#define SUPPORTTREE FullTree
//#define SUPPORTTREE BinaryTree
#define SUPPORTTREE MultiTree

namespace _4ti2_
{

template <class T>
class SupportAlgorithm : public QSolveAlgorithm<T>
{
public:
    SupportAlgorithm();
    SupportAlgorithm(QSolveConsOrder o);
    virtual ~SupportAlgorithm();

protected:
    virtual void compute(const ConeT<T>& cone, 
                VectorArrayT<T>& rays, std::vector<Index>& ray_ineqs, 
                VectorArrayT<T>& cirs, std::vector<Index>& cir_ineqs);

    template <class IndexSet>
    void compute(const ConeT<T>& cone, VectorArrayT<T>& rays, std::vector<IndexSet>& supps, 
                Index& cons_added, std::vector<int>& ineqs);
    template <class IndexSet>
    void compute(const ConeT<T>& cone, VectorArrayT<T>& rays, std::vector<IndexSet>& supps,
                Index& cons_added, VectorArrayT<T>& cirs, std::vector<int>& dbls);
};

template <class IndexSet>
class SupportSubAlgorithmBase
{
public:
    SupportSubAlgorithmBase(RayStateAPI<IndexSet>& helper, std::vector<IndexSet>& supps, const IndexSet& rel, const Index& cons_added,
                             const Index& next, const IndexSet& ray_mask, const SUPPORTTREE<IndexSet>& tree);
    virtual ~SupportSubAlgorithmBase();

    void compute_rays(Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end);
    void compute_cirs(Index r1_start, Index r1_end, Index r2_start, Index r2_end);
    void transfer();

protected:
    RayStateAPI<IndexSet>& helper;
    std::vector<IndexSet>& supps;
    //std::vector<IndexSet> new_supps;

    const IndexSet& rel;
    const Index& cons_added;
    const Index& next;
    const IndexSet& ray_mask;
    const SUPPORTTREE<IndexSet>& tree;

    //Index r1;
    //void set_r1_index(Index r1);
    //void create_ray(Index r2);
    //void create_circuit(Index r1, Index r2);
};

template <class IndexSet>
class SupportRayAlgorithm : public SupportSubAlgorithmBase<IndexSet>, public ThreadedAlgorithm
{
public:
    SupportRayAlgorithm(RayStateAPI<IndexSet>& helper, std::vector<IndexSet>& supps, 
                const IndexSet& rel, const Index& cons_added, const Index& next, const SUPPORTTREE<IndexSet>& tree,
                IndexRanges& indices);
    SupportRayAlgorithm* clone();

    virtual void compute();

protected:
    IndexRanges& indices;
};

template <class IndexSet>
class SupportCirAlgorithm : public SupportSubAlgorithmBase<IndexSet>, public ThreadedAlgorithm
{
public:
    SupportCirAlgorithm(RayStateAPI<IndexSet>& helper, std::vector<IndexSet>& supps, 
                const IndexSet& rel, const Index& cons_added, const Index& next, const SUPPORTTREE<IndexSet>& tree,
                const IndexSet& ray_mask, IndexRanges& indices);
    SupportCirAlgorithm* clone();

    virtual void compute();

protected:
    IndexRanges& indices;
};


} // namespace _4ti2_

// Definitions of template class functions.
#include "qsolve/SupportAlgorithm.hpp"

#endif
