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

//#include "qsolve/SupportTree.h"
//#include "qsolve/BinaryTree.h"
#include "qsolve/PlusTree.h"
#include "qsolve/FullTree.h"

#define SUPPORTTREE FullTree
//#define SUPPORTTREE PlusTree 

namespace _4ti2_
{

template <class T, class IndexSet>
class SupportAlgorithm : public QSolveAlgorithm<T,IndexSet>
{
public:
    SupportAlgorithm();
    SupportAlgorithm(QSolveConsOrder o);
    virtual ~SupportAlgorithm();

protected:
    virtual void compute(const ConeT<T>& cone, 
                VectorArrayT<T>& rays, std::vector<Index>& ray_ineqs, 
                VectorArrayT<T>& cirs, std::vector<Index>& cir_ineqs);

    void compute(const ConeT<T>& cone, VectorArrayT<T>& rays, std::vector<IndexSet>& supps, 
                Index& cons_added, std::vector<int>& ineqs);
    void compute1(const ConeT<T>& cone, VectorArrayT<T>& rays, std::vector<IndexSet>& supps, 
                Index& cons_added, std::vector<int>& ineqs);

    void compute(const ConeT<T>& cone, VectorArrayT<T>& rays, std::vector<IndexSet>& supps,
                Index& cons_added, VectorArrayT<T>& cirs, std::vector<int>& dbls);

    void compute(const ConeT<T>& cone, const SUPPORTTREE<IndexSet>& tree, const IndexSet& cir_mask,
                const char* buffer, VectorArrayT<T>& rays, std::vector<IndexSet>& supps,
                int next_col, int m1, int r1_start, int r1_end, int r2_start, int r2_end);
};

template <class T, class IndexSet>
class SupportRayAlgorithm : public ThreadedAlgorithm
{
public:
    SupportRayAlgorithm(const ConeT<T>& cone, const VectorArrayT<T>& rays, const std::vector<IndexSet>& supps, 
                const IndexSet& rel, const Index& cons_added, const Index& next, const SUPPORTTREE<IndexSet>& tree);
    virtual ~SupportRayAlgorithm();

    virtual void compute();
    void compute(Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end);
    void compute1(Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end);
    void threaded_compute(Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end);
    void transfer(VectorArrayT<T>& rays, std::vector<IndexSet>& supps);

protected:
    const ConeT<T>& cone;
    const VectorArrayT<T>& rays;
    const std::vector<IndexSet>& supps;
    const IndexSet& rel;
    const Index& cons_added;
    const Index& next;
    const SUPPORTTREE<IndexSet>& tree;

    VectorArrayT<T> new_rays;
    std::vector<IndexSet> new_supps;

    Index _r1_start;
    Index _r1_end;
    Index _r2_start;
    Index _r2_index;
    Index _r2_end;

    VectorT<T> temp;
    void create_ray(Index r1, const T& s1, Index r2);
};

template <class T, class IndexSet>
class SupportCirAlgorithm : public SupportRayAlgorithm<T,IndexSet>
{
public:
    SupportCirAlgorithm(const ConeT<T>& cone, const VectorArrayT<T>& rays, const std::vector<IndexSet>& supps, 
                const IndexSet& rel, const Index& cons_added, const Index& next, const SUPPORTTREE<IndexSet>& tree, const IndexSet& ray_mask);
    virtual ~SupportCirAlgorithm();

    virtual void compute();
    void compute(Index r1_start, Index r1_end, Index r2_start, Index r2_end);
    void threaded_compute(Index r1_start, Index r1_end, Index r2_start, Index r2_end);

protected:
    const IndexSet& ray_mask;

    void create_circuit(Index r1, const T& s1, Index r2);
};


} // namespace _4ti2_

// Definitions of template class functions.
#include "qsolve/SupportAlgorithm.hpp"

#endif
