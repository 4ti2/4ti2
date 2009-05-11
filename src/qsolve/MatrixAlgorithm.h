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

namespace _4ti2_
{

template <class T, class IndexSet>
class MatrixAlgorithm : public QSolveAlgorithm<T,IndexSet>
{
public:
    MatrixAlgorithm();
    MatrixAlgorithm(QSolveConsOrder o);
    ~MatrixAlgorithm();

protected:
    virtual void compute(const ConeT<T>& cone, 
                VectorArrayT<T>& rays, std::vector<Index>& ray_ineqs, 
                VectorArrayT<T>& cirs, std::vector<Index>& cir_ineqs);

    void compute(const ConeT<T>& cone, VectorArrayT<T>& rays, std::vector<IndexSet>& supps, 
                Index& cons_added, std::vector<int>& ineqs);
    void compute(const ConeT<T>& cone, VectorArrayT<T>& rays, std::vector<IndexSet>& supps,
                Index& cons_added, VectorArrayT<T>& cirs, std::vector<int>& ineqs);

    void check(const ConeT<T>& cone, const IndexSet& rem, const VectorArrayT<T>& rays,
                const std::vector<IndexSet>& supps,
                const std::vector<IndexSet>& cir_supps);
};

template <class T, class IndexSet>
class MatrixRayAlgorithm : public ThreadedAlgorithm
{
public:
    MatrixRayAlgorithm(const ConeT<T>& cone, const VectorArrayT<T>& rays, const std::vector<IndexSet>& supps, 
                const IndexSet& rel, const Index& cons_added);
    virtual ~MatrixRayAlgorithm();

    virtual void compute();
    void compute(Index next, Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end);
    void threaded_compute(Index next, Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end);
    void transfer(VectorArrayT<T>& rays, std::vector<IndexSet>& supps);

protected:
    const ConeT<T>& cone;
    const VectorArrayT<T>& rays;
    const std::vector<IndexSet>& supps;
    const IndexSet& rel;
    const Index& cons_added;

    VectorArrayT<T> new_rays;
    std::vector<IndexSet> new_supps;

    Index _next;
    Index _r1_start;
    Index _r1_end;
    Index _r2_start;
    Index _r2_index;
    Index _r2_end;

    VectorT<T> temp;
    MatrixT<T> matrix;

    void project_cone(
                const ConeT<T>& cone,
                const IndexSet& zero_supp,
                MatrixT<T>& trans,
                std::vector<Index>& con_map);
    void zero_cols(
                const MatrixT<T>& matrix,
                const std::vector<Index>& con_map,
                IndexSet& zeros);
    bool is_two_dimensional_face(
                const MatrixT<T>& trans,
                const std::vector<Index>& con_map,
                const IndexSet& diff);

    void create_ray(Index r1, const T& s1, Index r2);
};

template <class T, class IndexSet>
class MatrixCirAlgorithm : public MatrixRayAlgorithm<T,IndexSet>
{
public:
    MatrixCirAlgorithm(const ConeT<T>& cone, const VectorArrayT<T>& rays, 
                const std::vector<IndexSet>& supps, const std::vector<IndexSet>& cir_supps, 
                const IndexSet& rel, const Index& cons_added);
    virtual ~MatrixCirAlgorithm();

    virtual void compute();
    void compute(Index next, Index r1_start, Index r1_end, Index r2_start, Index r2_end);
    void threaded_compute(Index next, Index r1_start, Index r1_end, Index r2_start, Index r2_end);
    void transfer(VectorArrayT<T>& rays, std::vector<IndexSet>& supps, std::vector<IndexSet>& cir_supps);

protected:
    const std::vector<IndexSet>& cir_supps;
    std::vector<IndexSet> new_cir_supps;

    void create_circuit(Index r1, const T& s1, Index r2);
};

} // namespace _4ti2_

#include "qsolve/MatrixAlgorithm.hpp"

#endif
