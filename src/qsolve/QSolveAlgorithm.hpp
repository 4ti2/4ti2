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

#include "qsolve/QSolveAlgorithm.h"
#include "qsolve/Globals.h"
#include "qsolve/HermiteAlgorithm.h"
#include "qsolve/VectorStream.h"
#include "qsolve/Debug.h"
#include "qsolve/Stream.h"
#include "qsolve/IndexSetD.h"
#include "qsolve/IndexSetStream.h"
#include "4ti2/4ti2.h"
#include "qsolve/MatrixAlgorithm.h"
#include "qsolve/SupportAlgorithm.h"
#include "qsolve/RayState.h"

#undef DEBUG_4ti2
#define DEBUG_4ti2(X) //X

using namespace _4ti2_;

template <class T>
QSolveAlgorithm<T>::QSolveAlgorithm()
{
}

template <class T>
QSolveAlgorithm<T>::QSolveAlgorithm(QSolveVariant v, QSolveConsOrder o)
    : variant(v)
{
    order.set_constraint_order(o);
}

template <class T>
QSolveAlgorithm<T>::~QSolveAlgorithm()
{
}

template <class T>
void
QSolveAlgorithm<T>::compute(
                            const ConeT<T>& cone,
                            VectorArrayT<T>& rays,
                            VectorArrayT<T>& cirs,
                            VectorArrayT<T>& subspace)
{
    DEBUG_4ti2(*out << "MATRIX:\n" << cone.get_matrix() << "\n";)
    //DEBUG_4ti2(*out << "RELS:\n" << rels << "\n";)
    //DEBUG_4ti2(*out << "SIGN:\n" << sign << "\n";)

    Size n = cone.num_vars();
    Size m = cone.num_cons();
    Size num_cons = n+m;

    IndexSetD full_rs(num_cons,0);
    cone.constraint_set(_4ti2_LB, full_rs);
    IndexSetD full_cir(num_cons,0);
    cone.constraint_set(_4ti2_DB, full_cir);
    IndexSetD full_eq(num_cons,0);
    cone.constraint_set(_4ti2_EQ, full_eq);

    DEBUG_4ti2(*out << "RS:\n" << full_rs << "\n";)
    DEBUG_4ti2(*out << "CIR:\n" << full_cir << "\n";)
    DEBUG_4ti2(*out << "EQ:\n" << full_eq << "\n";)

    // If there are only ray components...
    // TODO: This could be better.
    if (full_rs.full()) {
        // Construct initial rays.
        std::vector<Index> ray_ineqs;
        for (Index i = 0; i < cone.num_vars(); ++i) {
            VectorT<T> ray(n,0);
            ray.set(i,1);
            rays.insert(ray);
            ray_ineqs.push_back(i);
        }
        std::vector<Index> cir_ineqs;
        compute(cone, rays, ray_ineqs, cirs, cir_ineqs);
        return;
    }

    // If there are only circuit components...
    // TODO: This could be better.
    if (full_cir.full()) {
        // Construct initial circuits.
        std::vector<Index> cir_ineqs;
        for (Index i = 0; i < cone.num_vars(); ++i) {
            VectorT<T> cir(n,0);
            cir.set(i,1);
            cirs.insert(cir);
            cir_ineqs.push_back(i);
        }
        std::vector<Index> ray_ineqs;
        compute(cone, rays, ray_ineqs, cirs, cir_ineqs);
        return;
    }

    ConeT<T> proj_cone;
    VectorArrayT<T> map;
    // Construct projected cone.
    cone.canonize(proj_cone, subspace, map);

    // Construct initial rays.
    VectorArrayT<T> proj_rays(0, proj_cone.num_vars());
    VectorArrayT<T> proj_cirs(0, proj_cone.num_vars());
    std::vector<Index> proj_ray_ineqs;
    std::vector<Index> proj_cir_ineqs;
    for (Index i = 0; i < proj_cone.num_vars(); ++i) {
        if (proj_cone.get_constraint_type(i) == _4ti2_LB) {
            VectorT<T> ray(proj_cone.num_vars(),0);
            ray.set(i,1);
            proj_rays.insert(ray);
            proj_ray_ineqs.push_back(i);
        }
        else {
            VectorT<T> cir(proj_cone.num_vars(),0);
            cir.set(i,1);
            proj_cirs.insert(cir);
            proj_cir_ineqs.push_back(i);
        }
    }
    DEBUG_4ti2(*out << "CONE CONS:\n" << proj_cone.get_constraint_types() << "\n";)

    compute(proj_cone, proj_rays, proj_ray_ineqs, proj_cirs, proj_cir_ineqs);

    // Lift rays back to original space.
    // TODO: This could be done faster.
    rays.init(proj_rays.get_number(), n);
    VectorArrayT<T>::dot(map, proj_rays, rays);
    rays.normalise();

    // Lift circuits back to original space.
    // TODO: This could be done faster.
    cirs.init(proj_cirs.get_number(), n);
    VectorArrayT<T>::dot(map, proj_cirs, cirs);
    cirs.normalise();

    return;
}

template <class T>
void
QSolveAlgorithm<T>::compute(const ConeT<T>& cone, 
            VectorArrayT<T>& rays, std::vector<Index>& ray_ineqs, 
            VectorArrayT<T>& cirs, std::vector<Index>& cir_ineqs)
{
    // The number of constraints added so far.
    Index cons_added = 0;
    Index next = -1;
    IndexSetD lbs(cone.num_vars()+cone.num_cons());
    cone.constraint_set(_4ti2_LB, lbs);
    IndexSetD dbs(cone.num_vars()+cone.num_cons());
    cone.constraint_set(_4ti2_DB, dbs);

    if (lbs.count()+2*dbs.count() <= IndexSetDS::max_size) {
        std::vector<IndexSetDS> supps;
        RayState<T,IndexSetDS> state(cone, rays, supps, next);
        if (variant == MATRIX) {
            MatrixAlgorithm<IndexSetDS> alg(order);
            // Compute ray only constraints first.
            alg.compute_rays(cone, state, supps, next, cons_added, ray_ineqs);
            // Compute circuits.
            rays.transfer(cirs, 0, cirs.get_number(), rays.get_number());
            alg.compute_cirs(cone, state, supps, next, cons_added, cir_ineqs);
        } else {
            SupportAlgorithm<IndexSetDS> alg(order);
            // Compute ray only constraints first.
            alg.compute_rays(cone, state, supps, next, cons_added, ray_ineqs);
            // Compute circuits.
            rays.transfer(cirs, 0, cirs.get_number(), rays.get_number());
            alg.compute_cirs(cone, state, supps, next, cons_added, cir_ineqs);
        }
        // Separate the rays from the circuits.
        split_rays(cone, supps, rays, cirs);
    } else {
        std::vector<IndexSetD> supps;
        RayState<T,IndexSetD> state(cone, rays, supps, next);
        if (variant == MATRIX) {
            MatrixAlgorithm<IndexSetD> alg(order);
            // Compute ray only constraints first.
            alg.compute_rays(cone, state, supps, next, cons_added, ray_ineqs);
            // Compute circuits.
            rays.transfer(cirs, 0, cirs.get_number(), rays.get_number());
            alg.compute_cirs(cone, state, supps, next, cons_added, cir_ineqs);
        } else {
            SupportAlgorithm<IndexSetD> alg(order);
            // Compute ray only constraints first.
            alg.compute_rays(cone, state, supps, next, cons_added, ray_ineqs);
            // Compute circuits.
            rays.transfer(cirs, 0, cirs.get_number(), rays.get_number());
            alg.compute_cirs(cone, state, supps, next, cons_added, cir_ineqs);
        }
        // Separate the rays from the circuits.
        split_rays(cone, supps, rays, cirs);
    }
}

#if 0
    template <class IndexSet>
    Index next_constraint(
            const ConeT<T>& cone, const IndexSet& rem, 
            VectorArrayT<T>& rays, std::vector<IndexSet>& supps,
            Index& pos_start, Index& pos_end, Index& neg_start, Index& neg_end);

    template <class IndexSet>
    void sort_nonzeros(
                const ConeT<T>& cone, VectorArrayT<T>& vs, std::vector<IndexSet>& supps,
                Index start, Index end, Index next_col, Index& middle);
    template <class IndexSet>
    void sort_positives(
                const ConeT<T>& cone, VectorArrayT<T>& vs, std::vector<IndexSet>& supps,
                Index start, Index end, Index next_col, Index& middle);
    template <class IndexSet>
    void sort_negatives(
                const ConeT<T>& cone, VectorArrayT<T>& vs, std::vector<IndexSet>& supps,
                Index start, Index end, Index next_col, Index& middle);
    template <class IndexSet>
    void sort_rays(
                VectorArrayT<T>& vs, std::vector<IndexSet>& supps,
                const IndexSet& ray_mask, Index start, Index end, Index& middle);

    template <class IndexSet>
    Index next_constraint(
                const ConeT<T>& cone, VectorArrayT<T>& vs, std::vector<IndexSet>& supps,
                const IndexSet& rem, const IndexSet& ray_mask,
                Index& ray_start, Index& ray_end, Index& cir_start, Index& cir_end,
                Index& pos_ray_start, Index& pos_ray_end, Index& neg_ray_start, Index& neg_ray_end,
                Index& pos_cir_start, Index& pos_cir_end, Index& neg_cir_start, Index& neg_cir_end);

    template <class IndexSet>
    Index next_constraint(const ConeT<T>& cone, const VectorArrayT<T>& vs, const IndexSet& rem);

    template <class IndexSet>
    Index next_circuit_constraint(const ConeT<T>& cone, const VectorArrayT<T>& vs, const IndexSet& rem);

    template <class IndexSet>
    void sort_nonzeros(
                const ConeT<T>& cone, VectorArrayT<T>& vs, std::vector<IndexSet>& supps,
                std::vector<IndexSet>& cir_supps,
                Index start, Index end, Index next_col, Index& middle);
    template <class IndexSet>
    void sort_positives(
                const ConeT<T>& cone, VectorArrayT<T>& vs, std::vector<IndexSet>& supps,
                std::vector<IndexSet>& cir_supps,
                Index start, Index end, Index next_col, Index& middle);
    template <class IndexSet>
    void sort_negatives(
                const ConeT<T>& cone, VectorArrayT<T>& vs, std::vector<IndexSet>& supps,
                std::vector<IndexSet>& cir_supps,
                Index start, Index end, Index next_col, Index& middle);
    template <class IndexSet>
    void sort_rays(
                VectorArrayT<T>& vs, std::vector<IndexSet>& supps,
                std::vector<IndexSet>& cir_supps,
                const IndexSet& ray_mask, Index start, Index end, Index& middle);

    template <class IndexSet>
    Index next_constraint(
                const ConeT<T>& cone, VectorArrayT<T>& vs, std::vector<IndexSet>& supps,
                std::vector<IndexSet>& cir_supps,
                const IndexSet& rem, const IndexSet& ray_mask,
                Index& ray_start, Index& ray_end, Index& cir_start, Index& cir_end,
                Index& pos_ray_start, Index& pos_ray_end, Index& neg_ray_start, Index& neg_ray_end,
                Index& pos_cir_start, Index& pos_cir_end, Index& neg_cir_start, Index& neg_cir_end);

    template <class IndexSet>
    void flip(std::vector<IndexSet>& supps, Index start, Index end);

    template <class IndexSet>
    void update_supports(
                    std::vector<IndexSet>& supps,
                    Index index, Index start, Index end);
    template <class IndexSet>
    void resize_supports(std::vector<IndexSet>& supps, Size size);

    template <class IndexSet>
    void print_debug_diagnostics(const ConeT<T>& cone, 
                const VectorArrayT<T>& rays, const std::vector<IndexSet>& supps, Index next);
#endif


#if 0
// Pushes zeros to the beginning.
template <class T> template <class IndexSet>
void
QSolveAlgorithm<T>::sort_nonzeros(
                const ConeT<T>& cone, VectorArrayT<T>& rays,
                std::vector<IndexSet>& supps,
                Index start, Index end,
                Index next_col, Index& middle)
{
    assert(start >= 0 && start <= end && end <= rays.get_number());
    T slack;
    Index index = start;
    for (Index i = start; i < end; ++i) {
        cone.get_slack(rays[i], next_col, slack);
        if (slack != 0) {
            rays.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes positives to the beginning.
template <class T> template <class IndexSet>
void
QSolveAlgorithm<T>::sort_positives(
                const ConeT<T>& cone, VectorArrayT<T>& rays,
                std::vector<IndexSet>& supps,
                Index start, Index end,
                Index next_col, Index& middle)
{
    assert(start >= 0 && start <= end && end <= rays.get_number());
    T slack;
    Index index = start;
    for (Index i = start; i < end; ++i) {
        cone.get_slack(rays[i], next_col, slack);
        if (slack > 0) {
            rays.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes positives to the beginning.
template <class T> template <class IndexSet>
void
QSolveAlgorithm<T>::sort_negatives(
                const ConeT<T>& cone, VectorArrayT<T>& rays,
                std::vector<IndexSet>& supps,
                Index start, Index end,
                Index next_col, Index& middle)
{
    assert(start >= 0 && start <= end && end <= rays.get_number());
    T slack;
    Index index = start;
    for (Index i = start; i < end; ++i) {
        cone.get_slack(rays[i], next_col, slack);
        if (slack < 0) {
            rays.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes rays to the beginning.
template <class T> template <class IndexSet>
void
QSolveAlgorithm<T>::sort_rays(
                VectorArrayT<T>& rays,
                std::vector<IndexSet>& supps,
                const IndexSet& ray_mask,
                Index start, Index end,
                Index& middle)
{
    Index index = start;
    for (Index i = start; i < end; ++i) {
        if (!ray_mask.set_disjoint(supps[i])) {
            rays.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes nonzeros to the beginning.
template <class T> template <class IndexSet>
void
QSolveAlgorithm<T>::sort_nonzeros(
                const ConeT<T>& cone, VectorArrayT<T>& rays,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& cir_supps,
                Index start, Index end,
                Index next_col, Index& middle)
{
    assert(start >= 0 && start <= end && end <= rays.get_number());
    T slack;
    Index index = start;
    for (Index i = start; i < end; ++i) {
        cone.get_slack(rays[i], next_col, slack);
        if (slack != 0) {
            rays.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            IndexSet::swap(cir_supps[i], cir_supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes positives to the beginning.
template <class T> template <class IndexSet>
void
QSolveAlgorithm<T>::sort_positives(
                const ConeT<T>& cone, VectorArrayT<T>& rays,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& cir_supps,
                Index start, Index end,
                Index next_col, Index& middle)
{
    assert(start >= 0 && start <= end && end <= rays.get_number());
    T slack;
    Index index = start;
    for (Index i = start; i < end; ++i) {
        cone.get_slack(rays[i], next_col, slack);
        if (slack > 0) {
            rays.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            IndexSet::swap(cir_supps[i], cir_supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes positives to the beginning.
template <class T> template <class IndexSet>
void
QSolveAlgorithm<T>::sort_negatives(
                const ConeT<T>& cone, VectorArrayT<T>& rays,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& cir_supps,
                int start, int end,
                int next_col, int& middle)
{
    assert(start >= 0 && start <= end && end <= rays.get_number());
    T slack;
    int index = start;
    for (int i = start; i < end; ++i) {
        cone.get_slack(rays[i], next_col, slack);
        if (slack < 0) {
            rays.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            IndexSet::swap(cir_supps[i], cir_supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes rays to the beginning.
template <class T> template <class IndexSet>
void
QSolveAlgorithm<T>::sort_rays(
                VectorArrayT<T>& vs,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& cir_supps,
                const IndexSet& ray_mask,
                int start, int end,
                int& middle)
{
    int index = start;
    for (int i = start; i < end; ++i) {
        if (!ray_mask.set_disjoint(supps[i])) {
            vs.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            IndexSet::swap(cir_supps[i], cir_supps[index]);
            ++index;
        }
    }
    middle = index;
}

template <class T> template <class IndexSet>
Index
QSolveAlgorithm<T>::next_constraint(
                const ConeT<T>& cone, VectorArrayT<T>& vs, std::vector<IndexSet>& supps,
                const IndexSet& rem, const IndexSet& ray_mask,
                Index& ray_start, Index& ray_end, Index& cir_start, Index& cir_end,
                Index& pos_ray_start, Index& pos_ray_end, Index& neg_ray_start, Index& neg_ray_end,
                Index& pos_cir_start, Index& pos_cir_end, Index& neg_cir_start, Index& neg_cir_end)
{
    assert(!rem.empty());
    // First, we choose the next constraint to add.
    Index next_col = next_circuit_constraint(cone, vs, rem);
    
    Index start = 0; Index end = vs.get_number(); Index middle;
    // We sort the vectors into nonzeros and then zeros.
    sort_nonzeros(cone, vs, supps, start, end, next_col, middle);
    Index nonzero_start = start, nonzero_end = middle;
    // We sort the nonzeros into rays and circuits.
    sort_rays(vs, supps, ray_mask, nonzero_start, nonzero_end, middle);
    ray_start = nonzero_start; ray_end = middle;
    cir_start = middle; cir_end = nonzero_end;
    // We sort the rays into positives and then negatives.
    sort_positives(cone, vs, supps, ray_start, ray_end, next_col, middle);
    pos_ray_start = ray_start; pos_ray_end = middle;
    neg_ray_start = middle; neg_ray_end = ray_end;
    // We sort the circuits into positives and then negatives.
    sort_positives(cone, vs, supps, cir_start, cir_end, next_col, middle);
    pos_cir_start = cir_start; pos_cir_end = middle;
    neg_cir_start = middle; neg_cir_end = cir_end;

    return next_col;
}

template <class T> template <class IndexSet>
Index
QSolveAlgorithm<T>::next_constraint(
                const ConeT<T>& cone, VectorArrayT<T>& vs,
                std::vector<IndexSet>& supps, std::vector<IndexSet>& cir_supps,
                const IndexSet& rem, const IndexSet& ray_mask,
                Index& ray_start, Index& ray_end, Index& cir_start, Index& cir_end,
                Index& pos_ray_start, Index& pos_ray_end, Index& neg_ray_start, Index& neg_ray_end,
                Index& pos_cir_start, Index& pos_cir_end, Index& neg_cir_start, Index& neg_cir_end)
{
    assert(!rem.empty());
    // First, we choose the next constraint to add.
    Index next_col = next_circuit_constraint(cone, vs, rem);

    Index start = 0; Index end = vs.get_number(); Index middle;
    // We sort the vectors into nonzeros and then zeros.
    sort_nonzeros(cone, vs, supps, cir_supps, start, end, next_col, middle);
    Index nonzero_start = start, nonzero_end = middle;
    // We sort the nonzeros into rays and circuits.
    sort_rays(vs, supps, cir_supps, ray_mask, nonzero_start, nonzero_end, middle);
    ray_start = nonzero_start; ray_end = middle;
    cir_start = middle; cir_end = nonzero_end;
    // We sort the rays into positives and then negatives.
    sort_positives(cone, vs, supps, cir_supps, ray_start, ray_end, next_col, middle);
    pos_ray_start = ray_start; pos_ray_end = middle;
    neg_ray_start = middle; neg_ray_end = ray_end;
    // We sort the circuits into positives and the negatives.
    sort_positives(cone, vs, supps, cir_supps, cir_start, cir_end, next_col, middle);
    pos_cir_start = cir_start; pos_cir_end = middle;
    neg_cir_start = middle; neg_cir_end = cir_end;

    return next_col;
}

template <class T> template <class IndexSet>
Index
QSolveAlgorithm<T>::next_constraint(
                const ConeT<T>& cone,
                const VectorArrayT<T>& vs,
                const IndexSet& rem)
{
    Index next_con = *rem.begin();
    if (order.get_constraint_order() == MININDEX) { return next_con; }

    typename IndexSet::Iter it = rem.begin(); ++it;
    if (it == rem.end()) { return next_con; }

    Size next_pos_count, next_neg_count, next_zero_count;
    cone.slack_count(vs, next_con, next_pos_count, next_neg_count, next_zero_count);
    while (it != rem.end()) {
        Size pos_count, neg_count, zero_count;
        cone.slack_count(vs, *it, pos_count, neg_count, zero_count);
        if ((*order.compare)(next_pos_count, next_neg_count, next_zero_count,
                        pos_count, neg_count, zero_count)) {
            next_con = *it;
            next_pos_count = pos_count;
            next_neg_count = neg_count;
            next_zero_count = zero_count;
        }
        ++it;
    }
    DEBUG_4ti2(*out << "Next Constraint is " << next_con << "\n";)
    return next_con;
}

template <class T> template <class IndexSet>
Index
QSolveAlgorithm<T>::next_circuit_constraint(
                const ConeT<T>& cone,
                const VectorArrayT<T>& vs,
                const IndexSet& rem)
{
    Index next_con = *rem.begin();
    if (order.get_constraint_order() == MININDEX) { return next_con; }

    typename IndexSet::Iter it = rem.begin(); ++it;
    if (it == rem.end()) { return next_con; }

    Size next_pos_count, next_neg_count, next_zero_count;
    cone.slack_count(vs, next_con, next_pos_count, next_neg_count, next_zero_count);
    while (it != rem.end()) {
        Size pos_count, neg_count, zero_count;
        cone.slack_count(vs, *it, pos_count, neg_count, zero_count);
        if ((*order.circuit_compare)(next_pos_count, next_neg_count, next_zero_count,
                        pos_count, neg_count, zero_count)) {
            next_con = *it;
            next_pos_count = pos_count;
            next_neg_count = neg_count;
            next_zero_count = zero_count;
        }
        ++it;
    }
    DEBUG_4ti2(*out << "Next Constraint is " << next_con << "\n";)
    return next_con;
}

template <class T> template <class IndexSet>
void
QSolveAlgorithm<T>::flip(std::vector<IndexSet>& supps, int start, int end)
{
    for (Index i = start; i < end; ++i) { supps[i].swap_odd_n_even(); }
}

template <class T> template <class IndexSet>
void
QSolveAlgorithm<T>::update_supports(
                std::vector<IndexSet>& supps,
                Index index, Index start, Index end)
{
    for (Index i = start; i < end; ++i) { supps[i].set(index); }
}

template <class T> template <class IndexSet>
void
QSolveAlgorithm<T>::resize_supports(
                std::vector<IndexSet>& supps, Size size)
{
    for (Index i = 0; i < (Index) supps.size(); ++i) { supps[i].resize(size); }
}

template <class T> template <class IndexSet>
Index
QSolveAlgorithm<T>::next_constraint(
            const ConeT<T>& cone,
            const IndexSet& rem,
            VectorArrayT<T>& rays,
            std::vector<IndexSet>& supps,
            Index& pos_start, Index& pos_end,
            Index& neg_start, Index& neg_end)
{
    assert(!rem.empty());
    // First, we choose the next constraint to add.
    Index next_con = next_constraint(cone, rays, rem);

    // TODO: Should we use a vector of slacks?
    Index start = 0; Index end = rays.get_number(); Index middle;
    // We sort the vectors into nonzeros and then zeros.
    sort_nonzeros(cone, rays, supps, start, end, next_con, middle);
    Index nonzero_start = start, nonzero_end = middle;
    // We sort the rays into positives and then negatives.
    sort_positives(cone, rays, supps, nonzero_start, nonzero_end, next_con, middle);
    pos_start = nonzero_start; pos_end = middle;
    neg_start = middle; neg_end = nonzero_end;

    return next_con;
}

template <class T> template <class IndexSet>
void
QSolveAlgorithm<T>::print_debug_diagnostics(
            const ConeT<T>& cone, 
            const VectorArrayT<T>& rays,
            const std::vector<IndexSet>& supps,
            Index next)
{
    *out << "Diagnostics.\n";
    T s;
    for (Index i = 0; i < rays.get_number(); ++i) {
        *out << "Ray" << i << " " << rays[i] << "\t";
        *out << "Sup " << supps[i] << "\t";
        cone.get_slack(rays[i], next, s);
        *out << "Sla " << s << "\n";
    }
}
#endif

// Splits rays into rays and circuits.
// NOTE: Does not update supports.
// TODO: cannot assume anything about the size of the supports.
template <class T> template <class IndexSet>
void
QSolveAlgorithm<T>::split_rays(
                const ConeT<T>& cone,
                const std::vector<IndexSet>& supps,
                VectorArrayT<T>& rays,
                VectorArrayT<T>& cirs)
{
    IndexSet ray_mask(cone.num_vars()+cone.num_cons());
    cone.constraint_set(_4ti2_LB, ray_mask);
    int index = 0;
    for (int i = 0; i < rays.get_number(); ++i) {
        if (!ray_mask.set_disjoint(supps[i])) {
            rays.swap_vectors(i,index);
            ++index;
        }
    }
    cirs.transfer(rays, index, rays.get_number(), 0);
}

inline
IndexRanges::IndexRanges()
{
    pthread_mutex_init(&mutex1, 0);
    r1_start = r1_index = r1_end = r2_start = r2_index = r2_end = increment = 0;
}

inline void
IndexRanges::init(Index _r1_start, Index _r1_end, Index _r2_start, Index _r2_index, Index _r2_end)
{
    r1_start = r1_index = _r1_start;
    r1_end = _r1_end;
    r2_start = _r2_start;
    r2_index = _r2_index;
    r2_end = _r2_end;
    if (Globals::num_threads==1) { increment = r1_end-r1_start; }
    else { increment = 1000; } // TODO
    //increment = (r1_end-r1_start)/Globals::num_threads+1;
}

inline void
IndexRanges::next(Index& _r1_start, Index& _r1_end, Index& _r2_start, Index& _r2_index, Index& _r2_end)
{
    pthread_mutex_lock(&mutex1);
    _r1_start = r1_index;
    r1_index = std::min(r1_index+increment, r1_end);
    _r1_end = r1_index;
    _r2_start = r2_start;
    _r2_index = r2_index;
    _r2_end = r2_end;
    pthread_mutex_unlock(&mutex1);
}

#undef DEBUG_4ti2

