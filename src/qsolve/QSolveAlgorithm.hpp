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

#undef DEBUG_4ti2
#define DEBUG_4ti2(X) //X

using namespace _4ti2_;

template <class T, class IndexSet>
QSolveAlgorithm<T,IndexSet>::QSolveAlgorithm()
{
    // The default constraint ordering is MININDEX.
    // TODO: Is this what we want?
    compare = &minindex;
}

template <class T, class IndexSet>
QSolveAlgorithm<T,IndexSet>::QSolveAlgorithm(QSolveConsOrder o)
{
    set_constraint_order(o);
}

template <class T, class IndexSet>
QSolveAlgorithm<T,IndexSet>::~QSolveAlgorithm()
{
}

#if 1
template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::compute(
            const ConeT<T>& cone,
            VectorArrayT<T>& rays,
            VectorArrayT<T>& cirs,
            VectorArrayT<T>& subspace)
{
    DEBUG_4ti2(*out << "MATRIX:\n" << _orig_matrix << "\n";)
    DEBUG_4ti2(*out << "RELS:\n" << rels << "\n";)
    DEBUG_4ti2(*out << "SIGN:\n" << sign << "\n";)

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
        for (Index i = 0; i < n; ++i) {
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
        for (Index i = 0; i < n; ++i) {
            VectorT<T> cir(n,0);
            cir.set(i,1);
            cirs.insert(cir);
            cir_ineqs.push_back(i);
        }
        std::vector<Index> ray_ineqs;
        compute(cone, rays, ray_ineqs, cirs, cir_ineqs);
        return;
    }

    VectorArrayT<T> trans(n, num_cons);
    // Add an identity matrix at the beginning of the transpose.
    trans.assignT(0, IndexSetR(0,n), IndexSetR(0,n));
    for (Index i = 0; i < n; ++i) { trans[i][i] = 1; }
    // Add transpose after the identity matrix.
    trans.assign_trans(cone.get_matrix(), IndexSetR(0,m), IndexSetR(0,n), IndexSetR(0,n), IndexSetR(n,n+m));
    DEBUG_4ti2(*out << "TRANS:\n" << trans << "\n";)

    // Process the equality constraints.
    Index eq_row = upper_triangle(trans, 0, n, full_eq.begin(), full_eq.end());
    DEBUG_4ti2(*out << "EQ TRANS:\n" << trans << "\n";)
    trans.remove(0, eq_row);

    // Process the inequality constraints.
    IndexSetD rs_pivots(num_cons, false);
    Index rs_dim = diagonal(trans, 0, trans.get_number(), full_rs.begin(), full_rs.end(), rs_pivots);
    DEBUG_4ti2(*out << "RS PIVOTS:\n" << rs_pivots << "\n";)
    DEBUG_4ti2(*out << "RS TRANS:\n" << trans << "\n";)

    // Process the circuit constraints.
    IndexSetD cir_pivots(num_cons, false);
    Index cir_dim=diagonal(trans, rs_dim, trans.get_number(), full_cir.begin(), full_cir.end(), cir_pivots);
    DEBUG_4ti2(*out << "CIR PIVOTS:\n" << cir_pivots << "\n";)
    DEBUG_4ti2(*out << "CIR TRANS:\n" << trans << "\n";)

    // Extract a linear subspace basis.
    // TODO: Is the subspace basis necessarily independent?
    subspace.init(trans.get_number()-cir_dim, n);
    subspace.assign(trans, IndexSetR(cir_dim, trans.get_number()), IndexSetR(0,n));
    DEBUG_4ti2(*out << "SUB BASIS:\n" << subspace << "\n";)

    // Construct the projected constraint matrix.
    full_rs.set_difference(rs_pivots);
    full_cir.set_difference(cir_pivots);
    VectorArrayT<T> proj_matrix(full_rs.count()+full_cir.count(), cir_dim);
    // First add rs constraints.
    proj_matrix.assign_trans(trans, IndexSetR(0,cir_dim), full_rs, 
                    IndexSetR(0,full_rs.count()), IndexSetR(0,cir_dim));
    // Then add cir constraints.
    proj_matrix.assign_trans(trans, IndexSetR(0, cir_dim), full_cir, 
                    IndexSetR(full_rs.count(),full_rs.count()+full_cir.count()), IndexSetR(0,cir_dim));
    proj_matrix.normalise();
    DEBUG_4ti2(*out << "PROJ MATRIX:\n" << proj_matrix << "\n";)

    // Construct projected cone and projected initial rays.
    ConeT<T> proj_cone(proj_matrix);
    VectorArrayT<T> proj_rays(0, proj_matrix.get_size());
    std::vector<Index> proj_ray_ineqs;
    for (Index i = 0; i < rs_dim; ++i) {
        VectorT<T> ray(cir_dim,0);
        ray.set(i,1);
        proj_rays.insert(ray);
        proj_ray_ineqs.push_back(i);
    }
    VectorArrayT<T> proj_cirs(0, proj_matrix.get_size());
    std::vector<Index> proj_cir_ineqs;
    for (Index i = rs_dim; i < cir_dim; ++i) {
        proj_cone.set_constraint_type(i, _4ti2_DB);
        VectorT<T> cir(cir_dim,0);
        cir.set(i,1);
        proj_cirs.insert(cir);
        proj_cir_ineqs.push_back(i);
    }
    for (Index i = cir_dim+full_rs.count(); i < cir_dim+full_rs.count()+full_cir.count(); ++i) {
        proj_cone.set_constraint_type(i, _4ti2_DB);
    }
    DEBUG_4ti2(*out << "CONE CONS:\n" << proj_cone.get_constraint_types() << "\n";)

    compute(proj_cone, proj_rays, proj_ray_ineqs, proj_cirs, proj_cir_ineqs);

    // The map to lift the rays and circuits back into the original space.
    VectorArrayT<T> map(n, cir_dim);
    map.assign_trans(trans, IndexSetR(0,cir_dim), IndexSetR(0,n));
    DEBUG_4ti2(*out << "MAP:\n" << map << std::endl;)

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
#else
template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::compute(
            const ConeT<T>& cone,
            VectorArrayT<T>& rays,
            VectorArrayT<T>& cirs,
            VectorArrayT<T>& subspace)
{
    DEBUG_4ti2(*out << "MATRIX:\n" << matrix << "\n";)
    DEBUG_4ti2(*out << "RELS:\n" << rels << "\n";)
    DEBUG_4ti2(*out << "SIGN:\n" << sign << "\n";)

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
        for (Index i = 0; i < n; ++i) {
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
        for (Index i = 0; i < n; ++i) {
            VectorT<T> cir(n,0);
            cir.set(i,1);
            cirs.insert(cir);
            cir_ineqs.push_back(i);
        }
        std::vector<Index> ray_ineqs;
        compute(cone, rays, ray_ineqs, cirs, cir_ineqs);
        return;
    }

    VectorArrayT<T> trans(n, num_cons);
    // Add an identity matrix at the beginning of the transpose.
    for (Index i = 0; i < n; ++i) {
        for (Index j = 0; j < n; ++j) {
            if (i == j) { trans[i][j] = 1; }
            else { trans[i][j] = 0; }
        }
    } 
    // Add transpose after the identity matrix.
    trans.assign_trans(cone.get_matrix(), IndexSetR(0,m), IndexSetR(0,n), IndexSetR(0,n), IndexSetR(n,n+m));
    DEBUG_4ti2(*out << "TRANS:\n" << trans << "\n";)

    // Process the equality constraints.
    Index eq_row = upper_triangle(trans, 0, n, full_eq.begin(), full_eq.end());
    DEBUG_4ti2(*out << "EQ TRANS:\n" << trans << "\n";)
    trans.remove(0, eq_row);

    // Process the inequality constraints.
    IndexSetD rs_pivots(num_cons, false);
    Index rs_dim = diagonal(trans, 0, trans.get_number(), full_rs.begin(), full_rs.end(), rs_pivots);
    DEBUG_4ti2(*out << "RS PIVOTS:\n" << rs_pivots << "\n";)
    DEBUG_4ti2(*out << "RS TRANS:\n" << trans << "\n";)

    // Process the circuit constraints.
    IndexSetD cir_pivots(num_cons, false);
    Index cir_dim=diagonal(trans, rs_dim, trans.get_number(), full_cir.begin(), full_cir.end(), cir_pivots);
    DEBUG_4ti2(*out << "CIR PIVOTS:\n" << cir_pivots << "\n";)
    DEBUG_4ti2(*out << "CIR TRANS:\n" << trans << "\n";)

    // Extract a linear subspace basis.
    // TODO: Is the subspace basis necessarily independent?
    subspace.init(trans.get_number()-cir_dim, n);
    subspace.assign(trans, IndexSetR(cir_dim, trans.get_number()), IndexSetR(0,n));
    DEBUG_4ti2(*out << "SUB BASIS:\n" << subspace << "\n";)

    // Construct initial rays.
    rays.init(rs_dim, n);
    rays.assign(trans, IndexSetR(0,rs_dim), IndexSetR(0,n));
    DEBUG_4ti2(*out << "RAY BASIS:\n" << rays << "\n";)
    std::vector<Index> ray_ineqs;
    for (typename IndexSet::Iter it = rs_pivots.begin(); it != rs_pivots.end(); ++it) {
        ray_ineqs.push_back(*it);
    }
    DEBUG_4ti2(*out << "RAY INEQS:\n" << ray_ineqs << "\n";)

    // Construct initial circuits.
    cirs.init(cir_dim-rs_dim, n);
    cirs.assign(trans, IndexSetR(rs_dim,cir_dim), IndexSetR(0,n));
    DEBUG_4ti2(*out << "CIR BASIS:\n" << cirs << "\n";)
    std::vector<Index> cir_ineqs;
    for (typename IndexSet::Iter it = cir_pivots.begin(); it != cir_pivots.end(); ++it) {
        cir_ineqs.push_back(*it);
    }
    DEBUG_4ti2(*out << "CIR INEQS:\n" << cir_ineqs << "\n";)

    compute(cone, rays, ray_ineqs, cirs, cir_ineqs);

    return;
}
#endif


template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::set_constraint_order(QSolveConsOrder o)
{
    order = o;
    switch (o)
    {
    case MININDEX:
        compare = &minindex;
        break;
    case MAXCUTOFF:
        compare = &maxcutoff;
        break;
    case MINCUTOFF:
        compare = &mincutoff;
        break;
    case MAXINTER:
        compare = &maxinter;
        break;
    default:
        //compare = &maxinter;
        compare = &maxcutoff;
        //compare = &minindex;
        break;
    }
}


template <class T, class IndexSet>
bool
QSolveAlgorithm<T,IndexSet>::minindex(
                Index next_pos_count, Index next_neg_count, Index next_zero_count,
                Index pos_count, Index neg_count, Index zero_count)
{
    return false;
}

template <class T, class IndexSet>
bool
QSolveAlgorithm<T,IndexSet>::maxinter(
                Index next_pos_count, Index next_neg_count, Index next_zero_count,
                Index pos_count, Index neg_count, Index zero_count)
{
    return (zero_count > next_zero_count);
}

template <class T, class IndexSet>
bool
QSolveAlgorithm<T,IndexSet>::maxcutoff(
                Index next_pos_count, Index next_neg_count, Index next_zero_count,
                Index pos_count, Index neg_count, Index zero_count)
{
    return (neg_count > next_neg_count);
}

template <class T, class IndexSet>
bool
QSolveAlgorithm<T,IndexSet>::mincutoff(
                Index next_pos_count, Index next_neg_count, Index next_zero_count,
                Index pos_count, Index neg_count, Index zero_count)
{
    return (neg_count < next_neg_count);
}

// Count how many zero, positive and negative entries there are in a column.
template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::slack_count(
        const ConeT<T>& cone, const VectorArrayT<T>& rays, Index next_con,
        Index& pos_count, Index& neg_count, Index& zero_count)
{
    zero_count = 0; pos_count = 0; neg_count = 0;
    T slack;
    for (Index r = 0; r < rays.get_number(); ++r) {
        cone.get_slack(rays[r], next_con, slack);
        if (slack == 0) { ++zero_count; }
        else if (slack > 0) { ++pos_count; }
        else { ++neg_count; }
    }
}

template <class T, class IndexSet>
Index
QSolveAlgorithm<T,IndexSet>::next_circuit_constraint(
                const ConeT<T>& cone,
                const VectorArrayT<T>& vs,
                const IndexSet& rem)
{
    // Sanity Check
    assert(vs.get_size() == remaining.get_size());

    Index next_con = *rem.begin();
    if (order == MININDEX) { return next_con; }

    typename IndexSet::Iter it = rem.begin(); ++it;
    if (it == rem.end()) { return next_con; }

    Index next_pos_count, next_neg_count, next_zero_count;
    slack_count(cone, vs, next_con, next_pos_count, next_neg_count, next_zero_count);
    while (it != rem.end()) {
        Index pos_count, neg_count, zero_count;
        slack_count(cone, vs, *it, pos_count, neg_count, zero_count);
        if (maxinter(next_pos_count, next_neg_count, next_zero_count,
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


template <class T, class IndexSet>
Index
QSolveAlgorithm<T,IndexSet>::next_constraint(
                const ConeT<T>& cone,
                const VectorArrayT<T>& vs,
                const IndexSet& rem)
{
    // Sanity Check
    assert(vs.get_size() == remaining.get_size());

    Index next_con = *rem.begin();
    if (order == MININDEX) { return next_con; }

    typename IndexSet::Iter it = rem.begin(); ++it;
    if (it == rem.end()) { return next_con; }

    Size next_pos_count, next_neg_count, next_zero_count;
    slack_count(cone, vs, next_con, next_pos_count, next_neg_count, next_zero_count);
    while (it != rem.end()) {
        Size pos_count, neg_count, zero_count;
        slack_count(cone, vs, *it, pos_count, neg_count, zero_count);
        if ((*compare)(next_pos_count, next_neg_count, next_zero_count,
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

// Pushes zeros to the beginning.
template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::sort_nonzeros(
                const ConeT<T>& cone, 
                VectorArrayT<T>& rays,
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
template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::sort_positives(
                const ConeT<T>& cone, 
                VectorArrayT<T>& rays,
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
template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::sort_negatives(
                const ConeT<T>& cone, 
                VectorArrayT<T>& rays,
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
template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::sort_rays(
                VectorArrayT<T>& vs,
                std::vector<IndexSet>& supps,
                const IndexSet& ray_mask,
                Index start, Index end,
                Index& middle)
{
    Index index = start;
    for (Index i = start; i < end; ++i) {
        if (!ray_mask.set_disjoint(supps[i])) {
            vs.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes nonzeros to the beginning.
template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::sort_nonzeros(
                const ConeT<T>& cone, 
                VectorArrayT<T>& rays,
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
template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::sort_positives(
                const ConeT<T>& cone, 
                VectorArrayT<T>& rays,
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
template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::sort_negatives(
                const ConeT<T>& cone, 
                VectorArrayT<T>& rays,
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
template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::sort_rays(
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

template <class T, class IndexSet>
Index
QSolveAlgorithm<T,IndexSet>::next_constraint(
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

template <class T, class IndexSet>
Index
QSolveAlgorithm<T,IndexSet>::next_constraint(
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

template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::flip(std::vector<IndexSet>& supps, int start, int end)
{
    for (Index i = start; i < end; ++i) { supps[i].swap_odd_n_even(); }
}

// Splits rays into rays and circuits.
// NOTE: Does not update supports.
// TODO: cannot assume anything about the size of the supports.
template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::split_rays(
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

template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::update_supports(
                std::vector<IndexSet>& supps,
                Index index, Index start, Index end)
{
    for (Index i = start; i < end; ++i) { supps[i].set(index); }
}

template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::resize_supports(
                std::vector<IndexSet>& supps, Size size)
{
    for (Index i = 0; i < (Index) supps.size(); ++i) { supps[i].resize(size); }
}

template <class T, class IndexSet>
Index
QSolveAlgorithm<T,IndexSet>::next_constraint(
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

template <class T, class IndexSet>
void
QSolveAlgorithm<T,IndexSet>::print_debug_diagnostics(
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

#undef DEBUG_4ti2

