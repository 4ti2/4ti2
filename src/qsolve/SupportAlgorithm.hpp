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

#include <pthread.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "qsolve/SupportAlgorithm.h"
#include "qsolve/Globals.h"
#include "qsolve/Timer.h"
#include "qsolve/Debug.h"
#include "qsolve/Cone.h"

#include "qsolve/VectorStream.h"
#include "qsolve/VectorArrayStream.h"
#include "qsolve/IndexSetStream.h"
#include "qsolve/Stream.h"

#undef DEBUG_4ti2
#define DEBUG_4ti2(X) //X

using namespace _4ti2_;

template <class T>
SupportAlgorithm<T>::SupportAlgorithm()
    : QSolveAlgorithm<T>()
{
}

template <class T>
SupportAlgorithm<T>::SupportAlgorithm(QSolveConsOrder o)
    : QSolveAlgorithm<T>(o)
{
}

template <class T>
SupportAlgorithm<T>::~SupportAlgorithm()
{
}

template <class T>
void
SupportAlgorithm<T>::compute(const ConeT<T>& cone, 
            VectorArrayT<T>& rays, std::vector<Index>& ray_ineqs, 
            VectorArrayT<T>& cirs, std::vector<Index>& cir_ineqs)
{
    // The number of constraints added so far.
    Index cons_added = 0;
    IndexSetD lbs(cone.num_vars()+cone.num_cons());
    cone.constraint_set(_4ti2_LB, lbs);
    IndexSetD dbs(cone.num_vars()+cone.num_cons());
    cone.constraint_set(_4ti2_DB, dbs);

    if (lbs.count()+2*dbs.count() <= IndexSetDS::max_size) {
        std::vector<IndexSetDS> supps;
        // Compute ray only constraints first.
        compute(cone, rays, supps, cons_added, ray_ineqs);
        // Compute circuits.
        compute(cone, rays, supps, cons_added, cirs, cir_ineqs);
        // Separate the rays from the circuits.
        split_rays(cone, supps, rays, cirs);
    } else {
        std::vector<IndexSetD> supps;
        // Compute ray only constraints first.
        compute(cone, rays, supps, cons_added, ray_ineqs);
        // Compute circuits.
        compute(cone, rays, supps, cons_added, cirs, cir_ineqs);
        // Separate the rays from the circuits.
        split_rays(cone, supps, rays, cirs);
    }
}

// Computes extreme rays of the cone Ax>=0, x>=0.
template <class T> template <class IndexSet>
void
SupportAlgorithm<T>::compute(
            const ConeT<T>& cone,
            VectorArrayT<T>& rays,
            std::vector<IndexSet>& supps,
            Index& cons_added,
            std::vector<int>& ineqs)
{
    Timer t;
    *out << "Ray Support Algorithm.\n";
    DEBUG_4ti2(*out << "CONSTRAINT MATRIX:\n" << cone.get_matrix() << "\n";)

    // The number of variables.
    Size n = cone.num_vars();
    // The number of constraints.
    Size m = cone.num_cons();

    // The dimension
    Size dim = rays.get_number();

    // The set of constraints to be processed.
    IndexSet rem(n+m, 0);
    cone.constraint_set(_4ti2_LB, rem);

#if 0
    // Construct the initial set supports.
    for (Index i = 0; i != dim; ++i) { 
        IndexSet supp(n+m, 0);
        supp.set(ineqs[i]);
        supps.push_back(supp);
        rem.unset(ineqs[i]);
    }
#endif 
#if 1
    // Construct the initial set supports.
    for (Index i = 0; i != dim; ++i) { 
        IndexSet supp(dim,0);
        supp.set(i);
        supps.push_back(supp);
        rem.unset(ineqs[i]);
    }
#endif
    DEBUG_4ti2(*out << "Initial Rays:\n" << rays << "\n";)
    DEBUG_4ti2(*out << "Initial Supps:\n" << supps << "\n";)

    // The total set of relaxed constraints.
    IndexSet rel(rem);
    cone.constraint_set(_4ti2_FR, rel);
    cone.constraint_set(_4ti2_DB, rel);

    SUPPORTTREE<IndexSet> tree;
    Index next = -1;
    IndexRanges index_ranges;
    IndexSet ray_mask(dim,true);
    // Construct main algorithm object.
    RayState<T,IndexSet> state(cone, rays, supps, rem, ray_mask, next);
    SupportRayAlgorithm<IndexSet> alg(state, supps, rel, cons_added, next, tree, index_ranges);
    // Construct threaded algorithm objects.
    std::vector<SupportRayAlgorithm<IndexSet>*> algs;
    for (Index i = 0; i < Globals::num_threads-1; ++i) { algs.push_back(alg.clone()); }

    // While there are still rows to choose from.
    while (rays.get_number()>0 && !rem.empty()) {
        DEBUG_4ti2(*out << "RAYS:\n" << rays << "\n";)
        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)

        // Choose the next constraint and sort rays.
        Index pos_start, pos_end, neg_start, neg_end;
        next = state.next_constraint(QSolveAlgorithm<T>::order, pos_start, pos_end, neg_start, neg_end);

        // Check to see whether the next constraint is redundant. If so, ignore it.
        if (neg_start == neg_end) { rem.unset(next); continue; }

        char buffer[256];
        sprintf(buffer, "%cLeft %3d  Col %3d  Size %8d", ENDL, rem.count(), next, rays.get_number());
        *out << buffer << std::flush;

        if (pos_start != pos_end) {
            // Next, we assign index ranges for the inner and outer for loops so
            // that the outer loop is smaller than the inner loop.
            Index r1_start, r1_end, r2_start, r2_end;
            if (pos_end-pos_start <= neg_end-neg_start) {
                r1_start = pos_start; r1_end = pos_end;
                r2_start = neg_start; r2_end = neg_end;
            }
            else {
                r2_start = pos_start; r2_end = pos_end;
                r1_start = neg_start; r1_end = neg_end;
            }

            // We sort the r2's into vectors where r2_supp.count()==cons_added+1.
            Index r2_index = state.sort_count(cons_added+1, r2_start, r2_end);

            // Put supports into tree structure.
            tree.insert(supps);
            DEBUG_4ti2(tree.dump();)
            //tree.print_statistics();

            // Run threads.
            index_ranges.init(r1_start, r1_end, r2_start, r2_index, r2_end);
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->threaded_compute(); }
            // Run primary algorithm.
            alg.compute();

            // Wait for threads to finish.
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }

            // Clear tree.
            tree.clear();
        }

        // Delete all the vectors with a negative entry in the column next.
        state.remove(neg_start, neg_end);

        // Add new rays and supps.
        alg.transfer();
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->transfer(); }

        // Update the support vectors for the next_col.
        assert(pos_end <= neg_start); // ASSUMING pos is before neg. 
        //update_supports(supps, next, pos_start, pos_end);
        state.resize(dim+cons_added+1);
        state.update(dim+cons_added, pos_start, pos_end);

        ray_mask.resize(dim+cons_added+1);
        ray_mask.set(dim+cons_added);

        rem.unset(next);
        rel.unset(next);
        ++cons_added;

        // Output statistics.
        sprintf(buffer, "%cLeft %3d  Col %3d  Size %8d  Time %8.2fs", ENDL, rem.count(), next, 
                rays.get_number(), t.get_elapsed_time());
        *out << buffer << std::endl;

        DEBUG_4ti2(*out << "RAYS:\n" << rays << "\n";)
        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)
    }

    // Clean up threaded algorithms objects.
    for (Index i = 0; i < Globals::num_threads-1; ++i) { delete algs[i]; }
}

template <class T> template <class IndexSet>
void
SupportAlgorithm<T>::compute(
        const ConeT<T>& cone,
        VectorArrayT<T>& rays,
        std::vector<IndexSet>& supps,
        Index& cons_added,
        VectorArrayT<T>& cirs,
        std::vector<int>& dbls)
{
    Timer t;

    // The number of variables.
    Size n = cone.num_vars();
    // The number of constraints.
    Size m = cone.num_cons();

    // The remaining columns to process.
    IndexSet rem(n+m, 0); 
    // We only process circuit constraints.
    cone.constraint_set(_4ti2_DB, rem);
    if (rem.empty()) { return; }
    // We have already processed the initial constraints.
    for (Index i = 0; i < (Index) dbls.size(); ++i) { rem.unset(dbls[i]); }

    *out << "Circuit Support Algorithm.\n";

    // The number of circuit components.
    Size dim_cirs = cirs.get_number();

    // Extend the support of the rays to include the circuit components.
    Index ray_supp_size = 0;
    if (!supps.empty()) { ray_supp_size = supps[0].get_size(); }
    if (ray_supp_size%2 == 1) { ++ray_supp_size; } // Make ray support size even.
    Index cir_supp_size = ray_supp_size + 2*dim_cirs;

    // The mask for the circuit constraints.
    IndexSet ray_mask(cir_supp_size, false);
    for (Index i = 0; i < ray_supp_size; ++i) { ray_mask.set(i); }

    Index next = -1;
    RayState<T,IndexSet> state(cone, rays, supps, rem, ray_mask, next);
    state.resize(cir_supp_size);

    // Construct the circuit supports.
    rays.transfer(cirs, 0, cirs.get_number(), rays.get_number());
    for (Index i = 0; i < dim_cirs; ++i) {
        IndexSet supp(cir_supp_size, false);
        supp.set(ray_supp_size+2*i);
        supps.push_back(supp);
    }
    DEBUG_4ti2(*out << "Initial Rays + Circuits:\n" << rays << "\n";)
    DEBUG_4ti2(*out << "Initial Supports:\n" << supps << "\n";)

    SUPPORTTREE<IndexSet> tree;

    IndexRanges index_ranges;
    SupportCirAlgorithm<IndexSet> alg(state, supps, rem, cons_added, next, tree, ray_mask, index_ranges);
    std::vector<SupportCirAlgorithm<IndexSet>*> algs;
    for (Index i = 0; i < Globals::num_threads-1; ++i) { algs.push_back(alg.clone()); }

    while (rays.get_number() > 0 && !rem.empty()) {
        // Find the next constraint.
        int pos_ray_start, pos_ray_end, neg_ray_start, neg_ray_end;
        int pos_cir_start, pos_cir_end, neg_cir_start, neg_cir_end;
        next = state.next_constraint(QSolveAlgorithm<T>::order,
                    pos_ray_start, pos_ray_end, neg_ray_start, neg_ray_end,
                    pos_cir_start, pos_cir_end, neg_cir_start, neg_cir_end);
        DEBUG_4ti2(*out << pos_ray_start << " " << pos_ray_end << " " << neg_ray_start << " " << neg_ray_end << " ";)
        DEBUG_4ti2(*out << pos_cir_start << " " << pos_cir_end << " " << neg_cir_start << " " << neg_cir_end << std::endl;)

        //DEBUG_4ti2(print_debug_diagnostics(cone, rays, supps, next);)

        // Ouput statistics.
        char buffer[256];
        sprintf(buffer, "%cLeft %3d  Col %3d", ENDL, rem.count(), next);
        *out << buffer << std::flush;
        DEBUG_4ti2(*out << "\nRAYS:\n" << rays << "\n";)
        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)

        // Note that the tree needs the ordering of the current vectors to be constant.
        DEBUG_4ti2(*out << "\nBuilding Tree ... " << std::endl;)
        tree.insert(supps);
        // Insert the negatives of the circuit supports.
        for (int i = 0; i < rays.get_number(); ++i) {
            if (ray_mask.set_disjoint(supps[i])) {
                supps[i].swap_odd_n_even();
                tree.insert(supps[i], i);
                supps[i].swap_odd_n_even();
            }
        }
        //tree.print_statistics();
        DEBUG_4ti2(*out << "done." << std::endl;)
        DEBUG_4ti2(tree.dump();)

        // TODO: Make threaded.
        state.flip(pos_cir_start, pos_cir_end);
        alg.compute_cirs(pos_ray_start, neg_cir_end, pos_cir_start, neg_ray_end);
        state.flip(neg_cir_start, neg_cir_end);

        alg.transfer();
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->transfer(); }

        tree.clear();

        // Update the supp vectors for the next_col.
        state.resize(cir_supp_size+2);
        state.update(cir_supp_size, pos_ray_start, pos_ray_end);
        state.update(cir_supp_size, pos_cir_start, pos_cir_end);
        state.update(cir_supp_size+1, neg_ray_start, neg_ray_end);
        state.update(cir_supp_size+1, neg_cir_start, neg_cir_end);
        ray_mask.resize(cir_supp_size+2);
        cir_supp_size+=2;

        sprintf(buffer, "%cLeft %3d  Col %3d  Size %8d  Time %8.2fs", ENDL, rem.count(), next, 
                rays.get_number(), t.get_elapsed_time());
        *out << buffer << std::endl;
        DEBUG_4ti2(*out << "\nRAYS:\n" << rays << "\n";)
        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)

        rem.unset(next);
        ++cons_added;
    }
}

#if 0
template <class T> template <class IndexSet>
void
SupportAlgorithm<T>::check(
                const ConeT<T>& cone,
                const IndexSet& rem,
                const VectorArrayT<T>& rays,
                const std::vector<IndexSet>& supps,
                const std::vector<IndexSet>& cir_supps)
{
    Index n = cone.num_vars();
    Index m = cone.num_cons();
    VectorT<T> slacks(n+m);
    for (Index i = 0; i < rays.get_number(); ++i) {
        cone.get_slacks(rays[i], slacks);
        for (Index j = 0; j < n+m; ++j) {
            if (!rem[j]) {
                if ((slacks[j] != 0) != supps[i][j]) { *out << "Support Check failed.\n"; }
            }
        }
    }
}
#endif


template <class IndexSet>
SupportSubAlgorithmBase<IndexSet>::SupportSubAlgorithmBase(
                RayStateAPI<IndexSet>& _helper, std::vector<IndexSet>& _supps, const IndexSet& _rel, const Index& _cons_added, 
                const Index& _next, const IndexSet& _ray_mask, const SUPPORTTREE<IndexSet>& _tree)
        : helper(*_helper.clone()), supps(_supps), rel(_rel), cons_added(_cons_added), next(_next), ray_mask(_ray_mask), tree(_tree)
{
}

template <class IndexSet>
SupportSubAlgorithmBase<IndexSet>::~SupportSubAlgorithmBase()
{
    delete &helper;
}

#if 0
template <class T, class IndexSet>
void
SupportSubAlgorithmBase<T,IndexSet>::compute(
                Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }
    const std::vector<IndexSet>& local_supps = supps;
    Size supps_size = supps[r1_start].get_size();
    _next = next; //TODO: This is a hack.

    // Temporary variables.
    IndexSet temp_supp(supps_size);
    IndexSet temp_union(supps_size);
    IndexSet zeros(supps_size);
    IndexSet r1_supp(supps_size);
    Size r1_count;
    T s1;

    char buffer[256];
    sprintf(buffer, "  Left = %3d  Col = %3d", rel.count(), next);

    Index index_count = 0;
    // Check negative and positive combinations for adjacency.
    for (Index r1 = r1_start; r1 < r1_end; ++r1) { // Outer loop.
        // Output statistics.
        if (index_count % Globals::output_freq == 0) {
            *out << ENDL << buffer;
            *out << "  Size = " << std::setw(8) << rays.get_number() << ", ";
            *out << "  Index = " << std::setw(8) << r1-r1_start << "/" << r1_end-r1_start;
            *out << std::endl;
        }
        ++index_count;

        r1_supp = local_supps[r1];
        r1_count = r1_supp.count();
        cone.get_slack(rays[r1], next, s1);
        if (r1_count == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                temp_supp.set_difference(local_supps[r2], r1_supp);
                if (temp_supp.singleton()) { helper.create_ray(r2);  }
            }
            continue;
        }

        for (Index r2 = r2_start; r2 < r2_index; ++r2) {
            temp_supp.set_difference(r1_supp, local_supps[r2]);
            if (temp_supp.singleton()) { helper.create_ray(r2);  }
        }

        for (Index r2 = r2_index; r2 < r2_end; ++r2) { // Inner loop.
            // Quick sufficient check whether the two rays are adjacent.
            temp_union.set_union(local_supps[r2], r1_supp);
            if (temp_union.count() > cons_added+2) { continue; }
            // Check whether the two rays r1 and r2 are adjacent.
            temp_supp.set_difference(r1_supp, local_supps[r2]);
            if (temp_supp.singleton()) {
                helper.create_ray(r2); 
                continue;
            }
            temp_supp.set_difference(local_supps[r2], r1_supp);
            if (temp_supp.singleton()) {
                helper.create_ray(r2); 
                continue;
            }
            STATS_4ti2(++num_checks;)
            if (tree.dominated(temp_union, r1, r2)) { STATS_4ti2(++num_dominated;) continue; }
            helper.create_ray(r2); 
        }
    }
}

#endif

#if 1
template <class IndexSet>
void
SupportSubAlgorithmBase<IndexSet>::compute_rays(
                Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }
    const std::vector<IndexSet>& local_supps = supps;
    Size supps_size = supps[r1_start].get_size();

    // Temporary variables.
    IndexSet temp_supp(supps_size);
    IndexSet temp_union(supps_size);
    IndexSet temp_diff(supps_size);
    IndexSet zeros(supps_size);
    IndexSet r1_supp(supps_size);

    std::vector<int> indices;

    char buffer[256];
    sprintf(buffer, "  Left = %3d  Col = %3d", 0, 0);

    SUPPORTTREE<IndexSet> neg_tree;
    neg_tree.insert(local_supps, r2_index, r2_end);
    DEBUG_4ti2(neg_tree.dump());
    std::vector<Index> neg_indices;

    Index index_count = 0;
    // Check negative and positive combinations for adjacency.
    for (Index r1 = r1_start; r1 < r1_end; ++r1) { // Outer loop.
        // Output statistics.
        if (index_count % Globals::output_freq == 0) {
            sprintf(buffer, "%cLeft %3d  Col %3d  Index %8d/%-8d", ENDL, rel.count(), next, index_count, r1_end-r1_start);
            *out << buffer << std::flush;
        }
        ++index_count;

        r1_supp = local_supps[r1];
        helper.set_r1_index(r1);
        if (r1_supp.count() == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (local_supps[r2].singleton_diff(r1_supp)) { helper.create_ray(r2);  }
            }
            continue;
        }

        for (Index r2 = r2_start; r2 < r2_index; ++r2) {
            if (r1_supp.singleton_diff(local_supps[r2])) { helper.create_ray(r2);  }
        }

        zeros.zero();
        indices.clear();
        temp_supp.set_complement(r1_supp);

        // Find the rays whose support differs by one from the current ray's support.
        // These rays are adjacent to r1.
        tree.find_singleton_diff(indices, r1_supp);
        DEBUG_4ti2(std::cout << "Singleton Indices:\n" << indices << "\n";)
        for (unsigned int i = 0; i < indices.size(); ++i) {
            Index r2 = indices[i];
            zeros.set_union(supps[r2]);
            if (r2 >= r2_index && r2 < r2_end) { helper.create_ray(r2); }
        }
        zeros.set_difference(r1_supp);

#if 0
        for (Index r2 = r2_index; r2 < r2_end; ++r2) { // Inner loop.
            if (zeros.set_disjoint(local_supps[r2])) {
                if (r1_supp.count_union(local_supps[r2]) <= cons_added+2) {
                    if (r1_supp.singleton_diff(local_supps[r2])) { helper.create_ray(r2);  continue; }
                    if (local_supps[r2].count_lte_diff(2, r1_supp)) { helper.create_ray(r2); continue; }
                    STATS_4ti2(++num_checks;)
                    temp_union.set_union(local_supps[r2], r1_supp);
                    if (!tree.dominated(temp_union, r1, r2)) { helper.create_ray(r2); }
                }
            }
        }
#endif

#if 1
        neg_indices.clear();
        neg_tree.find(neg_indices, zeros, r1_supp, cons_added+2);
        DEBUG_4ti2(std::cout << "Neg Indices:\n" << neg_indices << "\n";)
        for (Index i = 0; i < (Index) neg_indices.size(); ++i) { // Inner loop.
            Index r2 = neg_indices[i];
            if (r1_supp.singleton_diff(local_supps[r2])) { helper.create_ray(r2);  continue; }
            if (local_supps[r2].count_lte_diff(2, r1_supp)) { helper.create_ray(r2); continue; }
            STATS_4ti2(++num_checks;)
            temp_union.set_union(local_supps[r2], r1_supp);
            //for (Index j = 0; j < (Index) neg_indices.size(); ++j) {
            //    if (local_supps[j].set_subset(temp_union) && i != j) { continue; }
            //}
            if (!tree.dominated(temp_union, r1, r2)) { helper.create_ray(r2); }
        }
#endif
    }
}
#endif

#if 0
template <class IndexSet>
void
SupportSubAlgorithmBase<IndexSet>::compute_rays(
                Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }
    const std::vector<IndexSet>& local_supps = supps;
    Size supps_size = supps[r1_start].get_size();

    // Temporary variables.
    IndexSet r1_supp(supps_size);
    IndexSet temp_union(supps_size);

    SUPPORTTREE<IndexSet> neg_tree;
    neg_tree.insert(local_supps, r2_index, r2_end);
    DEBUG_4ti2(*out << "\nNeg Tree.\n";)
    DEBUG_4ti2(neg_tree.dump());

    SUPPORTTREE<IndexSet> pos_tree;
    for (Index r1 = r1_start; r1 < r1_end; ++r1) { // Outer loop.
        r1_supp = local_supps[r1];
        helper.set_r1_index(r1);
        if (r1_supp.count() == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (local_supps[r2].singleton_diff(r1_supp)) { helper.create_ray(r2);  }
            }
        } else {
            pos_tree.insert(r1_supp, r1);
            for (Index r2 = r2_start; r2 < r2_index; ++r2) {
                if (r1_supp.singleton_diff(local_supps[r2])) { helper.create_ray(r2);  }
            }
        }
    }
    DEBUG_4ti2(*out << "\nPosTree.\n";)
    DEBUG_4ti2(neg_tree.dump());

    std::vector<std::pair<Index,Index> > indices;
    pos_tree.find(indices, neg_tree, cons_added+2);
    for (Index i = 0; i < (Index) indices.size(); ++i) { // Inner loop
        Index r1 = indices[i].first;
        Index r2 = indices[i].second;
        if (local_supps[r2].singleton_diff(local_supps[r1])) { create_ray(r1, r2);  continue; }
        if (local_supps[r1].singleton_diff(local_supps[r2])) { create_ray(r1, r2);  continue; }
        temp_union.set_union(local_supps[r1], local_supps[r2]);
        if (!tree.dominated(temp_union, r1, r2)) { create_ray(r1, r2); }
    }
}
#endif    


#if 0
template <class IndexSet>
void
SupportSubAlgorithmBase<IndexSet>::compute_rays(
                Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }
    const std::vector<IndexSet>& local_supps = supps;
    Size supps_size = supps[r1_start].get_size();

    // Temporary variables.
    IndexSet r1_supp(supps_size);
    IndexSet temp_union(supps_size);
    IndexSet temp_diff(supps_size);
    IndexSet temp_empty(supps_size,0);

    SUPPORTTREE<IndexSet> neg_tree;
    neg_tree.insert(local_supps, r2_index, r2_end);
    DEBUG_4ti2(*out << "\nNeg Tree.\n";)
    DEBUG_4ti2(neg_tree.dump());

    SUPPORTTREE<IndexSet> pos_tree;
    std::vector<IndexSet> diff_supps;
    std::vector<Index> diffs;
    for (Index r1 = r1_start; r1 < r1_end; ++r1) { // Outer loop.
        r1_supp = local_supps[r1];
        helper.set_r1_index(r1);
        if (r1_supp.count() == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (local_supps[r2].singleton_diff(r1_supp)) { helper.create_ray(r2);  }
            }
            diff_supps.push_back(temp_empty);
            continue;
        }

        pos_tree.insert(r1_supp, r1);
        for (Index r2 = r2_start; r2 < r2_index; ++r2) {
            if (r1_supp.singleton_diff(local_supps[r2])) { helper.create_ray(r2);  }
        }

        // Find the rays whose support differs by one from the current ray's support.
        // These rays are adjacent to r1.
        diffs.clear();
        temp_diff.zero();
        tree.find_singleton_diff(diffs, r1_supp);
        DEBUG_4ti2(std::cout << "Singleton Indices:\n" << indices << "\n";)
        for (unsigned int i = 0; i < diffs.size(); ++i) {
            Index r2 = diffs[i];
            temp_diff.set_union(supps[r2]);
            if (r2 >= r2_index && r2 < r2_end) { helper.create_ray(r2); }
        }
        temp_diff.set_difference(r1_supp);
        diff_supps.push_back(temp_diff);
    }
    DEBUG_4ti2(*out << "\nPosTree.\n";)
    DEBUG_4ti2(neg_tree.dump());

    std::vector<std::pair<Index,Index> > indices;
    pos_tree.find(indices, neg_tree, cons_added+2);
    Index r1, r2;
    for (Index i = 0; i < (Index) indices.size(); ++i) { // Inner loop
        if (indices[i].first >= r1_start && indices[i].first < r1_end) { 
            r1 = indices[i].first;
            r2 = indices[i].second;
        } else {
            r2 = indices[i].first;
            r1 = indices[i].second;
        }
        //if (local_supps[r2].singleton_diff(local_supps[r1])) { create_ray(r1, r2);  continue; }
        if (!local_supps[r2].set_disjoint(diff_supps[r1-r1_start])) { continue; }
        if (local_supps[r1].singleton_diff(local_supps[r2])) { create_ray(r1, r2);  continue; }
        temp_union.set_union(local_supps[r1], local_supps[r2]);
        if (!tree.dominated(temp_union, r1, r2)) { create_ray(r1, r2); }
    }
}
#endif    

template <class IndexSet>
void
SupportSubAlgorithmBase<IndexSet>::compute_cirs(Index r1_start, Index r1_end, Index r2_start, Index r2_end)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }
    // TODO: Copy class variables onto the stack???

    char buffer[256];

    DEBUG_4ti2(*out << "\nComputing circuits for ranges ";)
    DEBUG_4ti2(*out << "R1 [" << r1_start << ",...," << r1_end << "] and ";)
    DEBUG_4ti2(*out << "R2 [" << r2_start << ",...," << r2_end << "]\n";)

    Index cir_supp_size = ray_mask.get_size();
    IndexSet temp_supp(cir_supp_size);
    IndexSet cir_mask(ray_mask);
    cir_mask.set_complement();

    IndexSet r1_supp(cir_supp_size);
    IndexSet r1_neg_supp(cir_supp_size);
    Size r1_count;

    int index_count = 0;
    for (int r1 = r1_start; r1 < r1_end; ++r1) {
        if (r2_start <= r1) { r2_start = r1+1; supps[r1].swap_odd_n_even(); }

        r1_supp = supps[r1];
        r1_count = r1_supp.count();
        r1_neg_supp.set_intersection(r1_supp, cir_mask);
        r1_neg_supp.swap_odd_n_even();
        helper.set_r1_index(r1);
        //if (r2_start <= r1) { r2_start = r1+1; IndexSet::swap(r1_supp, r1_neg_supp); }


        if (r1_count == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (r1_neg_supp.set_disjoint(supps[r2])
                    && supps[r2].singleton_diff(r1_supp)) {
                    helper.create_circuit(r2);
                }
            }
            continue;
        }

        for (Index r2 = r2_start; r2 < r2_end; ++r2) {
            if (r1_neg_supp.set_disjoint(supps[r2])
                && r1_supp.count_union(supps[r2]) <= cons_added+2) {
                temp_supp.set_union(r1_supp, supps[r2]);
                if (!tree.dominated(temp_supp, r1, r2)) { helper.create_circuit(r2); }
            }
        }

        if (index_count % Globals::output_freq == 0) {
            sprintf(buffer, "%cLeft %3d  Col %3d  Index %8d/%-8d", ENDL, rel.count(), next, index_count, r1_end-r1_start);
            *out << buffer << std::flush;
        }
        ++index_count;
    }
    sprintf(buffer, "%cLeft %3d  Col %3d  Index %8d/%-8d", ENDL, rel.count(), next, r1_end, r2_end);
    *out << buffer << std::flush;
}

template <class IndexSet>
inline
void
SupportSubAlgorithmBase<IndexSet>::transfer()
{
    //supps.insert(supps.end(), new_supps.begin(), new_supps.end());
    //new_supps.clear();
    helper.transfer();
}

#if 0
template <class IndexSet>
inline
void
SupportSubAlgorithmBase<IndexSet>::set_r1_index(Index _r1)
{
    r1 = _r1;
    helper.set_r1_index(_r1);
}

template <class IndexSet>
inline
void
SupportSubAlgorithmBase<IndexSet>::create_ray(Index r2)
{
    IndexSet temp_supp(supps[r1]);
    temp_supp.set_union(supps[r2]);
    new_supps.push_back(temp_supp);
    helper.create_ray(r2);

    DEBUG_4ti2(
    *out << "\nADDING VECTOR.\n";
    *out << "R1: " << r1 << "\n";
    *out << "R2: " << r2 << "\n";
    *out << temp_supp << "\n";
    )
}

template <class IndexSet>
void
SupportSubAlgorithmBase<IndexSet>::create_circuit(Index i1, Index i2)
{
    DEBUG_4ti2(*out << "Creating new circuit.\n";)
    bool is_pos = helper.create_circuit(i1, i2);

    IndexSet tmp_union(supps[i2]);
    tmp_union.swap_odd_n_even();
    tmp_union.set_union(supps[i1]);
    if (ray_mask.set_disjoint(tmp_union)) {
        if (is_pos) { tmp_union.swap_odd_n_even(); }
    }
    new_supps.push_back(tmp_union);

    DEBUG_4ti2(
        *out << "Cir1 " << supps[i1] << "\n";
        *out << "Cir2 " << supps[i2] << "\n";
        *out << "Cir0 " << tmp_union << "\n";
    )
}
#endif


template <class IndexSet>
SupportRayAlgorithm<IndexSet>::SupportRayAlgorithm(
                RayStateAPI<IndexSet>& _helper, std::vector<IndexSet>& _supps, 
                const IndexSet& _rel, const Index& _cons_added, const Index& _next, const SUPPORTTREE<IndexSet>& _tree,
                IndexRanges& _indices)
        : SupportSubAlgorithmBase<IndexSet>(_helper, _supps, _rel, _cons_added, _next, IndexSet(_rel.get_size(),true), _tree),
          indices(_indices)
{
}

template <class IndexSet>
inline
void
SupportRayAlgorithm<IndexSet>::compute()
{
    Index r1_start, r1_end, r2_start, r2_index, r2_end;
    indices.next(r1_start, r1_end, r2_start, r2_index, r2_end);
    while (r1_start != r1_end) {
        SupportSubAlgorithmBase<IndexSet>::compute_rays(r1_start, r1_end, r2_start, r2_index, r2_end);
        indices.next(r1_start, r1_end, r2_start, r2_index, r2_end);
    }
}

template <class IndexSet>
SupportRayAlgorithm<IndexSet>*
SupportRayAlgorithm<IndexSet>::clone()
{
    return new SupportRayAlgorithm(
                SupportSubAlgorithmBase<IndexSet>::helper,
                SupportSubAlgorithmBase<IndexSet>::supps, 
                SupportSubAlgorithmBase<IndexSet>::rel,
                SupportSubAlgorithmBase<IndexSet>::cons_added,
                SupportSubAlgorithmBase<IndexSet>::next,
                SupportSubAlgorithmBase<IndexSet>::tree,
                indices);
}

template <class IndexSet>
SupportCirAlgorithm<IndexSet>::SupportCirAlgorithm(
                RayStateAPI<IndexSet>& _helper, std::vector<IndexSet>& _supps, 
                const IndexSet& _rel, const Index& _cons_added, const Index& _next, const SUPPORTTREE<IndexSet>& _tree,
                const IndexSet& _ray_mask, IndexRanges& _indices)
        : SupportSubAlgorithmBase<IndexSet>(_helper, _supps, _rel, _cons_added, _next, _ray_mask, _tree),
          indices(_indices)
{
}

template <class IndexSet>
inline
void
SupportCirAlgorithm<IndexSet>::compute()
{
    Index r1_start, r1_end, r2_start, r2_index, r2_end;
    indices.next(r1_start, r1_end, r2_start, r2_index, r2_end);
    while (r1_start != r1_end) {
        SupportSubAlgorithmBase<IndexSet>::compute_cirs(r1_start, r1_end, r2_start, r2_end);
        indices.next(r1_start, r1_end, r2_start, r2_index, r2_end);
    }
}

template <class IndexSet>
SupportCirAlgorithm<IndexSet>*
SupportCirAlgorithm<IndexSet>::clone()
{
    return new SupportCirAlgorithm(
                SupportSubAlgorithmBase<IndexSet>::helper,
                SupportSubAlgorithmBase<IndexSet>::supps, 
                SupportSubAlgorithmBase<IndexSet>::rel,
                SupportSubAlgorithmBase<IndexSet>::cons_added,
                SupportSubAlgorithmBase<IndexSet>::next,
                SupportSubAlgorithmBase<IndexSet>::tree,
                SupportSubAlgorithmBase<IndexSet>::ray_mask,
                indices);
}

#undef DEBUG_4ti2
