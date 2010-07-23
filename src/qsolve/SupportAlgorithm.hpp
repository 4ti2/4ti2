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
#undef STATS_4ti2
#define STATS_4ti2(X) //X

using namespace _4ti2_;

template <class IndexSet>
SupportAlgorithm<IndexSet>::SupportAlgorithm()
{
}

template <class IndexSet>
SupportAlgorithm<IndexSet>::SupportAlgorithm(ConsOrder o)
    : order(o)
{
}

template <class IndexSet>
SupportAlgorithm<IndexSet>::~SupportAlgorithm()
{
}

// Computes extreme rays of the cone Ax>=0, x>=0.
template <class IndexSet>
void
SupportAlgorithm<IndexSet>::compute_rays(
            const ConeAPI& cone,
            RayStateAPI<IndexSet>& state,
            std::vector<int>& ineqs)
{
    Timer t;
    *out << "Ray Support Algorithm.\n";
    DEBUG_4ti2(*out << "CONSTRAINT MATRIX:\n" << cone.get_matrix() << "\n";)

    std::vector<IndexSet>& supps = state.supps;
    Index& next = state.next;
    Index& cons_added = state.cons_added;
    IndexSet& ray_mask = state.ray_mask;

    // The number of variables.
    Size n = cone.num_vars();
    // The number of constraints.
    Size m = cone.num_cons();
    // The dimension
    Size dim = state.num_gens();

    // The set of constraints to be processed.
    IndexSet& rem = state.rem;
    cone.get_constraint_set(_4ti2_LB, rem);

    // The total set of relaxed constraints.
    IndexSet& rel = state.rel;
    cone.get_constraint_set(_4ti2_EQ, rel);
    rel.set_complement();

    IndexRanges index_ranges;
    ray_mask.resize(dim);
    ray_mask.one();

    // Construct the initial set supports.
    state.supp_types.clear();
    state.supp_types.resize(dim, _4ti2_LB);
    state.supps_to_cons = ineqs;
    state.cons_to_supps.clear();
    state.cons_to_supps.resize(n+m,-1);
    for (Index i = 0; i != dim; ++i) { 
        IndexSet supp(dim,0);
        supp.set(i);
        supps.push_back(supp);
        rem.unset(ineqs[i]);
        rel.unset(ineqs[i]);
        state.cons_to_supps[ineqs[i]] = i;
    }

    DEBUG_4ti2(*out << "Initial Supps:\n" << supps << "\n";)

    // Construct main algorithm object.
    SupportRayAlgorithm<IndexSet> alg(state, index_ranges);
    // Construct threaded algorithm objects.
    std::vector<SupportRayAlgorithm<IndexSet>*> algs;
    for (Index i = 0; i < Globals::num_threads-1; ++i) { algs.push_back(alg.clone()); }

    // While there are still rows to choose from.
    while (state.num_gens() > 0 && !rem.empty()) {
        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)

        // Choose the next constraint and sort rays.
        Index pos_start, pos_end, neg_start, neg_end;
        next = state.next_constraint(order, rem, pos_start, pos_end, neg_start, neg_end);

        // Check to see whether the next constraint is redundant. If so, ignore it.
        // TODO: We should check this before sorting the supports and rays.
        if (neg_start == neg_end) { rem.unset(next); continue; }

        char buffer[256];
        sprintf(buffer, "%cLeft %3d  Col %3d  Size %8d", ENDL, rem.count(), next, state.num_gens());
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
            state.tree.insert(supps);
            DEBUG_4ti2(state.tree.dump();)
            //tree.print_statistics();

            // Run threads.
            index_ranges.init(r1_start, r1_end, r2_start, r2_index, r2_end);
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->threaded_compute(); }
            // Run primary algorithm.
            alg.compute();

            // Wait for threads to finish.
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }

            // Clear tree.
            state.tree.clear();
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

        state.supps_to_cons.push_back(next);
        state.cons_to_supps[next] = dim+cons_added;
        state.supp_types.push_back(_4ti2_LB);

        rem.unset(next);
        rel.unset(next);
        ++cons_added;

        // Output statistics.
        sprintf(buffer, "%cLeft %3d  Col %3d  Size %8d  Time %8.2fs", ENDL, rem.count(), next, 
                state.num_gens(), t.get_elapsed_time());
        *out << buffer << std::endl;

        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)
        DEBUG_4ti2(state.check());
    }

    // Clean up threaded algorithms objects.
    for (Index i = 0; i < Globals::num_threads-1; ++i) { delete algs[i]; }
}

template <class IndexSet>
void
SupportAlgorithm<IndexSet>::compute_cirs(
        const ConeAPI& cone,
        RayStateAPI<IndexSet>& state,
        std::vector<Index>& dbls)
{
    Timer t;

    std::vector<IndexSet>& supps = state.supps;
    Index& next = state.next;
    Index& cons_added = state.cons_added;
    IndexSet& ray_mask = state.ray_mask;
    IndexSet& rem = state.rem;
    IndexSet& rel = state.rel;

    // We only process circuit constraints.
    rem.zero();
    cone.get_constraint_set(_4ti2_DB, rem);
    if (rem.empty()) { return; }
    // We have already processed the initial constraints.
    for (Index i = 0; i < (Index) dbls.size(); ++i) { rem.unset(dbls[i]); rel.unset(dbls[i]); }

    *out << "Circuit Support Algorithm.\n";

    // The number of circuit components.
    Size dim_cirs = dbls.size();

    // Extend the support of the rays to include the circuit components.
    Index ray_supp_size = 0;
    if (!supps.empty()) { ray_supp_size = supps[0].get_size(); }
    // Make ray support size even.
    if (ray_supp_size%2 == 1) { 
        ++ray_supp_size; 
        state.supps_to_cons.push_back(-1); 
        state.supp_types.push_back(_4ti2_FR);
    }
    Index cir_supp_size = ray_supp_size + 2*dim_cirs;

    // The mask for the circuit constraints.
    ray_mask.resize(cir_supp_size);
    ray_mask.zero();
    for (Index i = 0; i < ray_supp_size; ++i) { ray_mask.set(i); }
    state.resize(cir_supp_size);

    // Construct the circuit supports.
    for (Index i = 0; i < dim_cirs; ++i) {
        IndexSet supp(cir_supp_size, false);
        supp.set(ray_supp_size+2*i);
        supps.push_back(supp);

        state.supps_to_cons.push_back(dbls[i]);
        state.supps_to_cons.push_back(dbls[i]);
        state.supp_types.push_back(_4ti2_LB);
        state.supp_types.push_back(_4ti2_UB);
        state.cons_to_supps[dbls[i]] = ray_supp_size + 2*i;
    }
    DEBUG_4ti2(*out << "Initial Rays + Circuits:\n" << rays << "\n";)
    DEBUG_4ti2(*out << "Initial Supports:\n" << supps << "\n";)

    IndexRanges index_ranges;
    SupportCirAlgorithm<IndexSet> alg(state, index_ranges);
    std::vector<SupportCirAlgorithm<IndexSet>*> algs;
    for (Index i = 0; i < Globals::num_threads-1; ++i) { algs.push_back(alg.clone()); }

    while (state.num_gens() > 0 && !rem.empty()) {
        // Find the next constraint.
        int pos_ray_start, pos_ray_end, neg_ray_start, neg_ray_end;
        int pos_cir_start, pos_cir_end, neg_cir_start, neg_cir_end;
        next = state.next_constraint(order, rem, ray_mask,
                    pos_ray_start, pos_ray_end, neg_ray_start, neg_ray_end,
                    pos_cir_start, pos_cir_end, neg_cir_start, neg_cir_end);

        // Ouput statistics.
        char buffer[256];
        sprintf(buffer, "%cLeft %3d  Col %3d", ENDL, rem.count(), next);
        *out << buffer << std::flush;
        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)

        // Note that the tree needs the ordering of the current vectors to be constant.
        DEBUG_4ti2(*out << "\nBuilding Tree ... " << std::endl;)
        state.tree.insert(supps);
        // Insert the negatives of the circuit supports.
        for (int i = 0; i < state.num_gens(); ++i) {
            if (ray_mask.set_disjoint(supps[i])) {
                supps[i].swap_odd_n_even();
                state.tree.insert(supps[i], i);
                supps[i].swap_odd_n_even();
            }
        }
        //state.tree.print_statistics();
        DEBUG_4ti2(*out << "done." << std::endl;)
        DEBUG_4ti2(state.ree.dump();)

        // TODO: Make threaded. // Change algorithm so that supps are constant.
        state.flip(pos_cir_start, pos_cir_end);
        alg.compute_cirs(pos_ray_start, neg_cir_end, pos_cir_start, neg_ray_end);
        state.flip(neg_cir_start, neg_cir_end);

        alg.transfer();
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->transfer(); }

        state.tree.clear();

        // Update the supp vectors for the next_col.
        state.resize(cir_supp_size+2);
        state.update(cir_supp_size, pos_ray_start, pos_ray_end);
        state.update(cir_supp_size, pos_cir_start, pos_cir_end);
        state.update(cir_supp_size+1, neg_ray_start, neg_ray_end);
        state.update(cir_supp_size+1, neg_cir_start, neg_cir_end);
        ray_mask.resize(cir_supp_size+2);

        state.supps_to_cons.push_back(next);
        state.supps_to_cons.push_back(next);
        state.supp_types.push_back(_4ti2_LB);
        state.supp_types.push_back(_4ti2_UB);
        state.cons_to_supps[next] = cir_supp_size;

        cir_supp_size+=2;

        sprintf(buffer, "%cLeft %3d  Col %3d  Size %8d  Time %8.2fs", ENDL, rem.count(), next, 
                state.num_gens(), t.get_elapsed_time());
        *out << buffer << std::endl;

        rem.unset(next);
        rel.unset(next);
        ++cons_added;

        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)
        DEBUG_4ti2(state.check());
    }
}

template <class IndexSet>
SupportSubAlgorithmBase<IndexSet>::SupportSubAlgorithmBase(RayStateAPI<IndexSet>& _state)
        : state(_state), helper(*_state.clone())
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
    const std::vector<IndexSet>& supps = supps;
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
    sprintf(buffer, "  Left = %3d  Col = %3d", state.rem.count(), state.next);

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

        r1_supp = supps[r1];
        r1_count = r1_supp.count();
        cone.get_slack(rays[r1], next, s1);
        if (r1_count == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                temp_supp.set_difference(supps[r2], r1_supp);
                if (temp_supp.singleton()) { helper.create_ray(r2);  }
            }
            continue;
        }

        for (Index r2 = r2_start; r2 < r2_index; ++r2) {
            temp_supp.set_difference(r1_supp, supps[r2]);
            if (temp_supp.singleton()) { helper.create_ray(r2);  }
        }

        for (Index r2 = r2_index; r2 < r2_end; ++r2) { // Inner loop.
            // Quick sufficient check whether the two rays are adjacent.
            temp_union.set_union(supps[r2], r1_supp);
            if (temp_union.count() > cons_added+2) { continue; }
            // Check whether the two rays r1 and r2 are adjacent.
            temp_supp.set_difference(r1_supp, supps[r2]);
            if (temp_supp.singleton()) {
                helper.create_ray(r2); 
                continue;
            }
            temp_supp.set_difference(supps[r2], r1_supp);
            if (temp_supp.singleton()) {
                helper.create_ray(r2); 
                continue;
            }
            STATS_4ti2(++num_checks;)
            if (state.tree.dominated(temp_union, r1, r2)) { STATS_4ti2(++num_dominated;) continue; }
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
    const std::vector<IndexSet>& supps = state.supps;
    SUPPORTTREE<IndexSet>& tree = state.tree;

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
    neg_tree.insert(supps, r2_index, r2_end);
    DEBUG_4ti2(neg_tree.dump());
    std::vector<Index> neg_indices;

    Index index_count = 0;
    // Check negative and positive combinations for adjacency.
    for (Index r1 = r1_start; r1 < r1_end; ++r1) { // Outer loop.
        // Output statistics.
        if (index_count % Globals::output_freq == 0) {
            sprintf(buffer, "%cLeft %3d  Col %3d  Index %8d/%-8d", ENDL, state.rem.count(), state.next, index_count, r1_end-r1_start);
            *out << buffer << std::flush;
        }
        ++index_count;

        r1_supp = supps[r1];
        helper.set_r1_index(r1);
        if (r1_supp.count() == state.cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (supps[r2].singleton_diff(r1_supp)) { helper.create_ray(r2);  }
            }
            continue;
        }

        for (Index r2 = r2_start; r2 < r2_index; ++r2) {
            if (r1_supp.singleton_diff(supps[r2])) { helper.create_ray(r2);  }
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
            if (zeros.set_disjoint(supps[r2])) {
                if (r1_supp.count_union(supps[r2]) <= cons_added+2) {
                    if (r1_supp.singleton_diff(supps[r2])) { helper.create_ray(r2);  continue; }
                    if (supps[r2].count_lte_diff(2, r1_supp)) { helper.create_ray(r2); continue; }
                    STATS_4ti2(++num_checks;)
                    temp_union.set_union(supps[r2], r1_supp);
                    if (!tree.dominated(temp_union, r1, r2)) { helper.create_ray(r2); }
                }
            }
        }
#endif

#if 1
        neg_indices.clear();
        neg_tree.find(neg_indices, zeros, r1_supp, state.cons_added+2);
        DEBUG_4ti2(std::cout << "Neg Indices:\n" << neg_indices << "\n";)
        for (Index i = 0; i < (Index) neg_indices.size(); ++i) { // Inner loop.
            Index r2 = neg_indices[i];
            if (r1_supp.singleton_diff(supps[r2])) { helper.create_ray(r2);  continue; }
            if (supps[r2].count_lte_diff(2, r1_supp)) { helper.create_ray(r2); continue; }
            STATS_4ti2(++num_checks;)
            temp_union.set_union(supps[r2], r1_supp);
            //for (Index j = 0; j < (Index) neg_indices.size(); ++j) {
            //    if (supps[j].set_subset(temp_union) && i != j) { continue; }
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
    const std::vector<IndexSet>& supps = state.supps;
    SUPPORTTREE<IndexSet>& tree = state.tree;
    Size supps_size = supps[r1_start].get_size();

    // Temporary variables.
    IndexSet r1_supp(supps_size);
    IndexSet temp_union(supps_size);

    SUPPORTTREE<IndexSet> neg_tree;
    neg_tree.insert(supps, r2_index, r2_end);
    DEBUG_4ti2(*out << "\nNeg Tree.\n";)
    DEBUG_4ti2(neg_tree.dump());

    SUPPORTTREE<IndexSet> pos_tree;
    for (Index r1 = r1_start; r1 < r1_end; ++r1) { // Outer loop.
        r1_supp = supps[r1];
        helper.set_r1_index(r1);
        if (r1_supp.count() == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (supps[r2].singleton_diff(r1_supp)) { helper.create_ray(r2);  }
            }
        } else {
            pos_tree.insert(r1_supp, r1);
            for (Index r2 = r2_start; r2 < r2_index; ++r2) {
                if (r1_supp.singleton_diff(supps[r2])) { helper.create_ray(r2);  }
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
        if (supps[r2].singleton_diff(supps[r1])) { create_ray(r1, r2);  continue; }
        if (supps[r1].singleton_diff(supps[r2])) { create_ray(r1, r2);  continue; }
        temp_union.set_union(supps[r1], supps[r2]);
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
    const std::vector<IndexSet>& supps = state.supps;
    SUPPORTTREE<IndexSet>& tree = state.tree;
    Size supps_size = supps[r1_start].get_size();

    // Temporary variables.
    IndexSet r1_supp(supps_size);
    IndexSet temp_union(supps_size);
    IndexSet temp_diff(supps_size);
    IndexSet temp_empty(supps_size,0);

    SUPPORTTREE<IndexSet> neg_tree;
    neg_tree.insert(supps, r2_index, r2_end);
    DEBUG_4ti2(*out << "\nNeg Tree.\n";)
    DEBUG_4ti2(neg_tree.dump());

    SUPPORTTREE<IndexSet> pos_tree;
    std::vector<IndexSet> diff_supps;
    std::vector<Index> diffs;
    for (Index r1 = r1_start; r1 < r1_end; ++r1) { // Outer loop.
        r1_supp = supps[r1];
        helper.set_r1_index(r1);
        if (r1_supp.count() == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (supps[r2].singleton_diff(r1_supp)) { helper.create_ray(r2);  }
            }
            diff_supps.push_back(temp_empty);
            continue;
        }

        pos_tree.insert(r1_supp, r1);
        for (Index r2 = r2_start; r2 < r2_index; ++r2) {
            if (r1_supp.singleton_diff(supps[r2])) { helper.create_ray(r2);  }
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
        //if (supps[r2].singleton_diff(supps[r1])) { create_ray(r1, r2);  continue; }
        if (!supps[r2].set_disjoint(diff_supps[r1-r1_start])) { continue; }
        if (supps[r1].singleton_diff(supps[r2])) { create_ray(r1, r2);  continue; }
        temp_union.set_union(supps[r1], supps[r2]);
        if (!tree.dominated(temp_union, r1, r2)) { create_ray(r1, r2); }
    }
}
#endif    

template <class IndexSet>
void
SupportSubAlgorithmBase<IndexSet>::compute_cirs(Index r1_start, Index r1_end, Index r2_start, Index r2_end)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }
    std::vector<IndexSet>& supps = state.supps;
    SUPPORTTREE<IndexSet>& tree = state.tree;

    char buffer[256];

    DEBUG_4ti2(*out << "\nComputing circuits for ranges ";)
    DEBUG_4ti2(*out << "R1 [" << r1_start << ",...," << r1_end << "] and ";)
    DEBUG_4ti2(*out << "R2 [" << r2_start << ",...," << r2_end << "]\n";)

    Index cir_supp_size = state.ray_mask.get_size();
    IndexSet temp_supp(cir_supp_size);
    IndexSet cir_mask(state.ray_mask);
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

        if (r1_count == state.cons_added+1) {
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
                && r1_supp.count_union(supps[r2]) <= state.cons_added+2) {
                temp_supp.set_union(r1_supp, supps[r2]);
                if (!tree.dominated(temp_supp, r1, r2)) { helper.create_circuit(r2); }
            }
        }

        if (index_count % Globals::output_freq == 0) {
            sprintf(buffer, "%cLeft %3d  Col %3d  Index %8d/%-8d", ENDL, state.rem.count(), state.next, index_count, r1_end-r1_start);
            *out << buffer << std::flush;
        }
        ++index_count;
    }
    sprintf(buffer, "%cLeft %3d  Col %3d  Index %8d/%-8d", ENDL, state.rem.count(), state.next, r1_end, r2_end);
    *out << buffer << std::flush;
}

template <class IndexSet>
inline
void
SupportSubAlgorithmBase<IndexSet>::transfer()
{
    helper.transfer();
}

template <class IndexSet>
SupportRayAlgorithm<IndexSet>::SupportRayAlgorithm(RayStateAPI<IndexSet>& _state, IndexRanges& _indices)
        : SupportSubAlgorithmBase<IndexSet>(_state), indices(_indices)
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
    return new SupportRayAlgorithm(SupportSubAlgorithmBase<IndexSet>::state, indices);
}

template <class IndexSet>
SupportCirAlgorithm<IndexSet>::SupportCirAlgorithm(RayStateAPI<IndexSet>& _state, IndexRanges& _indices)
        : SupportSubAlgorithmBase<IndexSet>(_state), indices(_indices)
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
    return new SupportCirAlgorithm(SupportSubAlgorithmBase<IndexSet>::state, indices);
}

#undef DEBUG_4ti2
