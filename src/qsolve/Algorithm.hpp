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

#include "qsolve/Algorithm.h"
#include "qsolve/Globals.h"
#include "qsolve/Timer.h"
#include "qsolve/Debug.h"
#include "qsolve/Cone.h"

#include "qsolve/VectorStream.h"
#include "qsolve/VectorArrayStream.h"
#include "qsolve/IndexSetStream.h"
#include "qsolve/Stream.h"
#include "qsolve/MatrixAlgorithm.h"
#include "qsolve/SupportAlgorithm.h"

#undef DEBUG_4ti2
#define DEBUG_4ti2(X) //X

using namespace _4ti2_;

template <class IndexSet>
Algorithm<IndexSet>::Algorithm()
{
}

template <class IndexSet>
Algorithm<IndexSet>::Algorithm(ConsOrder o, QSolveVariant v)
    : order(o), variant(v)
{
}

template <class IndexSet>
Algorithm<IndexSet>::~Algorithm()
{
}

// Computes extreme rays of the cone Ax>=0, x>=0.
template <class IndexSet>
void
Algorithm<IndexSet>::compute_rays(
            const ConeAPI& cone,
            RayStateAPI<IndexSet>& state,
            std::vector<int>& ineqs)
{
    Timer t;
    if (variant == SUPPORT) { *out << "Ray Support Algorithm.\n"; }
    else  { *out << "Ray Matrix Algorithm.\n"; }

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
    SubAlgorithm* alg = 0;
    if (variant == SUPPORT) { alg = new SupportRayAlgorithm<IndexSet>(state, index_ranges); }
    else { alg = new MatrixRayAlgorithm<IndexSet>(state, index_ranges); }

    // Construct threaded algorithm objects.
    std::vector<SubAlgorithm*> algs;
    for (Index i = 0; i < Globals::num_threads-1; ++i) { algs.push_back(alg->clone()); }

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

            if (variant == SUPPORT) {
                // Put supports into tree structure.
                state.tree.insert(supps);
                DEBUG_4ti2(state.tree.dump();)
                DEBUG_4ti2(tree.print_statistics();)
            }

            // Run threads.
            index_ranges.init(r1_start, r1_end, r2_start, r2_index, r2_end);
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->threaded_compute(); }
            // Run primary algorithm.
            alg->compute();

            // Wait for threads to finish.
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }

            // Clear tree.
            if (variant == SUPPORT) { state.tree.clear(); }
        }

        // Delete all the vectors with a negative entry in the column next.
        state.remove(neg_start, neg_end);

        // Add new rays and supps.
        alg->transfer();
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
    delete alg;
}

template <class IndexSet>
void
Algorithm<IndexSet>::compute_cirs(
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

    if (variant == SUPPORT) { *out << "Circuit Support Algorithm.\n"; }
    else  { *out << "Circuit Matrix Algorithm.\n"; }

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
    SubAlgorithm* alg = 0;
    if (variant == SUPPORT) { alg = new SupportCirAlgorithm<IndexSet>(state, index_ranges); }
    else { alg = new MatrixCirAlgorithm<IndexSet>(state, index_ranges); }
    std::vector<SubAlgorithm*> algs;
    for (Index i = 0; i < Globals::num_threads-1; ++i) { algs.push_back(alg->clone()); }

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

        // We sort the r2's into vectors where r2_supp.count()==cons_added+1.
        Index pos_cir_middle = state.sort_count(cons_added+1, pos_cir_start, pos_cir_end);

        if (variant == SUPPORT) {
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
            DEBUG_4ti2(*out << "done." << std::endl;)
            DEBUG_4ti2(state.ree.dump();)
            DEBUG_4ti2(state.tree.print_statistics();)
        }

        // Flip the positive circuits.
        state.flip(pos_cir_start, pos_cir_end);
        // Run threads.
        index_ranges.init(pos_ray_start, neg_cir_end, pos_cir_start, pos_cir_middle, neg_ray_end);
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->threaded_compute(); }
        // Run primary algorithm.
        alg->compute();

        // Wait for threads to finish.
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }

        // Unflip the positive circuits.
        state.flip(pos_cir_start, pos_cir_end);

        alg->transfer();
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->transfer(); }

        if (variant == SUPPORT) { state.tree.clear(); }

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

    // Clean up threaded algorithms objects.
    for (Index i = 0; i < Globals::num_threads-1; ++i) { delete algs[i]; }
    delete alg;
}


#undef DEBUG_4ti2
