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
    helper.r1_supp.resize(supps_size);
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
    helper.r1_supp.resize(supps_size);

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
    helper.r1_supp.resize(supps_size);

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
    helper.r1_supp.resize(supps_size);

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
    helper.r1_supp.resize(cir_supp_size);
    Size r1_count;

    int index_count = 0;
    for (int r1 = r1_start; r1 < r1_end; ++r1) {
        r1_supp = supps[r1];
        r1_count = r1_supp.count();
        r1_neg_supp.set_intersection(r1_supp, cir_mask);
        r1_neg_supp.swap_odd_n_even();
        helper.set_r1_index(r1);
        if (r2_start <= r1) { r2_start = r1+1; IndexSet::swap(r1_supp, r1_neg_supp); helper.r1_supp = r1_supp; }

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
