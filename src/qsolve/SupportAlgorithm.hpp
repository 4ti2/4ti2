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

#include "qsolve/SupportAlgorithm.h"
#include "qsolve/DiagonalAlgorithm.h"
#include "qsolve/HermiteAlgorithm.h"
#include "qsolve/Euclidean.h"
#include "qsolve/Globals.h"
#include "qsolve/Timer.h"
#include "qsolve/Debug.h"
#include "qsolve/Cone.h"

#include "qsolve/VectorStream.h"
#include "qsolve/VectorArrayStream.h"
#include "qsolve/IndexSetStream.h"
#include "qsolve/Stream.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

#undef DEBUG_4ti2
#define DEBUG_4ti2(X) //X

using namespace _4ti2_;

template <class T, class IndexSet>
SupportAlgorithm<T,IndexSet>::SupportAlgorithm()
    : QSolveAlgorithm<T,IndexSet>()
{
}

template <class T, class IndexSet>
SupportAlgorithm<T,IndexSet>::SupportAlgorithm(QSolveConsOrder o)
    : QSolveAlgorithm<T,IndexSet>(o)
{
}

template <class T, class IndexSet>
SupportAlgorithm<T,IndexSet>::~SupportAlgorithm()
{
}

template <class T, class IndexSet>
void
SupportAlgorithm<T,IndexSet>::compute(const ConeT<T>& cone, 
            VectorArrayT<T>& rays, std::vector<Index>& ray_ineqs, 
            VectorArrayT<T>& cirs, std::vector<Index>& cir_ineqs)
{
    // The number of constraints added so far.
    Index cons_added = 0;
    std::vector<IndexSet> supps;
    // Compute ray only constraints first.
    compute(cone, rays, supps, cons_added, ray_ineqs);

    // Compute circuits.
    compute(cone, rays, supps, cons_added, cirs, cir_ineqs);

    // Separate the rays from the circuits.
    split_rays(cone, supps, rays, cirs);
}

// Computes extreme rays of the cone Ax>=0, x>=0.
template <class T, class IndexSet>
void
SupportAlgorithm<T,IndexSet>::compute(
            const ConeT<T>& cone,
            VectorArrayT<T>& rays,
            std::vector<IndexSet>& supps,
            Index& cons_added,
            std::vector<int>& ineqs)
{
    Timer t;
    *out << "Ray Support Algorithm.\n";
    DEBUG_4ti2(*out << "CONSTRAINT MATRIX:\n" << cone.constraint_matrix() << "\n";)

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

    // Construct main algorithm object.
    SupportRayAlgorithm<T,IndexSet> alg(cone, rays, supps, rel, cons_added, next, tree);
    // Construct threaded algorithm objects.
    std::vector<SupportRayAlgorithm<T,IndexSet>*> algs;
    for (Index i = 0; i < Globals::num_threads-1; ++i) {
            algs.push_back(new SupportRayAlgorithm<T,IndexSet>(cone, rays, supps, rel, cons_added, next, tree));
    }

    // While there are still rows to choose from.
    while (rays.get_number()>0 && !rem.empty()) {
        DEBUG_4ti2(*out << "RAYS:\n" << rays << "\n";)
        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)

        // Choose the next constraint and sort rays.
        Index pos_start, pos_end, neg_start, neg_end;
        next = next_constraint(cone, rem, rays, supps,
                        pos_start, pos_end, neg_start, neg_end);

        char buffer[256];
        sprintf(buffer, "  Left = %3d  Col = %3d", rem.count(), next);
        *out << ENDL << buffer;
        *out << "  Size = " << std::setw(8) << rays.get_number() << "  Time: " << t;

        // Check to see whether the next constraint is redundant. If so, ignore it.
        if (neg_start == neg_end) { rem.unset(next); continue; }

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
            Index r2_index = r2_start;
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (supps[r2].count() == cons_added+1) {
                    rays.swap_rows(r2, r2_index);
                    IndexSet::swap(supps[r2], supps[r2_index]);
                    ++r2_index;
                }
            }
    
            // Put supports into tree structure.
            tree.insert(supps);
            //tree.print_statistics();
    
            // Run threads.
            Index increment = (r1_end-r1_start)/Globals::num_threads;
            for (Index i = 0; i < Globals::num_threads-1; ++i) {
                algs[i]->threaded_compute(r1_start, r1_start+increment, r2_start, r2_index, r2_end);
                r1_start += increment;
            }
            // Run primary algorithm.
            alg.compute(r1_start, r1_end, r2_start, r2_index, r2_end);
    
            // Wait for threads to finish.
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }
    
            // Clear tree.
            tree.clear();
        }

        // Delete all the vectors with a negative entry in the column next.
        supps.erase(supps.begin()+neg_start, supps.begin()+neg_end);
        rays.remove(neg_start, neg_end);

        // Add new rays and supps.
        alg.transfer(rays, supps);
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->transfer(rays, supps); }

        // Update the support vectors for the next_col.
        assert(pos_end <= neg_start); // ASSUMING pos is before neg. 
        //update_supports(supps, next, pos_start, pos_end);
        resize_supports(supps, dim+cons_added+1);
        update_supports(supps, dim+cons_added, pos_start, pos_end);

        rem.unset(next);
        rel.unset(next);
        ++cons_added;

        // Output statistics.
        *out << ENDL << buffer;
        *out << "  Size = " << std::setw(8) << rays.get_number() << ", ";
        *out << "  Time: " << t << "                \n";

        DEBUG_4ti2(*out << "RAYS:\n" << rays << "\n";)
        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)
    }

    // Clean up threaded algorithms objects.
    for (Index i = 0; i < Globals::num_threads-1; ++i) { delete algs[i]; }
}

template <class T, class IndexSet>
void
SupportAlgorithm<T,IndexSet>::compute(
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
    if (rays.get_number() != 0) { ray_supp_size = supps[0].get_size(); }
    if (ray_supp_size%2 == 1) { ++ray_supp_size; } // Make ray support size even.
    Index cir_supp_size = ray_supp_size + 2*dim_cirs;
    for (Index i = 0; i < (Index) supps.size(); ++i) { supps[i].resize(cir_supp_size); }

    // Construct the circuit supports.
    rays.transfer(cirs, 0, cirs.get_number(), rays.get_number());
    for (Index i = 0; i < dim_cirs; ++i) {
        IndexSet supp(cir_supp_size, false);
        supp.set(ray_supp_size+2*i);
        supps.push_back(supp);
    }

    // The mask for the circuit constraints.
    IndexSet ray_mask(cir_supp_size, false);
    for (Index i = 0; i < ray_supp_size; ++i) { ray_mask.set(i); }

    SUPPORTTREE<IndexSet> tree;
    Index next = -1;

    SupportCirAlgorithm<T,IndexSet> alg(cone, rays, supps, rem, cons_added, next, tree, ray_mask);

    while (rays.get_number() > 0 && !rem.empty()) {
        // Find the next constraint.
        int ray_start, ray_end, cir_start, cir_end;
        int pos_ray_start, pos_ray_end, neg_ray_start, neg_ray_end;
        int pos_cir_start, pos_cir_end, neg_cir_start, neg_cir_end;
        Index next_col = next_constraint(cone, rays, supps, rem, ray_mask,
                    ray_start, ray_end, cir_start, cir_end,
                    pos_ray_start, pos_ray_end, neg_ray_start, neg_ray_end,
                    pos_cir_start, pos_cir_end, neg_cir_start, neg_cir_end);

        DEBUG_4ti2(print_debug_diagnostics(cone, rays, supps, next_col);)

        // Ouput statistics.
        char buffer[256];
        sprintf(buffer, "  Left = %3d  Col = %3d", rem.count(), next_col);
        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << rays.get_number() << "  Time: " << t;
        DEBUG_4ti2(*out << std::endl;)

        // Note that the tree needs the ordering of the current vectors to be
        // constant.
        DEBUG_4ti2(*out << "Building Tree ... " << std::flush;)
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
        if (pos_ray_start != pos_ray_end) {
            flip(supps, pos_cir_start, pos_cir_end);
            alg.compute(pos_ray_start, pos_ray_end, neg_ray_start, cir_end);
            flip(supps, pos_cir_start, pos_cir_end);
        }
        if (neg_ray_start != neg_ray_end) {
            flip(supps, neg_cir_start, neg_cir_end);
            alg.compute(neg_ray_start, neg_ray_end, cir_start, cir_end);
            flip(supps, neg_cir_start, neg_cir_end);
        }
        if (cir_start != cir_end) {
            flip(supps, pos_cir_start, pos_cir_end);
            alg.compute(cir_start, cir_end, cir_start, cir_end);
            flip(supps, pos_cir_start, pos_cir_end);
        }

        alg.transfer(rays, supps);
        tree.clear();

        // Update the supp vectors for the next_col.
        for (Index i = 0; i < rays.get_number(); ++i) { supps[i].resize(cir_supp_size+2); }
        update_supports(supps, cir_supp_size, pos_ray_start, pos_ray_end);
        update_supports(supps, cir_supp_size, pos_cir_start, pos_cir_end);
        update_supports(supps, cir_supp_size+1, neg_ray_start, neg_ray_end);
        update_supports(supps, cir_supp_size+1, neg_cir_start, neg_cir_end);
        ray_mask.resize(cir_supp_size+2);
        cir_supp_size+=2;

        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << rays.get_number();
        *out  << ",  Time: " << t << "                \n";

        rem.unset(next_col);
        ++cons_added;
    }
}

template <class T, class IndexSet>
SupportRayAlgorithm<T,IndexSet>::SupportRayAlgorithm(
                const ConeT<T>& _cone, const VectorArrayT<T>& _rays, const std::vector<IndexSet>& _supps, 
                const IndexSet& _rel, const Index& _cons_added, const Index& _next, const SUPPORTTREE<IndexSet>& _tree)
        : cone(_cone), rays(_rays), supps(_supps), rel(_rel),
          cons_added(_cons_added), next(_next), tree(_tree), new_rays(0,_rays.get_size()), temp(_cone.num_vars())
{
    _r1_start = _r1_end = _r2_start = _r2_index = _r2_end = 0;
}

template <class T, class IndexSet>
SupportRayAlgorithm<T,IndexSet>::~SupportRayAlgorithm()
{
    new_rays.clear();
    new_supps.clear();
}

template <class T, class IndexSet>
inline
void
SupportRayAlgorithm<T,IndexSet>::compute()
{
    compute(_r1_start, _r1_end, _r2_start, _r2_index, _r2_end);
}

template <class T, class IndexSet>
inline
void
SupportRayAlgorithm<T,IndexSet>::threaded_compute(
                Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end)
{
    _r1_start = r1_start; _r1_end = r1_end; _r2_start = r2_start; _r2_index = r2_index; _r2_end = r2_end;
    ThreadedAlgorithm::threaded_compute();
}

#if 0
template <class T, class IndexSet>
void
SupportRayAlgorithm<T,IndexSet>::compute(
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
            *out << "\r" << buffer;
            *out << "  Size = " << std::setw(8) << rays.get_number() << ", ";
            *out << "  Index = " << std::setw(8) << r1-r1_start << "/" << r1_end-r1_start;
            *out << std::flush;
        }
        ++index_count;

        r1_supp = local_supps[r1];
        r1_count = r1_supp.count();
        cone.get_slack(rays[r1], next, s1);
        if (r1_count == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                temp_supp.set_difference(local_supps[r2], r1_supp);
                if (temp_supp.singleton()) { create_ray(r1, s1, r2);  }
            }
            continue;
        }

        for (Index r2 = r2_start; r2 < r2_index; ++r2) {
            temp_supp.set_difference(r1_supp, local_supps[r2]);
            if (temp_supp.singleton()) { create_ray(r1, s1, r2);  }
        }

        for (Index r2 = r2_index; r2 < r2_end; ++r2) { // Inner loop.
            // Quick sufficient check whether the two rays are adjacent.
            temp_union.set_union(local_supps[r2], r1_supp);
            if (temp_union.count() > cons_added+2) { continue; }
            // Check whether the two rays r1 and r2 are adjacent.
            temp_supp.set_difference(r1_supp, local_supps[r2]);
            if (temp_supp.singleton()) {
                create_ray(r1, s1, r2); 
                continue;
            }
            temp_supp.set_difference(local_supps[r2], r1_supp);
            if (temp_supp.singleton()) {
                create_ray(r1, s1, r2); 
                continue;
            }
            STATS_4ti2(++num_checks;)
            if (tree.dominated(temp_union, r1, r2)) { STATS_4ti2(++num_dominated;) continue; }
            create_ray(r1, s1, r2); 
        }
    }
}

#else

template <class T, class IndexSet>
void
SupportRayAlgorithm<T,IndexSet>::compute(
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
    Size r1_count;
    T s1;

    std::vector<int> indices;

    char buffer[256];
    sprintf(buffer, "  Left = %3d  Col = %3d", rel.count(), next);

    //PlusTree<IndexSet> tmp_tree;
    //tmp_tree.insert(supps);

    //tmp_tree.dump();
    //tree.dump();

    FullTree<IndexSet> neg_tree;
    neg_tree.insert(local_supps, r2_index, r2_end);
    std::vector<Index> neg_indices;
    std::vector<Index> neg_indices0;

    Index index_count = 0;
    // Check negative and positive combinations for adjacency.
    for (Index r1 = r1_start; r1 < r1_end; ++r1) { // Outer loop.
        // Output statistics.
        if (index_count % Globals::output_freq == 0) {
            *out << "\r" << buffer;
            *out << "  Size = " << std::setw(8) << rays.get_number() << ", ";
            *out << "  Index = " << std::setw(8) << r1-r1_start << "/" << r1_end-r1_start;
            *out << std::flush;
        }
        ++index_count;

        r1_supp = local_supps[r1];
        r1_count = r1_supp.count();
        cone.get_slack(rays[r1], next, s1);
        if (r1_count == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (local_supps[r2].singleton_diff(r1_supp)) { create_ray(r1, s1, r2);  }
            }
            continue;
        }

        for (Index r2 = r2_start; r2 < r2_index; ++r2) {
            if (r1_supp.singleton_diff(local_supps[r2])) { create_ray(r1, s1, r2);  }
        }

        zeros.zero();
        indices.clear();
        temp_supp.set_complement(r1_supp);

        // Find the rays whose support differs by one from the current ray's support.
        // These rays are adjacent to r1.
        tree.find_singleton_diff(indices, r1_supp);
        for (unsigned int i = 0; i < indices.size(); ++i) {
            Index r2 = indices[i];
            zeros.set_union(supps[r2]);
            if (r2 >= r2_index && r2 < r2_end) { create_ray(r1, s1, r2); }
        }
        zeros.set_difference(r1_supp);

#if 0
        for (Index r2 = r2_index; r2 < r2_end; ++r2) { // Inner loop.
            if (zeros.set_disjoint(local_supps[r2])) {
                if (r1_supp.count_union(local_supps[r2]) <= cons_added+2) {
                    if (r1_supp.singleton_diff(local_supps[r2])) { create_ray(r1, s1, r2);  continue; }
                    if (local_supps[r2].count_lte_diff(2, r1_supp)) { create_ray(r1, s1, r2); continue; }
                    STATS_4ti2(++num_checks;)
                    temp_union.set_union(local_supps[r2], r1_supp);
                    if (!tree.dominated(temp_union, r1, r2)) { create_ray(r1, s1, r2); }
                }
            }
        }
#endif

#if 1
        neg_indices.clear();
        neg_tree.find(neg_indices, zeros, r1_supp, cons_added+2);
        for (Index i = 0; i < (Index) neg_indices.size(); ++i) { // Inner loop.
            Index r2 = neg_indices[i];
            if (r1_supp.singleton_diff(local_supps[r2])) { create_ray(r1, s1, r2);  continue; }
            if (local_supps[r2].count_lte_diff(2, r1_supp)) { create_ray(r1, s1, r2); continue; }
            STATS_4ti2(++num_checks;)
            temp_union.set_union(local_supps[r2], r1_supp);
            if (!tree.dominated(temp_union, r1, r2)) { create_ray(r1, s1, r2); }
        }
#endif
    }
}
#endif

template <class T, class IndexSet>
void
SupportRayAlgorithm<T,IndexSet>::compute1(
                Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }
    const std::vector<IndexSet>& local_supps = supps;
    Size supps_size = supps[r1_start].get_size();

    // Temporary variables.
    IndexSet temp_supp(supps_size);
    IndexSet temp_union(supps_size);
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
            *out << "\r" << buffer;
            *out << "  Size = " << std::setw(8) << rays.get_number() << ", ";
            *out << "  Index = " << std::setw(8) << r1-r1_start << "/" << r1_end-r1_start;
            *out << std::flush;
        }
        ++index_count;

        r1_supp = local_supps[r1];
        r1_count = r1_supp.count();
        cone.get_slack(rays[r1], next, s1);
        if (r1_count == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                temp_supp.set_difference(local_supps[r2], r1_supp);
                if (temp_supp.singleton()) { create_ray(r1, s1, r2);  }
            }
            continue;
        }

        for (Index r2 = r2_start; r2 < r2_index; ++r2) {
            temp_supp.set_difference(r1_supp, local_supps[r2]);
            if (temp_supp.singleton()) { create_ray(r1, s1, r2);  }
        }

        for (Index r2 = r2_index; r2 < r2_end; ++r2) { // Inner loop.
            // Quick sufficient check whether the two rays are adjacent.
            temp_union.set_union(local_supps[r2], r1_supp);
            if (temp_union.count() > cons_added+2) { continue; }
            // Check whether the two rays r1 and r2 are adjacent.
            temp_supp.set_difference(r1_supp, local_supps[r2]);
            if (temp_supp.singleton()) {
                create_ray(r1, s1, r2); 
                continue;
            }
            temp_supp.set_difference(local_supps[r2], r1_supp);
            if (temp_supp.singleton()) {
                create_ray(r1, s1, r2); 
                continue;
            }
            STATS_4ti2(++num_checks;)
            bool is_adjacent = true;
            for (Index r3 = 0; r3 < rays.get_number(); ++r3) {
                if (IndexSet::set_subset(local_supps[r3], temp_union) && r3 != r1 && r3 != r2) {
                    is_adjacent = false;
                    break;
                }
            }
            if (is_adjacent) { create_ray(r1, s1, r2);  }
        }
    }
}

template <class T, class IndexSet>
inline
void
SupportRayAlgorithm<T,IndexSet>::transfer(VectorArrayT<T>& rays, std::vector<IndexSet>& supps)
{
    rays.transfer(new_rays, 0, new_rays.get_number(), rays.get_number());
    supps.insert(supps.end(), new_supps.begin(), new_supps.end());
    new_supps.clear();
}

template <class T, class IndexSet>
inline
void
SupportRayAlgorithm<T,IndexSet>::create_ray(Index r1, const T& s1, Index r2)
{
    //T s1; cone.get_slack(rays[r1], next, s1);
    T s2; cone.get_slack(rays[r2], next, s2);
    if (r1 < r2) {
        temp.add(rays[r2], s1, rays[r1], -s2);
    }
    else {
        temp.add(rays[r1], s2, rays[r2], -s1);
    }
    temp.normalise();
    new_rays.insert(temp);
    IndexSet temp_supp(supps[r1]);
    temp_supp.set_union(supps[r2]);
    new_supps.push_back(temp_supp);

    DEBUG_4ti2(
    *out << "\nADDING VECTOR.\n";
    *out << "R1: " << r1 << "\n";
    *out << rays[r1] << "\n";
    *out << "R2: " << r2 << "\n";
    *out << rays[r2] << "\n";
    *out << "NEW:\n";
    *out << temp << "\n";
    )
}

template <class T, class IndexSet>
SupportCirAlgorithm<T,IndexSet>::SupportCirAlgorithm(
                const ConeT<T>& _cone, const VectorArrayT<T>& _rays, const std::vector<IndexSet>& _supps, 
                const IndexSet& _rel, const Index& _cons_added, const Index& _next, const SUPPORTTREE<IndexSet>& _tree, const IndexSet& _ray_mask)
        : SupportRayAlgorithm<T,IndexSet>(_cone, _rays, _supps, _rel, _cons_added, _next, _tree), ray_mask(_ray_mask)
{
}

template <class T, class IndexSet>
SupportCirAlgorithm<T,IndexSet>::~SupportCirAlgorithm()
{
}

template <class T, class IndexSet>
inline
void
SupportCirAlgorithm<T,IndexSet>::compute()
{
    compute(SupportRayAlgorithm<T,IndexSet>::_r1_start, SupportRayAlgorithm<T,IndexSet>::_r1_end, 
            SupportRayAlgorithm<T,IndexSet>::_r2_start, SupportRayAlgorithm<T,IndexSet>::_r2_end);
}

template <class T, class IndexSet>
inline
void
SupportCirAlgorithm<T,IndexSet>::threaded_compute(Index r1_start, Index r1_end, Index r2_start, Index r2_end)
{
    SupportRayAlgorithm<T,IndexSet>::_r1_start = r1_start; 
    SupportRayAlgorithm<T,IndexSet>::_r1_end = r1_end; 
    SupportRayAlgorithm<T,IndexSet>::_r2_start = r2_start; 
    SupportRayAlgorithm<T,IndexSet>::_r2_end = r2_end;
    ThreadedAlgorithm::threaded_compute();
}

template <class T, class IndexSet>
void
SupportCirAlgorithm<T,IndexSet>::compute(Index r1_start, Index r1_end, Index r2_start, Index r2_end)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }
    // Copy class variables onto the stack.
    const ConeT<T>& cone = SupportRayAlgorithm<T,IndexSet>::cone;
    const VectorArrayT<T>& rays = SupportRayAlgorithm<T,IndexSet>::rays;
    const std::vector<IndexSet>& supps = SupportRayAlgorithm<T,IndexSet>::supps;
    const Index cons_added = SupportRayAlgorithm<T,IndexSet>::cons_added;
    const Index next = SupportRayAlgorithm<T,IndexSet>::next;
    const SUPPORTTREE<IndexSet>& tree = SupportRayAlgorithm<T,IndexSet>::tree;

    char buffer[256];
    //sprintf(buffer, "  Left = %3d  Col = %3d", rel.count(), next);

    DEBUG_4ti2(*out << "\nComputing circuits for ranges ";)
    DEBUG_4ti2(*out << "R1 [" << r1_start << "..." << r1_end << "] and ";)
    DEBUG_4ti2(*out << "R2 [" << r2_start << "..." << r2_end << "]\n";)

    Index cir_supp_size = ray_mask.get_size();
    IndexSet temp_supp(cir_supp_size);
    IndexSet cir_mask(ray_mask);
    cir_mask.set_complement();

    IndexSet r1_supp(cir_supp_size);
    IndexSet r1_neg_supp(cir_supp_size);
    Size r1_count;
    T s1;

    int index_count = 0;
    for (int r1 = r1_start; r1 < r1_end; ++r1) {
        r1_supp = supps[r1];
        r1_count = r1_supp.count();
        r1_neg_supp.set_intersection(r1_supp, cir_mask);
        r1_neg_supp.swap_odd_n_even();
        cone.get_slack(rays[r1], next, s1);

        if (r2_start == r1) { ++r2_start; IndexSet::swap(r1_supp, r1_neg_supp); }

        if (r1_count == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (r1_neg_supp.set_disjoint(supps[r2])
                    && supps[r2].singleton_diff(r1_supp)) {
                    create_circuit(r1, s1, r2);
                }
            }
            continue;
        }

        for (Index r2 = r2_start; r2 < r2_end; ++r2) {
            if (r1_neg_supp.set_disjoint(supps[r2])
                && r1_supp.count_union(supps[r2]) <= cons_added+2) {
                temp_supp.set_union(r1_supp, supps[r2]);
                if (!tree.dominated(temp_supp, r1, r2)) { create_circuit(r1, s1, r2); }
            }
        }

        if (index_count % Globals::output_freq == 0) {
           *out << "\r" << buffer;
           *out << "  Size = " << std::setw(8) << rays.get_number();
           *out << ",  Index = " << r1 << "/" << r2_end << std::flush;
            DEBUG_4ti2(*out << std::endl;)
        }
        ++index_count;
    }
    *out << "\r" << buffer;
    *out << "  Size = " << std::setw(8) << rays.get_number();
    *out << ",  Index = " << r1_end << "/" << r2_end << std::flush;
    DEBUG_4ti2(*out << std::endl;)
}

template <class T, class IndexSet>
void
SupportCirAlgorithm<T,IndexSet>::create_circuit(Index i1, const T& s1, Index i2)
{
    DEBUG_4ti2(*out << "Creating new circuit.\n";)
    const VectorR<T>& r1 = SupportRayAlgorithm<T,IndexSet>::rays[i1];
    const VectorR<T>& r2 = SupportRayAlgorithm<T,IndexSet>::rays[i2];
    const Index next = SupportRayAlgorithm<T,IndexSet>::next;
    VectorT<T>& temp = SupportRayAlgorithm<T,IndexSet>::temp;

    //T s1; cone.get_slack(r1, next, s1); 
    T s2; SupportRayAlgorithm<T,IndexSet>::cone.get_slack(r2, next, s2); 

    if (s2 > 0) { temp.add(r1,s2,r2,-s1); }
    else { temp.add(r1,-s2,r2,s1); }
    temp.normalise();
    SupportRayAlgorithm<T,IndexSet>::new_rays.insert(temp);

    IndexSet r1_supp(SupportRayAlgorithm<T,IndexSet>::supps[i1]);
    IndexSet r2_supp(SupportRayAlgorithm<T,IndexSet>::supps[i2]);
    IndexSet tmp_union(r1_supp.get_size());
    tmp_union.set_union(r1_supp, r2_supp);
    if (ray_mask.set_disjoint(tmp_union)) {
        if (s1 > 0) { r1_supp.swap_odd_n_even(); }
        else { r2_supp.swap_odd_n_even(); }
    }
    tmp_union.set_union(r1_supp, r2_supp);
    SupportRayAlgorithm<T,IndexSet>::new_supps.push_back(tmp_union);

    DEBUG_4ti2(*out << "Ray" << i1 << " " << r1 << "\t";)
    DEBUG_4ti2(*out << "Sup " << supps[i1] << "\t";)
    DEBUG_4ti2(*out << "Sla " << s1 << "\n";)
    DEBUG_4ti2(*out << "Ray" << i2 << " " << r2 << "\t";)
    DEBUG_4ti2(*out << "Sup " << supps[i2] << "\t";)
    DEBUG_4ti2(*out << "Sla " << s2 << "\n";)
    DEBUG_4ti2(*out << "Ray " << temp << "\t";)
    DEBUG_4ti2(*out << "Sup " << temp_supp << "\n";)
}

#undef DEBUG_4ti2
