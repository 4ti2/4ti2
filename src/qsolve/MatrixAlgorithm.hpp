#include "qsolve/VectorStream.h"
#include "qsolve/VectorArray.h"
#include "qsolve/VectorArrayStream.h"
#include "qsolve/Vector.h"
#include "qsolve/IndexSet.h"
#include "qsolve/IndexSetR.h"
#include "qsolve/IndexSetStream.h"
#include "qsolve/HermiteAlgorithm.h"
#include "qsolve/DiagonalAlgorithm.h"
#include "qsolve/Cone.h"
#include "qsolve/Globals.h"
#include "qsolve/Timer.h"
#include "qsolve/Debug.h"
#include "qsolve/Matrix.h"
#include "qsolve/MultiTree.h"
#include <iomanip>
#include <cmath>
#include <pthread.h>
#include <set>

#define STATS_4ti2(X) //X
#undef DEBUG_4ti2
#define DEBUG_4ti2(X) //X

using namespace _4ti2_;

template <class IndexSet>
MatrixAlgorithm<IndexSet>::MatrixAlgorithm()
{
}

template <class IndexSet>
MatrixAlgorithm<IndexSet>::MatrixAlgorithm(ConsOrder o)
    : order(o)
{
}

template <class IndexSet>
MatrixAlgorithm<IndexSet>::~MatrixAlgorithm()
{
}

#if 1
// Computes extreme rays of the cone Ax>=0, x>=0.
template <class IndexSet>
void
MatrixAlgorithm<IndexSet>::compute_rays(
            const ConeAPI& cone,
            RayStateAPI<IndexSet>& state,
            std::vector<int>& ineqs)
{
    Timer t;
    *out << "Ray Matrix Algorithm.\n";
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

    // Construct the initial set of supports.
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
    MatrixRayAlgorithm<IndexSet> alg(state, supps, rel, cons_added, next, index_ranges);
    // Construct threaded algorithm objects.
    std::vector<MatrixRayAlgorithm<IndexSet>*> algs;
    for (Index i = 0; i < Globals::num_threads-1; ++i) { algs.push_back(alg.clone()); }

    // While there are still rows to choose from.
    while (state.num_gens()>0 && !rem.empty()) {
        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)

        // Choose the next constraint and sort rays.
        Index pos_start, pos_end, neg_start, neg_end;
        next = state.next_constraint(order, rem, pos_start, pos_end, neg_start, neg_end);

        // Check to see whether the next constraint is redundant. If so, ignore it.
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

            // Run threads.
            index_ranges.init(r1_start, r1_end, r2_start, r2_index, r2_end);
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->threaded_compute(); }
            // Run primary algorithm.
            alg.compute();

            // Wait for threads to finish.
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }
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
#endif

#if 0
// Computes extreme rays of the cone Ax>=0, x>=0.
template <class IndexSet>
void
MatrixAlgorithm<IndexSet>::compute_rays(
            const ConeAPI& cone,
            RayStateAPI<IndexSet>& state,
            std::vector<int>& ineqs)
{
    Timer t;

    std::vector<IndexSet>& supps = state.supps;
    Index& cons_added = state.cons_added;
    IndexSet& rel = state.rel;
    IndexSet& rem = state.rem;

    *out << "Ray Matrix Algorithm.\n";
    DEBUG_4ti2(*out << "CONSTRAINT MATRIX:\n" << cone.get_matrix() << "\n";)

    // The number of variables.
    Size n = cone.num_vars();
    // The number of constraints.
    Size m = cone.num_cons();

    // The dimension
    Size dim = state.num_gens();

    // The set of constraints to be processed.
    rem.zero();
    cone.get_constraint_set(_4ti2_LB, rem);

    // Construct the initial set supports.
    for (Index i = 0; i != dim; ++i) { 
        IndexSet supp(n+m, 0);
        supp.set(ineqs[i]);
        supps.push_back(supp);
        rem.unset(ineqs[i]);
    }
    DEBUG_4ti2(*out << "Initial Supps:\n" << supps << "\n";)
    DEBUG_4ti2(*out << "Initial Rem:\n" << rem << "\n";)

    state.supp_types.clear();
    state.supp_types.resize(n+m, _4ti2_LB);
    state.supps_to_cons.resize(n+m, -1);
    state.cons_to_supps.resize(n+m, -1);
    for (Index i = 0; i < n+m; ++i) {
        state.supps_to_cons[i]=i;
        state.cons_to_supps[i]=i;
    }

    // The total set of relaxed constraints.
    rel = rem;
    cone.get_constraint_set(_4ti2_FR, rel);
    cone.get_constraint_set(_4ti2_DB, rel);

    IndexSet ray_mask(n+m);
    IndexRanges index_ranges;

    // Construct main algorithm object.
    MatrixRayAlgorithm<IndexSet> alg(state, supps, rel, cons_added, state.next, index_ranges);
    // Construct threaded algorithm objects.
    std::vector<MatrixRayAlgorithm<IndexSet>*> algs;
    for (Index i = 0; i < Globals::num_threads-1; ++i) { algs.push_back(alg.clone()); }

    // While there are still rows to choose from.
    while (state.num_gens()>0 && !rem.empty()) {
        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)

        // Choose the next constraint and sort rays.
        Index pos_start, pos_end, neg_start, neg_end;
        state.next = state.next_constraint(order, rem, pos_start, pos_end, neg_start, neg_end);

        // Check to see whether the next constraint is redundant. If so, ignore it.
        if (neg_start == neg_end) { rem.unset(state.next); continue; }

        char buffer[256];
        sprintf(buffer, "%cLeft %3d  Col %3d  Size %8d", ENDL, rem.count(), state.next, state.num_gens());
        *out << buffer << std::flush;

        // Next, we assign index ranges for the inner and outer for loops so
        // that the outer loop is smaller than the inner loop.
        Index r1_start, r1_end, r2_start, r2_end;
        Size pos_count = pos_end-pos_start;
        Size neg_count = neg_end-neg_start;
        if (pos_count <= neg_count) {
            r1_start = pos_start; r1_end = pos_end;
            r2_start = neg_start; r2_end = neg_end;
        }
        else {
            r2_start = pos_start; r2_end = pos_end;
            r1_start = neg_start; r1_end = neg_end;
        }

        // We sort the r2's into vectors where r2_supp.count()==cons_added+1.
        Index r2_index = state.sort_count(cons_added+1, r2_start, r2_end);

        // Run threads.
        index_ranges.init(r1_start, r1_end, r2_start, r2_index, r2_end);
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->threaded_compute(); }
        // Run primary algorithm.
        alg.compute();

        // Wait for threads to finish.
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }

        // Update the support vectors for the next_con.
        state.update(state.next, pos_start, pos_end);
        // Delete all the vectors with a negative entry in the column next.
        state.remove(neg_start, neg_end);

        // Add new rays and supps.
        alg.transfer();
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->transfer(); }

        rem.unset(state.next);
        rel.unset(state.next);
        ++cons_added;

        // Output statistics.
        sprintf(buffer, "%cLeft %3d  Col %3d  Size %8d  Time %8.2fs", ENDL, rem.count(), state.next, 
                state.num_gens(), t.get_elapsed_time());
        *out << buffer << std::endl;

        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)
        DEBUG_4ti2(state.check());
    }

    // Clean up threaded algorithms objects.
    for (Index i = 0; i < Globals::num_threads-1; ++i) { delete algs[i]; }
}
#endif

template <class IndexSet>
void
MatrixAlgorithm<IndexSet>::compute_cirs(
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

    *out << "Circuit Matrix Algorithm.\n";

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
    MatrixCirAlgorithm<IndexSet> alg(state, supps, rem, cons_added, next, index_ranges);
    std::vector<MatrixCirAlgorithm<IndexSet>*> algs;
    for (Index i = 0; i < Globals::num_threads-1; ++i) { algs.push_back(alg.clone()); }

    while (state.num_gens() > 0 && !rem.empty()) {
        // Find the next constraint.
        int pos_ray_start, pos_ray_end, neg_ray_start, neg_ray_end;
        int pos_cir_start, pos_cir_end, neg_cir_start, neg_cir_end;
        next = state.next_constraint(order, rem, ray_mask,
                    pos_ray_start, pos_ray_end, neg_ray_start, neg_ray_end,
                    pos_cir_start, pos_cir_end, neg_cir_start, neg_cir_end);

        //DEBUG_4ti2(print_debug_diagnostics(cone, rays, supps, next);)

        // Ouput statistics.
        char buffer[256];
        sprintf(buffer, "%cLeft %3d  Col %3d", ENDL, rem.count(), next);
        *out << buffer << std::flush;
        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)

        // TODO: Make threaded. // Change algorithm so that supps are constant.
        state.flip(pos_cir_start, pos_cir_end);
        alg.compute_cirs(pos_ray_start, neg_cir_end, pos_cir_start, neg_ray_end);
        state.flip(neg_cir_start, neg_cir_end);

        alg.transfer();
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->transfer(); }

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
MatrixSubAlgorithmBase<IndexSet>::MatrixSubAlgorithmBase(
                RayStateAPI<IndexSet>& _state, std::vector<IndexSet>& _supps, 
                const IndexSet& _rel, const Index& _cons_added, const Index& _next)
        : state(_state), helper(*_state.clone()), supps(_supps), rel(_rel), cons_added(_cons_added), next(_next)
{
}

template <class IndexSet>
MatrixSubAlgorithmBase<IndexSet>::~MatrixSubAlgorithmBase()
{
    delete &helper;
}

#if 1
template <class IndexSet>
void
MatrixSubAlgorithmBase<IndexSet>::compute_rays(
                Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }
    const std::vector<IndexSet>& local_supps = supps;
    Size supps_size = supps[r1_start].get_size();

    IndexSet temp_supp(supps_size);
    IndexSet temp_diff(supps_size);
    IndexSet temp_union(supps_size);
    IndexSet zeros(supps_size);
    IndexSet r1_supp(supps_size);
    Size r1_count;

    char buffer[256];

#if 0
    // Look for redundant constraints.
    temp_supp.zero();
    for (Index i = 0; i < (Index) supps.size(); ++i) {
        temp_supp.set_intersection(supps[i]);
    }
    if (!temp_supp.empty()) { *out << "All Intersection:\n" << temp_supp << "\n"; }
#endif
    
    std::vector<Index> supps_to_cons(supps_size);
    for (Index i = 0; i < supps_size; ++i) { supps_to_cons[i] = i; }
    std::vector<Index> cons_to_supps = supps_to_cons;

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
        r1_count = r1_supp.count();
        helper.set_r1_index(r1);
        if (r1_count == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (local_supps[r2].singleton_diff(r1_supp)) { helper.create_ray(r2); }
            }
            continue;
        }

        for (Index r2 = r2_start; r2 < r2_index; ++r2) {
            if (r1_supp.singleton_diff(local_supps[r2])) { helper.create_ray(r2);  }
        }
        if (r2_index == r2_end) { continue; }

        temp_supp = r1_supp;
        std::vector<Index> con_map;
        helper.project_cone(temp_supp, con_map, zeros);

        for (Index r2 = r2_index; r2 < r2_end; ++r2) { // Inner loop.
            if (zeros.singleton_intersection(local_supps[r2])) {
                if (local_supps[r2].count_lte_2_diff(r1_supp)) { helper.create_ray(r2); continue; }
                if (r1_supp.count_union(local_supps[r2]) <= cons_added+2) {
                    if (r1_supp.singleton_diff(local_supps[r2])) { helper.create_ray(r2);  continue; }
                    temp_supp.set_difference(local_supps[r2], r1_supp);
                    //if (temp_supp.count_lte(2)) { create_ray(r2); continue; }
                    if (helper.is_two_dimensional_face(con_map, temp_supp)) { helper.create_ray(r2); }
                }
            }
        }

#if 0
        for (Index r2 = r2_index; r2 < r2_end; ++r2) { // Inner loop.
            if (zeros.set_disjoint(local_supps[r2])) {
                temp_union.set_union(local_supps[r2], r1_supp);
                if (temp_union.count() <= cons_added+2) {
                    // Check whether the two rays r1 and r2 are adjacent.
                    temp_diff.set_difference(r1_supp, local_supps[r2]);
                    if (temp_diff.singleton()) { create_ray(r2); continue; }
                    temp_diff.set_difference(local_supps[r2], r1_supp);
                    if (temp_diff.count() == 2) { create_ray(r2); continue; }
                    if (is_two_dimensional_face(con_map, temp_diff)) {
                        create_ray(r2); 
                    }
                }
            }
            else {
                temp_diff.set_difference(local_supps[r2], r1_supp);
                if (temp_diff.singleton()) { create_ray(r2); }
            }
        }
#endif

    }
}
#endif

template <class IndexSet>
void
MatrixSubAlgorithmBase<IndexSet>::compute_cirs(Index r1_start, Index r1_end, Index r2_start, Index r2_end)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }
    // TODO: Copy class variables onto the stack???

    char buffer[256];

    DEBUG_4ti2(*out << "\nComputing circuits for ranges ";)
    DEBUG_4ti2(*out << "R1 [" << r1_start << ",...," << r1_end << "] and ";)
    DEBUG_4ti2(*out << "R2 [" << r2_start << ",...," << r2_end << "]\n";)

    Index cir_supp_size = state.ray_mask.get_size();
    IndexSet temp_supp(cir_supp_size);
    IndexSet temp_zeros(cir_supp_size);
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

        if (r1_count == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (r1_neg_supp.set_disjoint(supps[r2])
                    && supps[r2].singleton_diff(r1_supp)) {
                    helper.create_circuit(r2);
                }
            }
            continue;
        }

        temp_supp = r1_supp;
        std::vector<Index> con_map;
        helper.project_cone(temp_supp, con_map, temp_zeros);
        DEBUG_4ti2(*out << "NEW ZEROS:\n" << temp_zeros << "\n";)

        for (Index r2 = r2_start; r2 < r2_end; ++r2) {
            if (temp_zeros.singleton_intersection(supps[r2])
                && r1_neg_supp.set_disjoint(supps[r2])
                && r1_supp.count_union(supps[r2]) <= cons_added+2) {
                temp_supp.set_difference(supps[r2], r1_supp);
                //if (temp_supp.count_lte(2))) { helper.create_circuit(r2); continue; }
                if (helper.is_two_dimensional_face(con_map, temp_supp)) { helper.create_circuit(r2); }
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
MatrixSubAlgorithmBase<IndexSet>::transfer()
{
    helper.transfer();
}


template <class IndexSet>
MatrixRayAlgorithm<IndexSet>::MatrixRayAlgorithm(
                RayStateAPI<IndexSet>& state,
                std::vector<IndexSet>& supps,
                const IndexSet& rel, const Index& cons_added, const Index& next, IndexRanges& _indices)
        : MatrixSubAlgorithmBase<IndexSet>(state, supps, rel, cons_added, next), indices(_indices)
{
}

template <class IndexSet>
MatrixRayAlgorithm<IndexSet>*
MatrixRayAlgorithm<IndexSet>::clone()
{
    return new MatrixRayAlgorithm(
                MatrixSubAlgorithmBase<IndexSet>::state,
                MatrixSubAlgorithmBase<IndexSet>::supps, 
                MatrixSubAlgorithmBase<IndexSet>::rel,
                MatrixSubAlgorithmBase<IndexSet>::cons_added,
                MatrixSubAlgorithmBase<IndexSet>::next,
                indices);
}

template <class IndexSet>
inline
void
MatrixRayAlgorithm<IndexSet>::compute()
{
    Index r1_start, r1_end, r2_start, r2_index, r2_end;
    indices.next(r1_start, r1_end, r2_start, r2_index, r2_end);
    while (r1_start != r1_end) {
        MatrixSubAlgorithmBase<IndexSet>::compute_rays(r1_start, r1_end, r2_start, r2_index, r2_end);
        indices.next(r1_start, r1_end, r2_start, r2_index, r2_end);
    }
}

template <class IndexSet>
MatrixCirAlgorithm<IndexSet>::MatrixCirAlgorithm(
                RayStateAPI<IndexSet>& state, std::vector<IndexSet>& supps, 
                const IndexSet& rel, const Index& cons_added, const Index& next, IndexRanges& _indices)
        : MatrixSubAlgorithmBase<IndexSet>(state, supps, rel, cons_added, next), indices(_indices)
{
}

template <class IndexSet>
MatrixCirAlgorithm<IndexSet>*
MatrixCirAlgorithm<IndexSet>::clone()
{
    return new MatrixCirAlgorithm(
                MatrixSubAlgorithmBase<IndexSet>::state,
                MatrixSubAlgorithmBase<IndexSet>::supps,
                MatrixSubAlgorithmBase<IndexSet>::rel,
                MatrixSubAlgorithmBase<IndexSet>::cons_added,
                MatrixSubAlgorithmBase<IndexSet>::next,
                indices);
}

template <class IndexSet>
inline
void
MatrixCirAlgorithm<IndexSet>::compute()
{
    Index r1_start, r1_end, r2_start, r2_index, r2_end;
    indices.next(r1_start, r1_end, r2_start, r2_index, r2_end);
    while (r1_start != r1_end) {
        MatrixSubAlgorithmBase<IndexSet>::compute_cirs(r1_start, r1_end, r2_start, r2_end);
        indices.next(r1_start, r1_end, r2_start, r2_index, r2_end);
    }
}

#undef DEBUG_4ti2
#define DEBUG_4ti2(X) //X
