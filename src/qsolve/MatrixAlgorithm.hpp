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
#include <iomanip>
#include <cmath>
#include <pthread.h>

#define STATS_4ti2(X) //X
#undef DEBUG_4ti2
#define DEBUG_4ti2(X) //X

using namespace _4ti2_;

template <class T>
MatrixAlgorithm<T>::MatrixAlgorithm()
    : QSolveAlgorithm<T>()
{
}

template <class T>
MatrixAlgorithm<T>::MatrixAlgorithm(QSolveConsOrder o)
    : QSolveAlgorithm<T>(o)
{
}

template <class T>
MatrixAlgorithm<T>::~MatrixAlgorithm()
{
}

template <class T>
void
MatrixAlgorithm<T>::compute(const ConeT<T>& cone, 
            VectorArrayT<T>& rays, std::vector<Index>& ray_ineqs, 
            VectorArrayT<T>& cirs, std::vector<Index>& cir_ineqs)
{
    // The number of constraints added so far.
    Index cons_added = 0;
    // TODO: set size better. Must incorporate circuit constraints.
    if (cone.num_vars()+cone.num_cons() <= IndexSetDS::max_size) {
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


#if 1
// Computes extreme rays of the cone Ax>=0, x>=0.
template <class T> template <class IndexSet>
void
MatrixAlgorithm<T>::compute(
            const ConeT<T>& cone,
            VectorArrayT<T>& rays,
            std::vector<IndexSet>& supps,
            Index& cons_added,
            std::vector<int>& ineqs)
{
    Timer t;
    *out << "Ray Matrix Algorithm.\n";
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

    // Construct the initial set supports.
    for (Index i = 0; i != dim; ++i) { 
        IndexSet supp(n+m, 0);
        supp.set(ineqs[i]);
        supps.push_back(supp);
        rem.unset(ineqs[i]);
    }
    DEBUG_4ti2(*out << "Initial Rays:\n" << rays << "\n";)
    DEBUG_4ti2(*out << "Initial Supps:\n" << supps << "\n";)
    DEBUG_4ti2(*out << "Initial Rem:\n" << rem << "\n";)

    // The total set of relaxed constraints.
    IndexSet rel(rem);
    cone.constraint_set(_4ti2_FR, rel);
    cone.constraint_set(_4ti2_DB, rel);

    Index next = -1;
    IndexSet ray_mask(n+m);
    IndexRanges index_ranges;

    // Construct main algorithm object.
    RayState<T,IndexSet> state(cone, rays, supps, rem, ray_mask, next);
    MatrixRayAlgorithm<IndexSet> alg(state, supps, rel, cons_added, next, index_ranges);
    // Construct threaded algorithm objects.
    std::vector<MatrixRayAlgorithm<IndexSet>*> algs;
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
        state.update(next, pos_start, pos_end);
        // Delete all the vectors with a negative entry in the column next.
        state.remove(neg_start, neg_end);

        // Add new rays and supps.
        alg.transfer();
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->transfer(); }

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
#endif

template <class T> template <class IndexSet>
void
MatrixAlgorithm<T>::compute(
        const ConeT<T>& cone,
        VectorArrayT<T>& rays,
        std::vector<IndexSet>& supps,
        Index& cons_added,
        VectorArrayT<T>& cirs,
        std::vector<Index>& dbls)
{
#if 0
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

    *out << "Circuit Matrix Algorithm.\n";

    // The number of circuit components.
    Size dim_cirs = cirs.get_number();

    // Construct the circuit supports.
    rays.transfer(cirs, 0, cirs.get_number(), rays.get_number());
    std::vector<IndexSet> cir_supps;
    for (Index i = 0; i < dim_cirs; ++i) {
        IndexSet supp(n+m, false);
        supp.set(dbls[i]);
        supps.push_back(supp);
        IndexSet pos_supp(dim_cirs*2, false);
        pos_supp.set(i*2);
        cir_supps.push_back(pos_supp);
    }

    // The mask for the circuit constraints.
    IndexSet cir_mask(n+m, false);
    cone.constraint_set(_4ti2_DB, cir_mask);
    IndexSet ray_mask(cir_mask);
    ray_mask.set_complement();

    // Construct main algorithm object.
    Index next = -1;
    IndexRanges index_ranges;
    MatrixCirAlgorithm<IndexSet> alg(*(new RayState<T,IndexSet>(cone, rays, supps, rem, next)), supps, cir_supps, rem, cons_added, next, index_ranges);
    // Construct threaded algorithm objects.
    std::vector<MatrixCirAlgorithm<IndexSet>*> algs;
    for (Index i = 0; i < Globals::num_threads-1; ++i) { algs.push_back(alg.clone()); }

    while (rays.get_number() > 0 && !rem.empty()) {
        // Find the next column and sort the rays.
        int ray_start, ray_end, cir_start, cir_end;
        int pos_ray_start, pos_ray_end, neg_ray_start, neg_ray_end;
        int pos_cir_start, pos_cir_end, neg_cir_start, neg_cir_end;
        next = next_constraint(cone, rays, supps, cir_supps, rem, ray_mask,
                    ray_start, ray_end, cir_start, cir_end,
                    pos_ray_start, pos_ray_end, neg_ray_start, neg_ray_end,
                    pos_cir_start, pos_cir_end, neg_cir_start, neg_cir_end);

        DEBUG_4ti2(check(cone, rem, rays, supps, cir_supps);)

        // Ouput statistics.
        char buffer[256];
        sprintf(buffer, "%cLeft %3d  Col %3d", ENDL, rem.count(), next);
        *out << buffer << std::flush;

        // Switch the negative circuit supports, so that it is as if all
        // vectors have a positive entry in the next column.
        flip(cir_supps, neg_ray_start, neg_ray_end);
        flip(cir_supps, neg_cir_start, neg_cir_end);
#if 0
        alg.compute_cirs(pos_ray_start, pos_ray_end, neg_ray_start, cir_end);
        alg.compute_cirs(neg_ray_start, neg_ray_end, cir_start, cir_end);
        alg.compute_cirs(cir_start, cir_end, cir_start, cir_end);
#endif
#if 1
        if (pos_ray_start != pos_ray_end) {
            index_ranges.init(pos_ray_start, pos_ray_end, neg_ray_start, neg_ray_start, cir_end); 
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->threaded_compute(); }
            alg.compute();
            // Wait for threads to finish.
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }
        }

        if (neg_ray_start != neg_ray_end) {
            index_ranges.init(neg_ray_start, neg_ray_end, cir_start, cir_start, cir_end); 
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->threaded_compute(); }
            alg.compute();
            // Wait for threads to finish.
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }
        }

        if (cir_start != cir_end) {
            index_ranges.init(cir_start, cir_end, cir_start, cir_start, cir_end); 
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->threaded_compute(); }
            alg.compute();
            // Wait for threads to finish.
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }
        }
#endif

        // Switch back the negative circuit supports.
        flip(cir_supps, neg_ray_start, neg_ray_end);
        flip(cir_supps, neg_cir_start, neg_cir_end);

        // Add new rays and supps.
        alg.transfer();
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->transfer(); }

        // Update the supp vectors for the next_con.
        update_supports(supps, next, ray_start, cir_end);
        resize_supports(cir_supps, 2*dim_cirs+2);
        update_supports(cir_supps, 2*dim_cirs, pos_ray_start, pos_ray_end);
        update_supports(cir_supps, 2*dim_cirs, pos_cir_start, pos_cir_end);
        update_supports(cir_supps, 2*dim_cirs+1, neg_ray_start, neg_ray_end);
        update_supports(cir_supps, 2*dim_cirs+1, neg_cir_start, neg_cir_end);
        dim_cirs+=1;
        cir_mask.set(next);

        //DEBUG_4ti2(check(cone, rem, rays, supps, cir_supps);)

        sprintf(buffer, "%cLeft %3d  Col %3d  Size %8d  Time %8.2fs", ENDL, rem.count(), next, 
                rays.get_number(), t.get_elapsed_time());
        *out << buffer << std::endl;

        rem.unset(next);
        ++cons_added;
    }
#endif
}

// TODO: Check circuit supports.
template <class T> template <class IndexSet>
void
MatrixAlgorithm<T>::check(
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

template <class IndexSet>
MatrixSubAlgorithmBase<IndexSet>::MatrixSubAlgorithmBase(
                RayStateAPI<IndexSet>& _helper,
                std::vector<IndexSet>& _supps, std::vector<IndexSet>& _cir_supps,
                const IndexSet& _rel, const Index& _cons_added, const Index& _next)
        : helper(*_helper.clone()), supps(_supps), cir_supps(_cir_supps), rel(_rel), cons_added(_cons_added), next(_next)
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

    //std::vector<long> zero_count(supps_size, 0);
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

        temp_supp.set_union(r1_supp, rel);
        temp_supp.set_complement();

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
    char buffer[256];

    if (r1_start == r1_end || r2_start == r2_end) { return; }
    DEBUG_4ti2(*out << "\nComputing circuits for ranges ";)
    DEBUG_4ti2(*out << "R1 [" << r1_start << "..." << r1_end << "] and ";)
    DEBUG_4ti2(*out << "R2 [" << r2_start << "..." << r2_end << "]\n";)

    Size s = supps[0].get_size();
    IndexSet temp_supp(s);
    IndexSet temp_zeros(s);
    IndexSet r1_supp(s);
    Size r1_count;

    Size t = cir_supps[0].get_size();
    IndexSet r1_cir_supp(t);

    Index index_count = 0;
    for (Index r1 = r1_start; r1 < r1_end; ++r1) {
        // Output Statistics
        if (index_count % Globals::output_freq == 0) {
            sprintf(buffer, "%cLeft %3d  Col %3d  Index %8d/%-8d", ENDL, rel.count(), next, r1, r2_end);
           *out << buffer << std::flush;
        }
        ++index_count;

        if (r2_start <= r1) { r2_start = r1+1; }

        r1_supp = supps[r1];
        r1_count = r1_supp.count();
        r1_cir_supp = cir_supps[r1];
        helper.set_r1_index(r1);

        if (r1_count == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (supps[r2].singleton_diff(r1_supp)
                    && r1_cir_supp.set_disjoint(cir_supps[r2])) {
                    helper.create_circuit(r2);
                }
            }
            continue;
        }

        temp_supp.set_union(r1_supp, rel);
        temp_supp.set_complement();

        std::vector<Index> con_map;
        helper.project_cone(temp_supp, con_map, temp_zeros);
        DEBUG_4ti2(*out << "NEW ZEROS:\n" << temp_zeros << "\n";)

        for (Index r2 = r2_start; r2 < r2_end; ++r2) {
            if (temp_zeros.singleton_intersection(supps[r2])
                && r1_supp.count_union(supps[r2]) <= cons_added+2
                && r1_cir_supp.set_disjoint(cir_supps[r2])) { 
                //if (supps[r2].singleton_diff(r1_supp)) { helper.create_circuit(r2); continue; }
                //if (r1_supp.singleton_diff(supps[r2])) { helper.create_circuit(r2); continue; }
                temp_supp.set_difference(supps[r2], r1_supp);
                //if (temp_supp.count_lte(2))) { helper.create_circuit(r2); continue; }
                if (helper.is_two_dimensional_face(con_map, temp_supp)) { helper.create_circuit(r2); }
            }
        }
    }
    sprintf(buffer, "%cLeft %3d  Col %3d  Index %8d/%-8d", ENDL, rel.count(), next, r1_end, r2_end);
    *out << buffer << std::flush;
}

template <class IndexSet>
inline
void
MatrixSubAlgorithmBase<IndexSet>::transfer()
{
    //supps.insert(supps.end(), new_supps.begin(), new_supps.end());
    //new_supps.clear();
    //cir_supps.insert(cir_supps.end(), new_cir_supps.begin(), new_cir_supps.end());
    //new_cir_supps.clear();
    helper.transfer();
}

#if 0
template <class IndexSet>
inline  
void    
MatrixSubAlgorithmBase<IndexSet>::set_r1_index(Index r1)
{
    _r1 = r1;
    helper.set_r1_index(r1);
}

template <class IndexSet>
void
MatrixSubAlgorithmBase<IndexSet>::create_ray(Index r2)
{
    IndexSet temp_supp(supps[_r1]);
    temp_supp.set_union(supps[r2]);
    new_supps.push_back(temp_supp);

    helper.create_ray(r2);

    DEBUG_4ti2(*out << "SUPP:\n" << temp_supp << "\n";)
}

template <class IndexSet>
void
MatrixSubAlgorithmBase<IndexSet>::create_rays(Index i1, std::vector<Index>& r2s)
{
    if (r2s.empty()) { return; }
    IndexSet r1_supp(supps[i1]);
    IndexSet r2_supp(r1_supp.get_size());
    for (std::vector<Index>::iterator i2 = r2s.begin(); i2 != r2s.end(); ++i2) {
        r2_supp.set_union(r1_supp, supps[*i2]);
        new_supps.push_back(r2_supp);
    }

    helper.create_rays(i1, r2s);

    DEBUG_4ti2(*out << "SUPP:\n" << temp_supp << "\n";)
}

template <class IndexSet>
void
MatrixSubAlgorithmBase<IndexSet>::create_circuit(Index i2)
{
    DEBUG_4ti2(*out << "Creating new circuit.\n";)
    const Index i1 = _r1;

    bool is_i1_pos = helper.create_circuit(i2);

    IndexSet temp_supp(supps[i1]);
    temp_supp.set_union(supps[i2]);
    new_supps.push_back(temp_supp);

    IndexSet cir_supp(cir_supps[i2]);
    cir_supp.swap_odd_n_even();
    cir_supp.set_union(cir_supps[i1]);
    if (!is_i1_pos) { cir_supp.swap_odd_n_even(); }
    new_cir_supps.push_back(cir_supp);

    DEBUG_4ti2(
        *out << "Ray1 " << i1;
        *out << "Sup1 " << supps[i1] << "\n";
        *out << "Cir1 " << cir_supps[i1] << "\n";
        *out << "Ray2 " << i2;
        *out << "Sup2 " << supps[i2] << "\n";
        *out << "Cir2 " << cir_supps[i2] << "\n";
        *out << "Sup0 " << temp_supp << "\n";
        *out << "Cir0 " << cir_supp << "\n";
    )
}
#endif


template <class IndexSet>
MatrixRayAlgorithm<IndexSet>::MatrixRayAlgorithm(
                RayStateAPI<IndexSet>& helper,
                std::vector<IndexSet>& supps,
                const IndexSet& rel, const Index& cons_added, const Index& next, IndexRanges& _indices)
        : MatrixSubAlgorithmBase<IndexSet>(helper, supps, supps, rel, cons_added, next), indices(_indices)
{
}

template <class IndexSet>
MatrixRayAlgorithm<IndexSet>*
MatrixRayAlgorithm<IndexSet>::clone()
{
    return new MatrixRayAlgorithm(
                MatrixSubAlgorithmBase<IndexSet>::helper,
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
                RayStateAPI<IndexSet>& helper,
                std::vector<IndexSet>& supps, std::vector<IndexSet>& cir_supps, 
                const IndexSet& rel, const Index& cons_added, const Index& next, IndexRanges& _indices)
        : MatrixSubAlgorithmBase<IndexSet>(helper, supps, cir_supps, rel, cons_added, next), indices(_indices)
{
}

template <class IndexSet>
MatrixCirAlgorithm<IndexSet>*
MatrixCirAlgorithm<IndexSet>::clone()
{
    return new MatrixCirAlgorithm(
                MatrixSubAlgorithmBase<IndexSet>::helper,
                MatrixSubAlgorithmBase<IndexSet>::supps,
                MatrixSubAlgorithmBase<IndexSet>::cir_supps, 
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
