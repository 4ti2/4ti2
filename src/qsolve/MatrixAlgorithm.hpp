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
MatrixSubAlgorithmBase<IndexSet>::MatrixSubAlgorithmBase(RayStateAPI<IndexSet>& _state)
        : state(_state), helper(*_state.clone())
{
}

template <class IndexSet>
MatrixSubAlgorithmBase<IndexSet>::~MatrixSubAlgorithmBase()
{
    delete &helper;
}

template <class IndexSet>
void
MatrixSubAlgorithmBase<IndexSet>::compute_rays(
                Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }
    const std::vector<IndexSet>& supps = state.supps;
    Size supps_size = supps[r1_start].get_size();

    IndexSet temp_supp(supps_size);
    IndexSet zeros(supps_size);
    IndexSet r1_supp(supps_size);
    helper.r1_supp.resize(supps_size);
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

    Index index_count = 0;
    // Check negative and positive combinations for adjacency.
    for (Index r1 = r1_start; r1 < r1_end; ++r1) { // Outer loop.
        // Output statistics.
        if (index_count % Globals::output_freq == 0) {
            sprintf(buffer, "%cLeft %3d  Col %3d  Index %8d/%-8d", ENDL, state.rel.count(), state.next, index_count, r1_end-r1_start);
            *out << buffer << std::flush;
        }
        ++index_count;

        r1_supp = supps[r1];
        r1_count = r1_supp.count();
        helper.set_r1_index(r1);
        if (r1_count == state.cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (supps[r2].singleton_diff(r1_supp)) { helper.create_ray(r2); }
            }
            continue;
        }

        for (Index r2 = r2_start; r2 < r2_index; ++r2) {
            if (r1_supp.singleton_diff(supps[r2])) { helper.create_ray(r2);  }
        }
        if (r2_index == r2_end) { continue; }

        temp_supp = r1_supp;
        std::vector<Index> con_map;
        helper.project_cone(temp_supp, con_map, zeros);

        for (Index r2 = r2_index; r2 < r2_end; ++r2) { // Inner loop.
            if (zeros.singleton_intersection(supps[r2])) {
                if (supps[r2].count_lte_2_diff(r1_supp)) { helper.create_ray(r2); continue; }
                if (r1_supp.count_union(supps[r2]) <= state.cons_added+2) {
                    if (r1_supp.singleton_diff(supps[r2])) { helper.create_ray(r2);  continue; }
                    temp_supp.set_difference(supps[r2], r1_supp);
                    //if (temp_supp.count_lte(2)) { create_ray(r2); continue; }
                    if (helper.is_two_dimensional_face(con_map, temp_supp)) { helper.create_ray(r2); }
                }
            }
        }
    }
}


template <class IndexSet>
void
MatrixSubAlgorithmBase<IndexSet>::compute_cirs(Index r1_start, Index r1_end, Index r2_start, Index r2_end)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }
    std::vector<IndexSet>& supps = state.supps;

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
    helper.r1_supp.resize(cir_supp_size);
    Size r1_count;

    int index_count = 0;
    for (int r1 = r1_start; r1 < r1_end; ++r1) {
        //if (r2_start <= r1) { r2_start = r1+1; supps[r1].swap_odd_n_even(); }
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

        temp_supp = r1_supp;
        std::vector<Index> con_map;
        helper.project_cone(temp_supp, con_map, temp_zeros);
        DEBUG_4ti2(*out << "NEW ZEROS:\n" << temp_zeros << "\n";)

        for (Index r2 = r2_start; r2 < r2_end; ++r2) {
            if (temp_zeros.singleton_intersection(supps[r2])
                && r1_neg_supp.set_disjoint(supps[r2])
                && r1_supp.count_union(supps[r2]) <= state.cons_added+2) {
                temp_supp.set_difference(supps[r2], r1_supp);
                //if (temp_supp.count_lte(2))) { helper.create_circuit(r2); continue; }
                if (helper.is_two_dimensional_face(con_map, temp_supp)) { helper.create_circuit(r2); }
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
MatrixSubAlgorithmBase<IndexSet>::transfer()
{
    helper.transfer();
}


template <class IndexSet>
MatrixRayAlgorithm<IndexSet>::MatrixRayAlgorithm(
                RayStateAPI<IndexSet>& state, IndexRanges& _indices)
        : MatrixSubAlgorithmBase<IndexSet>(state), indices(_indices)
{
}

template <class IndexSet>
MatrixRayAlgorithm<IndexSet>*
MatrixRayAlgorithm<IndexSet>::clone()
{
    return new MatrixRayAlgorithm(MatrixSubAlgorithmBase<IndexSet>::state, indices);
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
                RayStateAPI<IndexSet>& state, 
                IndexRanges& _indices)
        : MatrixSubAlgorithmBase<IndexSet>(state), indices(_indices)
{
}

template <class IndexSet>
MatrixCirAlgorithm<IndexSet>*
MatrixCirAlgorithm<IndexSet>::clone()
{
    return new MatrixCirAlgorithm(MatrixSubAlgorithmBase<IndexSet>::state, indices);
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
