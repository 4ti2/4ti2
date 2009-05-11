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

#define ENDL '\r'

#define STATS_4ti2(X) //X
#undef DEBUG_4ti2
#define DEBUG_4ti2(X) //X

using namespace _4ti2_;

//template <class T, class IndexSet>
//pthread_mutex_t MatrixRayAlgorithm::m = PTHREAD_MUTEX_INITIALIZER;

template <class T, class IndexSet>
MatrixAlgorithm<T,IndexSet>::MatrixAlgorithm()
    : QSolveAlgorithm<T,IndexSet>()
{
}

template <class T, class IndexSet>
MatrixAlgorithm<T,IndexSet>::MatrixAlgorithm(QSolveConsOrder o)
    : QSolveAlgorithm<T,IndexSet>(o)
{
}

template <class T, class IndexSet>
MatrixAlgorithm<T,IndexSet>::~MatrixAlgorithm()
{
}

template <class T, class IndexSet>
void
MatrixAlgorithm<T,IndexSet>::compute(const ConeT<T>& cone, 
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
MatrixAlgorithm<T,IndexSet>::compute(
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

    // Construct main algorithm object.
    MatrixRayAlgorithm<T,IndexSet> alg(cone, rays, supps, rel, cons_added);
    // Construct threaded algorithm objects.
    std::vector<MatrixRayAlgorithm<T,IndexSet>*> algs;
    for (Index i = 0; i < Globals::num_threads-1; ++i) {
            algs.push_back(new MatrixRayAlgorithm<T,IndexSet>(cone, rays, supps, rel, cons_added));
    }

    // While there are still rows to choose from.
    while (rays.get_number()>0 && !rem.empty()) {
        DEBUG_4ti2(*out << "RAYS:\n" << rays << "\n";)
        DEBUG_4ti2(*out << "SUPPORTS:\n" << supps << "\n";)

        // Choose the next constraint and sort rays.
        Index pos_start, pos_end, neg_start, neg_end;
        Index next = next_constraint(cone, rem, rays, supps,
                        pos_start, pos_end, neg_start, neg_end);

        char buffer[256];
        sprintf(buffer, "  Left = %3d  Col = %3d", rem.count(), next);
        *out << ENDL << buffer;
        *out << "  Size = " << std::setw(8) << rays.get_number() << "  Time: " << t;

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
        Index r2_index = r2_start;
        for (Index r2 = r2_start; r2 < r2_end; ++r2) {
            if (supps[r2].count() == cons_added+1) {
                rays.swap_rows(r2, r2_index);
                IndexSet::swap(supps[r2], supps[r2_index]);
                ++r2_index;
            }
        }

        // Run threads.
        Index increment = (r1_end-r1_start)/Globals::num_threads;
        for (Index i = 0; i < Globals::num_threads-1; ++i) {
            algs[i]->threaded_compute(next, r1_start, r1_start+increment, 
                            r2_start, r2_index, r2_end);
            r1_start += increment;
        }
        // Run primary algorithm.
        alg.compute(next, r1_start, r1_end, r2_start, r2_index, r2_end);

        // Wait for threads to finish.
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }

        // Update the support vectors for the next_con.
        update_supports(supps, next, pos_start, pos_end);
        // Delete all the vectors with a negative entry in the column next.
        supps.erase(supps.begin()+neg_start, supps.begin()+neg_end);
        rays.remove(neg_start, neg_end);

        // Add new rays and supps.
        alg.transfer(rays, supps);
        for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->transfer(rays, supps); }

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
MatrixAlgorithm<T,IndexSet>::compute(
        const ConeT<T>& cone,
        VectorArrayT<T>& rays,
        std::vector<IndexSet>& supps,
        Index& cons_added,
        VectorArrayT<T>& cirs,
        std::vector<Index>& dbls)
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
    MatrixCirAlgorithm<T,IndexSet> alg(cone, rays, supps, cir_supps, rem, cons_added);
    // Construct threaded algorithm objects.
    std::vector<MatrixCirAlgorithm<T,IndexSet>*> algs;
    for (Index i = 0; i < Globals::num_threads-1; ++i) {
        algs.push_back(new MatrixCirAlgorithm<T,IndexSet>(cone, rays, 
                            supps, cir_supps, rem, cons_added));
    }

    while (rays.get_number() > 0 && !rem.empty()) {
        // Find the next column and sort the rays.
        int ray_start, ray_end, cir_start, cir_end;
        int pos_ray_start, pos_ray_end, neg_ray_start, neg_ray_end;
        int pos_cir_start, pos_cir_end, neg_cir_start, neg_cir_end;
        Index next_con = next_constraint(cone, rays, supps, cir_supps, rem, ray_mask,
                    ray_start, ray_end, cir_start, cir_end,
                    pos_ray_start, pos_ray_end, neg_ray_start, neg_ray_end,
                    pos_cir_start, pos_cir_end, neg_cir_start, neg_cir_end);

        DEBUG_4ti2(check(cone, rem, rays, supps, cir_supps);)

        // Ouput statistics.
        char buffer[256];
        sprintf(buffer, "  Left = %3d  Col = %3d", rem.count(), next_con);
        *out << ENDL << buffer;
        *out << "  Size = " << std::setw(8) << rays.get_number() << "  Time: " << t;
        DEBUG_4ti2(*out << std::endl;)

        // Switch the negative circuit supports, so that it is as if all
        // vectors have a positive entry in the next column.
        flip(cir_supps, neg_ray_start, neg_ray_end);
        flip(cir_supps, neg_cir_start, neg_cir_end);
#if 0
        alg.compute(next_con, pos_ray_start, pos_ray_end, neg_ray_start, cir_end);
        alg.compute(next_con, neg_ray_start, neg_ray_end, cir_start, cir_end);
        alg.compute(next_con, cir_start, cir_end, cir_start, cir_end);
#endif
        if (pos_ray_start != pos_ray_end) {
            Index increment = (pos_ray_end-pos_ray_start)/Globals::num_threads;
            Index start = pos_ray_start;
            for (Index i = 0; i < Globals::num_threads-1; ++i) {
                algs[i]->threaded_compute(next_con, start, start+increment, neg_ray_start, cir_end);
                start += increment;
            }
            alg.compute(next_con, start, pos_ray_end, neg_ray_start, cir_end);
            // Wait for threads to finish.
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }
        }

        if (neg_ray_start != neg_ray_end) {
            Index increment = (neg_ray_end-neg_ray_start)/Globals::num_threads;
            Index start = neg_ray_start;
            for (Index i = 0; i < Globals::num_threads-1; ++i) {
                algs[i]->threaded_compute(next_con, start, start+increment, cir_start, cir_end);
                start += increment;
            }
            alg.compute(next_con, start, neg_ray_end, cir_start, cir_end);
            // Wait for threads to finish.
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }
        }

        if (cir_start != cir_end) {
            double frac = 1-sqrt(1-1/((double) Globals::num_threads+1));
            Index start = cir_start;
            for (Index i = 0; i < Globals::num_threads-1; ++i) {
                Index end = (Index) ((cir_end-start)*frac+start);
                algs[i]->threaded_compute(next_con, start, end, start, cir_end);
                start = end;
            }
            alg.compute(next_con, start, cir_end, start, cir_end);
            // Wait for threads to finish.
            for (Index i = 0; i < Globals::num_threads-1; ++i) { algs[i]->wait(); }
        }

        // Switch back the negative circuit supports.
        flip(cir_supps, neg_ray_start, neg_ray_end);
        flip(cir_supps, neg_cir_start, neg_cir_end);

        // Add new rays and supps.
        alg.transfer(rays, supps, cir_supps);
        for (Index i = 0; i < Globals::num_threads-1; ++i) { 
                algs[i]->transfer(rays, supps, cir_supps);
        }

        // Update the supp vectors for the next_con.
        update_supports(supps, next_con, ray_start, cir_end);
        for (Index i = 0; i < rays.get_number(); ++i) { cir_supps[i].resize(2*dim_cirs+2); }
        update_supports(cir_supps, 2*dim_cirs, pos_ray_start, pos_ray_end);
        update_supports(cir_supps, 2*dim_cirs, pos_cir_start, pos_cir_end);
        update_supports(cir_supps, 2*dim_cirs+1, neg_ray_start, neg_ray_end);
        update_supports(cir_supps, 2*dim_cirs+1, neg_cir_start, neg_cir_end);
        dim_cirs+=1;
        cir_mask.set(next_con);

        //DEBUG_4ti2(check(cone, rem, rays, supps, cir_supps);)

        *out << ENDL << buffer;
        *out << "  Size = " << std::setw(8) << rays.get_number();
        *out  << ",  Time: " << t << "                \n";

        rem.unset(next_con);
        ++cons_added;
    }
}

// TODO: Check circuit supports.
template <class T, class IndexSet>
void
MatrixAlgorithm<T,IndexSet>::check(
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

template <class T, class IndexSet>
MatrixRayAlgorithm<T,IndexSet>::MatrixRayAlgorithm(
                const ConeT<T>& _cone, const VectorArrayT<T>& _rays, const std::vector<IndexSet>& _supps, 
                const IndexSet& _rel, const Index& _cons_added)
        : cone(_cone), rays(_rays), supps(_supps),
          rel(_rel), cons_added(_cons_added), new_rays(0,_rays.get_size()), temp(_cone.num_vars())
{
    _next = _r1_start = _r1_end = _r2_start = _r2_index = _r2_end = 0;
}

template <class T, class IndexSet>
MatrixRayAlgorithm<T,IndexSet>::~MatrixRayAlgorithm()
{
    new_rays.clear();
    new_supps.clear();
}

template <class T, class IndexSet>
void
MatrixRayAlgorithm<T,IndexSet>::project_cone(
                const ConeT<T>& cone,
                const IndexSet& zero_supp,
                MatrixT<T>& trans,
                std::vector<Index>& con_map)
{
    Index n = cone.num_vars();
    Index m = cone.num_cons();
    DEBUG_4ti2(*out << "Constraint Matrix:\n" << cone.get_matrix() << "\n";)
    DEBUG_4ti2(*out << "Zero Support:\n" << zero_supp << "\n";)

    IndexSet rows(n, 0);
    IndexSet cols(m, 0);
    con_map.reserve(n+m);
    for (Index i = 0; i < n; ++i) {
        if (zero_supp[i]) { con_map.push_back(i); }
        else { rows.set(i); con_map.push_back(-1); }
    }
    for (Index i = n; i < m+n; ++i) {
        if (zero_supp[i]) { 
            cols.set(i-n);
            con_map.push_back(i);
        }
    }
    DEBUG_4ti2(*out << "Col Map:\n" << con_map << "\n";)
    DEBUG_4ti2(*out << "Column projection:\n" << cols << "\n";)
    DEBUG_4ti2(*out << "Row relaxation:\n" << rows << "\n";)

    trans.init(n, cols.count());
    trans.assign_trans(cone.get_matrix(), cols, IndexSetR(0,n));
    DEBUG_4ti2(*out << "NEW TRANS:\n"; trans.print();)

    typename MatrixT<T>::Pivots pivots;
    trans.row_diagonalise(rows, IndexSetR(0,trans.get_n()), &pivots);
    DEBUG_4ti2(*out << "Triangle:\n"; trans.print();)
    rows.set_complement();
    trans.row_diagonalise(pivots, rows);
    trans.row_normalise();
    DEBUG_4ti2(*out << "Diagonal:\n"; trans.print();)

    for (Index i = 0; i < (Index) pivots.size(); ++i) {
        con_map[pivots[i].first] = con_map[pivots[i].second+n];
        con_map[pivots[i].second+n] = -1;
    } 
    DEBUG_4ti2(*out << "Col Map:\n" << con_map << "\n";)
}

template <class T, class IndexSet>
void
MatrixRayAlgorithm<T,IndexSet>::zero_cols(
                const MatrixT<T>& matrix,
                const std::vector<Index>& con_map,
                IndexSet& zeros)
{
    zeros.zero();
    Index n = matrix.get_m();
    Index m = matrix.get_n();
    for (Index i = 0; i < n; ++i) {
        if (con_map[i] < 0) { continue; }
        zeros.set(con_map[i]);
        for (Index j = 0; j < m; ++j) {
            // TODO: matrix row.
            if (matrix(i,j) != 0 && con_map[j+n] >= 0) { zeros.unset(con_map[i]); break; }
        }
    }
    DEBUG_4ti2(*out << "Zeros:\n" << zeros << "\n";)
}

// Checks whether the given support determines a two dimensional face of the cone.
template <class T, class IndexSet>
bool
MatrixRayAlgorithm<T,IndexSet>::is_two_dimensional_face(
            const MatrixT<T>& trans,
            const std::vector<Index>& con_map,
            const IndexSet& diff)
{
    DEBUG_4ti2(*out << "\nis_two_dimensional_face\n";)
    Index n = trans.get_m();
    Index m = trans.get_n();

    IndexSet cols(n,0);
    for (Index i = 0; i < n; ++i) { 
        if (con_map[i] != -1 && diff[con_map[i]]) { cols.set(i); }
    }
    IndexSet rows(m,0);
    for (Index i = n; i < m+n; ++i) {
        if (con_map[i] != -1 && !diff[con_map[i]]) { rows.set(i-n); }
    }
    DEBUG_4ti2(*out << "Cols:\n" << cols << "\nRows:\n" << rows << "\n";)

    Size proj_n = cols.count();
    Size proj_m = rows.count();
    if (proj_n <= 1) { return true; }
    matrix.init(proj_n, proj_m);
    matrix.assign(trans, cols, rows);
    DEBUG_4ti2(*out << "MATRIX:\n"; matrix.print();)
    Size rank = matrix.row_triangulate();
    DEBUG_4ti2(*out << "Rank: " << rank << "\n";)
    if (rank == matrix.get_m()-1) { return true; }

    return false;
}

template <class T, class IndexSet>
inline
void
MatrixRayAlgorithm<T,IndexSet>::compute()
{
    compute(_next, _r1_start, _r1_end, _r2_start, _r2_index, _r2_end);
}

template <class T, class IndexSet>
inline
void
MatrixRayAlgorithm<T,IndexSet>::threaded_compute(
                Index next, Index r1_start, Index r1_end, 
                Index r2_start, Index r2_index, Index r2_end)
{
    _next = next; _r1_start = r1_start; _r1_end = r1_end; _r2_start = r2_start; _r2_index = r2_index; _r2_end = r2_end;
    ThreadedAlgorithm::threaded_compute();
}

template <class T, class IndexSet>
void
MatrixRayAlgorithm<T,IndexSet>::compute(
                Index next, Index r1_start, Index r1_end, 
                Index r2_start, Index r2_index, Index r2_end)
{
    const std::vector<IndexSet>& local_supps = supps;
    _next = next; //TODO: This is a hack.

    Index n = cone.num_vars();
    Index m = cone.num_cons();

    // Temporary variables.
    MatrixT<T> trans;
    IndexSet temp_supp(n+m);
    IndexSet temp_diff(n+m);
    IndexSet temp_union(n+m);
    IndexSet zeros(n+m);
    IndexSet r1_supp(n+m);
    Size r1_count;
    T s1;

    char buffer[256];
    sprintf(buffer, "  Left = %3d  Col = %3d", rel.count(), next);

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
            *out << ENDL << buffer;
            //*out << "  Size = " << std::setw(8) << rays.get_number()-neg_count << ", ";
            *out << "  Size = " << std::setw(8) << rays.get_number()+new_rays.get_number() << ", ";
            *out << "  Index = " << index_count << "/" << r1_end-r1_start << std::flush;
        }
        ++index_count;

        r1_supp = local_supps[r1];
        r1_count = r1_supp.count();
        cone.get_slack(rays[r1], next, s1);
        if (r1_count == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                //if (local_supps[r2].singleton_diff(r1_supp)) { create_ray(r1, s1 ,r2); }
                temp_diff.set_difference(local_supps[r2], r1_supp);
                if (temp_diff.singleton()) { create_ray(r1, s1, r2); }
            }
            continue;
        }

        for (Index r2 = r2_start; r2 < r2_index; ++r2) {
            //if (r1_supp.singleton_diff(local_supps[r2])) { create_ray(r1, s1, r2);  }
            temp_diff.set_difference(r1_supp, local_supps[r2]);
            if (temp_diff.singleton()) { create_ray(r1, s1, r2); continue; }
        }
        if (r2_index == r2_end) { continue; }

        temp_supp.set_union(r1_supp, rel);
        temp_supp.set_complement();

        std::vector<Index> con_map;
        project_cone(cone, temp_supp, trans, con_map);
        zero_cols(trans, con_map, zeros);

#if 0
        for (Index r2 = r2_index; r2 < r2_end; ++r2) { // Inner loop.
            if (zeros.singleton_intersection(local_supps[r2])) {
                if (local_supps[r2].singleton_diff(r1_supp)) { create_ray(r1, s1, r2); continue; }
                if (r1_supp.count_union(local_supps[r2]) <= cons_added+2) {
                    if (r1_supp.singleton_diff(local_supps[r2])) { create_ray(r1, s1, r2);  continue; }
                    temp_supp.set_difference(local_supps[r2], r1_supp);
                    if (temp_supp.count_lte(2)) { create_ray(r1, s1, r2); continue; }
                    if (is_two_dimensional_face(trans, con_map, temp_supp)) { create_ray(r1, s1, r2); }
                }
            }
        }
#endif

#if 1
        for (Index r2 = r2_index; r2 < r2_end; ++r2) { // Inner loop.
            if (zeros.set_disjoint(local_supps[r2])) {
                temp_union.set_union(local_supps[r2], r1_supp);
                if (temp_union.count() <= cons_added+2) {
                    // Check whether the two rays r1 and r2 are adjacent.
                    temp_diff.set_difference(r1_supp, local_supps[r2]);
                    if (temp_diff.singleton()) { create_ray(r1, s1, r2); continue; }
                    temp_diff.set_difference(local_supps[r2], r1_supp);
                    if (temp_diff.count() == 2) { create_ray(r1, s1, r2); continue; }
                    if (is_two_dimensional_face(trans, con_map, temp_diff)) {
                        create_ray(r1, s1, r2); 
                    }
                }
            }
            else {
                temp_diff.set_difference(local_supps[r2], r1_supp);
                if (temp_diff.singleton()) { create_ray(r1, s1, r2); }
            }
        }
#endif

    }
}

template <class T, class IndexSet>
inline
void
MatrixRayAlgorithm<T,IndexSet>::transfer(VectorArrayT<T>& rays, std::vector<IndexSet>& supps)
{
    rays.transfer(new_rays, 0, new_rays.get_number(), rays.get_number());
    supps.insert(supps.end(), new_supps.begin(), new_supps.end());
    new_supps.clear();
}

template <class T, class IndexSet>
inline
void
MatrixRayAlgorithm<T,IndexSet>::create_ray(Index r1, const T& s1, Index r2)
{
    //T s1; cone.get_slack(rays[r1], _next, s1);
    T s2; cone.get_slack(rays[r2], _next, s2);
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
MatrixCirAlgorithm<T,IndexSet>::MatrixCirAlgorithm(
                const ConeT<T>& _cone, const VectorArrayT<T>& _rays,
                const std::vector<IndexSet>& _supps, const std::vector<IndexSet>& _cir_supps, 
                const IndexSet& _rel, const Index& _cons_added)
        : MatrixRayAlgorithm<T,IndexSet>(_cone, _rays, _supps, _rel, _cons_added),
          cir_supps(_cir_supps)
{
}

template <class T, class IndexSet>
MatrixCirAlgorithm<T,IndexSet>::~MatrixCirAlgorithm()
{
}

template <class T, class IndexSet>
inline
void
MatrixCirAlgorithm<T,IndexSet>::compute()
{
    compute(MatrixRayAlgorithm<T,IndexSet>::_next, MatrixRayAlgorithm<T,IndexSet>::_r1_start, 
            MatrixRayAlgorithm<T,IndexSet>::_r1_end, MatrixRayAlgorithm<T,IndexSet>::_r2_start, 
            MatrixRayAlgorithm<T,IndexSet>::_r2_end);
}

template <class T, class IndexSet>
inline
void
MatrixCirAlgorithm<T,IndexSet>::threaded_compute(Index next, Index r1_start, Index r1_end, Index r2_start, Index r2_end)
{
    MatrixRayAlgorithm<T,IndexSet>::_next = next; 
    MatrixRayAlgorithm<T,IndexSet>::_r1_start = r1_start; 
    MatrixRayAlgorithm<T,IndexSet>::_r1_end = r1_end; 
    MatrixRayAlgorithm<T,IndexSet>::_r2_start = r2_start; 
    MatrixRayAlgorithm<T,IndexSet>::_r2_end = r2_end;
    ThreadedAlgorithm::threaded_compute();
}

template <class T, class IndexSet>
void
MatrixCirAlgorithm<T,IndexSet>::compute(Index next, Index r1_start, Index r1_end, Index r2_start, Index r2_end)
{
    // Copy class variables onto the stack.
    const ConeT<T>& cone = MatrixRayAlgorithm<T,IndexSet>::cone;
    const VectorArrayT<T>& rays = MatrixRayAlgorithm<T,IndexSet>::rays;
    const std::vector<IndexSet>& supps = MatrixRayAlgorithm<T,IndexSet>::supps;
    const IndexSet rel = MatrixRayAlgorithm<T,IndexSet>::rel;
    const Index cons_added = MatrixRayAlgorithm<T,IndexSet>::cons_added;
    MatrixRayAlgorithm<T,IndexSet>::_next = next; //TODO: This is a hack.

    char buffer[256];
    sprintf(buffer, "  Left = %3d  Col = %3d", rel.count(), next);

    if (r1_start == r1_end || r2_start == r2_end) { return; }
    DEBUG_4ti2(*out << "\nComputing circuits for ranges ";)
    DEBUG_4ti2(*out << "R1 [" << r1_start << "..." << r1_end << "] and ";)
    DEBUG_4ti2(*out << "R2 [" << r2_start << "..." << r2_end << "]\n";)

    MatrixT<T> trans;
    Size s = supps[0].get_size();
    IndexSet temp_supp(s);
    IndexSet temp_zeros(s);
    IndexSet r1_supp(s);
    Size r1_count;
    T s1;

    Size t = cir_supps[0].get_size();
    IndexSet r1_cir_supp(t);

    Index index_count = 0;
    for (Index r1 = r1_start; r1 < r1_end; ++r1) {
        // Output Statistics
        if (index_count % Globals::output_freq == 0) {
           *out << ENDL << buffer;
           *out << "  Size = " << std::setw(8) << rays.get_number();
           *out << ",  Index = " << r1 << "/" << r2_end << std::flush;
            DEBUG_4ti2(*out << std::endl;)
        }
        ++index_count;

        if (r2_start == r1) { ++r2_start; }

        r1_supp = supps[r1];
        r1_count = r1_supp.count();
        r1_cir_supp = cir_supps[r1];
        cone.get_slack(rays[r1], next, s1);

        if (r1_count == cons_added+1) {
            for (Index r2 = r2_start; r2 < r2_end; ++r2) {
                if (supps[r2].singleton_diff(r1_supp)
                    && r1_cir_supp.set_disjoint(cir_supps[r2])) {
                    create_circuit(r1, s1, r2);
                }
            }
            continue;
        }

        temp_supp.set_union(r1_supp, rel);
        temp_supp.set_complement();

        std::vector<Index> con_map;
        project_cone(cone, temp_supp, trans, con_map);
        DEBUG_4ti2(*out << "NEW TRANS:\n" << trans << "\n";)
        zero_cols(trans, con_map, temp_zeros);
        DEBUG_4ti2(*out << "NEW ZEROS:\n" << temp_zeros << "\n";)

        for (Index r2 = r2_start; r2 < r2_end; ++r2) {
            if (temp_zeros.singleton_intersection(supps[r2])
                && r1_supp.count_union(supps[r2]) <= cons_added+2
                && r1_cir_supp.set_disjoint(cir_supps[r2])) { 
                //if (supps[r2].singleton_diff(r1_supp)) { create_circuit(r1, s1, r2); continue; }
                //if (r1_supp.singleton_diff(supps[r2])) { create_circuit(r1, s1, r2); continue; }
                temp_supp.set_difference(supps[r2], r1_supp);
                //if (temp_supp.count_lte(2))) { create_circuit(r1, s1, r2); continue; }
                if (is_two_dimensional_face(trans, con_map, temp_supp)) {
                    create_circuit(r1, s1, r2);
                }
            }
        }
    }
    *out << ENDL << buffer;
    *out << "  Size = " << std::setw(8) << rays.get_number();
    *out << ",  Index = " << r1_end << "/" << r2_end << std::flush;
    DEBUG_4ti2(*out << std::endl;)
}

template <class T, class IndexSet>
void
MatrixCirAlgorithm<T,IndexSet>::create_circuit(Index i1, const T& s1, Index i2)
{
    DEBUG_4ti2(*out << "Creating new circuit.\n";)
    const VectorR<T>& r1 = MatrixRayAlgorithm<T,IndexSet>::rays[i1];
    const VectorR<T>& r2 = MatrixRayAlgorithm<T,IndexSet>::rays[i2];
    const Index next = MatrixRayAlgorithm<T,IndexSet>::_next;
    VectorT<T>& temp = MatrixRayAlgorithm<T,IndexSet>::temp;

    //T s1; MatrixRayAlgorithm<T,IndexSet>::cone.get_slack(r1, next, s1); 
    T s2; MatrixRayAlgorithm<T,IndexSet>::cone.get_slack(r2, next, s2); 
    if (s2 > 0) { temp.add(r1,s2,r2,-s1); }
    else { temp.add(r1,-s2,r2,s1); }
    temp.normalise();
    MatrixRayAlgorithm<T,IndexSet>::new_rays.insert(temp);

    const IndexSet& r1_supp = MatrixRayAlgorithm<T,IndexSet>::supps[i1];
    const IndexSet& r2_supp = MatrixRayAlgorithm<T,IndexSet>::supps[i2];
    IndexSet temp_supp(r1_supp.get_size());
    temp_supp.set_union(r1_supp, r2_supp);
    MatrixRayAlgorithm<T,IndexSet>::new_supps.push_back(temp_supp);

    IndexSet cir_supp1(cir_supps[i1]);
    IndexSet cir_supp2(cir_supps[i2]);
    cir_supp2.swap_odd_n_even();
    cir_supp1.set_union(cir_supp2);
    if (s1 < 0) { cir_supp1.swap_odd_n_even(); }
    new_cir_supps.push_back(cir_supp1);

    DEBUG_4ti2(*out << "Ray" << i1 << " " << r1 << "\t";)
    DEBUG_4ti2(*out << "Sup " << supps[i1] << "\t";)
    DEBUG_4ti2(*out << "Sla " << s1 << "\n";)
    DEBUG_4ti2(*out << "Ray" << i2 << " " << r2 << "\t";)
    DEBUG_4ti2(*out << "Sup " << supps[i2] << "\t";)
    DEBUG_4ti2(*out << "Sla " << s2 << "\n";)
    DEBUG_4ti2(*out << "Ray " << temp << "\t";)
    DEBUG_4ti2(*out << "Sup " << temp_supp << "\n";)
}

template <class T, class IndexSet>
inline
void
MatrixCirAlgorithm<T,IndexSet>::transfer(
            VectorArrayT<T>& rays,
            std::vector<IndexSet>& supps,
            std::vector<IndexSet>& cir_supps)
{
    MatrixRayAlgorithm<T,IndexSet>::transfer(rays, supps);
    cir_supps.insert(cir_supps.end(), new_cir_supps.begin(), new_cir_supps.end());
    new_cir_supps.clear();
}

#undef DEBUG_4ti2
#define DEBUG_4ti2(X) //X
