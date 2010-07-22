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

#include "qsolve/Cone.h"
#include "qsolve/Globals.h"
#include "qsolve/HermiteAlgorithm.h"
#include "qsolve/VectorStream.h"
#include "qsolve/Debug.h"
#include "qsolve/Stream.h"
#include "qsolve/IndexSetD.h"
#include "qsolve/IndexSetStream.h"
#include "4ti2/4ti2.h"

#undef DEBUG_4ti2
#define DEBUG_4ti2 0

using namespace _4ti2_;

template <class IndexSet>
RayStateAPI<IndexSet>::RayStateAPI(const ConeAPI& _cone)
    : cone_api(_cone), next(-1), cons_added(0), 
      rem(_cone.num_vars()+_cone.num_cons(),0), rel(_cone.num_vars()+_cone.num_cons(),0), ray_mask(0)
{
}

template <class IndexSet>
RayStateAPI<IndexSet>::~RayStateAPI()
{
}

template <class IndexSet>
void
RayStateAPI<IndexSet>::update(Index index, Index start, Index end)
{
    // Update supports
    for (Index i = start; i < end; ++i) { supps[i].set(index); }
}

template <class IndexSet>
void
RayStateAPI<IndexSet>::flip(Index start, Index end)
{
    for (Index i = start; i < end; ++i) { supps[i].swap_odd_n_even(); }
}

template <class IndexSet>
void
RayStateAPI<IndexSet>::resize(Size size)
{
    for (Index i = 0; i < (Index) supps.size(); ++i) { supps[i].resize(size); }
}

template <class T, class IndexSet>
RayState<T,IndexSet>::RayState(const ConeT<T>& _cone, VectorArrayT<T>& _rays)
        : RayStateAPI<IndexSet>(_cone), cone(_cone), rays(_rays)
{
}

template <class T, class IndexSet>
RayState<T,IndexSet>::~RayState()
{
}

template <class T, class IndexSet>
RaySubState<T,IndexSet>*
RayState<T,IndexSet>::clone()
{
    return new RaySubState<T,IndexSet>(*this, cone, rays, Base::supps, Base::next);
}

template <class T, class IndexSet>
Size
RayState<T,IndexSet>::num_vars()
{
    return cone.num_vars();
}

template <class T, class IndexSet>
Size
RayState<T,IndexSet>::num_cons()
{
    return cone.num_cons();
}

template <class T, class IndexSet>
Size
RayState<T,IndexSet>::num_gens()
{
    return rays.get_number();
}

// Pushes zeros to the beginning.
template <class T, class IndexSet>
void
RayState<T,IndexSet>::sort_nonzeros(Index start, Index end, Index& middle)
{
    assert(start >= 0 && start <= end && end <= rays.get_number());
    T slack;
    Index index = start;
    for (Index i = start; i < end; ++i) {
        cone.get_slack(rays[i], Base::next, slack);
        if (slack != 0) {
            rays.swap_vectors(i,index);
            IndexSet::swap(Base::supps[i], Base::supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes positives to the beginning.
template <class T, class IndexSet>
void
RayState<T,IndexSet>::sort_positives(Index start, Index end, Index& middle)
{
    assert(start >= 0 && start <= end && end <= rays.get_number());
    T slack;
    Index index = start;
    for (Index i = start; i < end; ++i) {
        cone.get_slack(rays[i], Base::next, slack);
        if (slack > 0) {
            rays.swap_vectors(i,index);
            IndexSet::swap(Base::supps[i], Base::supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes positives to the beginning.
template <class T, class IndexSet>
void
RayState<T,IndexSet>::sort_negatives(Index start, Index end, Index& middle)
{
    assert(start >= 0 && start <= end && end <= rays.get_number());
    T slack;
    Index index = start;
    for (Index i = start; i < end; ++i) {
        cone.get_slack(rays[i], Base::next, slack);
        if (slack < 0) {
            rays.swap_vectors(i,index);
            IndexSet::swap(Base::supps[i], Base::supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes rays to the beginning.
template <class T, class IndexSet>
void
RayState<T,IndexSet>::sort_filter(const IndexSet& filter, bool pos, Index start, Index end, Index& middle)
{
    Index index = start;
    for (Index i = start; i < end; ++i) {
        if (filter.set_disjoint(Base::supps[i]) != pos) {
            rays.swap_vectors(i,index);
            IndexSet::swap(Base::supps[i], Base::supps[index]);
            ++index;
        }
    }
    middle = index;
}

template <class T, class IndexSet>
void
RayState<T,IndexSet>::remove(Index start, Index end)
{
    Base::supps.erase(Base::supps.begin()+start, Base::supps.begin()+end);
    rays.remove(start, end);
}

template <class T, class IndexSet>
Index
RayState<T,IndexSet>::next_constraint(const ConsOrder& order, const IndexSet& rem, Index& pos_start, Index& pos_end, Index& neg_start, Index& neg_end)
{
    assert(!rem.empty());
    // First, we choose the next constraint to add.
    Base::next = next_constraint(order, rem);

    // TODO: Should we use a vector of slacks?
    Index start = 0; Index end = rays.get_number(); Index middle;

    // We sort the vectors into nonzeros and then zeros.
    sort_nonzeros(start, end, middle);
    Index nonzero_start = start, nonzero_end = middle;

    // We sort the rays into positives and then negatives.
    sort_positives(nonzero_start, nonzero_end, middle);
    pos_start = nonzero_start; pos_end = middle;
    neg_start = middle; neg_end = nonzero_end;

    return Base::next;
}

template <class T, class IndexSet>
Index
RayState<T,IndexSet>::next_constraint(
                const ConsOrder& order, const IndexSet& rem, const IndexSet& ray_mask,
                Index& pos_ray_start, Index& pos_ray_end, Index& neg_ray_start, Index& neg_ray_end,
                Index& pos_cir_start, Index& pos_cir_end, Index& neg_cir_start, Index& neg_cir_end)
{
    assert(!rem.empty());
    // First, we choose the next constraint to add.
    //TODO Index Base::next = next_circuit_constraint(order, rem);
    Base::next = next_constraint(order, rem);
    
    Index start = 0; Index end = rays.get_number(); Index middle;

    // We sort the vectors into nonzeros and then zeros.
    sort_nonzeros(start, end, middle);
    Index nonzero_start = start, nonzero_end = middle;

    // We sort the nonzeros into positives and then negatives.
    sort_positives(nonzero_start, nonzero_end, middle);
    Index pos_start = nonzero_start; Index pos_end = middle;
    Index neg_start = middle; Index neg_end = nonzero_end;

    // We sort the positives into rays and then circuits.
    sort_filter(ray_mask, true, pos_start, pos_end, middle);
    pos_ray_start = pos_start; pos_ray_end = middle;
    pos_cir_start = middle; pos_cir_end = pos_end;

    // We sort the negatives into circuits then rays.
    sort_filter(ray_mask, false, neg_start, neg_end, middle);
    neg_cir_start = neg_start; neg_cir_end = middle;
    neg_ray_start = middle; neg_ray_end = neg_end;

#if DEBUG_4ti2 > 0
    *out << "\nr+ " << pos_ray_start << "," << pos_ray_end << " r- " << neg_ray_start << "," << neg_ray_end;
    *out << " c+ " << pos_cir_start << "," << pos_cir_end << " c-" << neg_cir_start << "," << neg_cir_end << std::endl;
#endif

    return Base::next;
}

template <class T, class IndexSet>
Index
RayState<T,IndexSet>::sort_count(Size count, Index start, Index end)
{
    // We sort the r2's into vectors where r2_supp.count()==cons_added+1.
    Index middle = start;
    for (Index i = start; i < end; ++i) {
        if (Base::supps[i].count() == count) {
            rays.swap_rows(i, middle);
            IndexSet::swap(Base::supps[i], Base::supps[middle]);
            ++middle;
        }
    }
    return middle;
}

template <class T, class IndexSet>
Index
RayState<T,IndexSet>::next_constraint(const ConsOrder& order, const IndexSet& rem)
{
    Index next_con = *rem.begin();
    if (order.get_constraint_order() == MININDEX) { return next_con; }

    typename IndexSet::Iter it = rem.begin(); ++it;
    if (it == rem.end()) { return next_con; }

    Size next_pos_count, next_neg_count, next_zero_count;
    cone.slack_count(rays, next_con, next_pos_count, next_neg_count, next_zero_count);
    while (it != rem.end()) {
        Size pos_count, neg_count, zero_count;
        cone.slack_count(rays, *it, pos_count, neg_count, zero_count);
        if ((*order.compare)(next_pos_count, next_neg_count, next_zero_count,
                        pos_count, neg_count, zero_count)) {
            next_con = *it;
            next_pos_count = pos_count;
            next_neg_count = neg_count;
            next_zero_count = zero_count;
        }
        ++it;
    }
    return next_con;
}

template <class T, class IndexSet>
void
RayState<T, IndexSet>::check()
{
    bool failed = false;
    for (Index i = 0; i < rays.get_number(); ++i) {
        failed |= check(rays[i], Base::supps[i]);
    }

    if (failed) { 
        *out << "supps_to_cons: ";
        for (Index j = 0; j < (Index) Base::supps_to_cons.size(); ++j) { *out << " " << Base::supps_to_cons[j]; }
        *out << "\ncons_to_supps: ";
        for (Index j = 0; j < (Index) Base::cons_to_supps.size(); ++j) { *out << " " << Base::cons_to_supps[j]; }
        *out << "\nsupp_types:    ";
        for (Index j = 0; j < (Index) Base::supp_types.size(); ++j) { *out << " " << Base::supp_types[j]; }
        *out << "\n";
    }
}

template <class T, class IndexSet>
bool
RayState<T, IndexSet>::check(const VectorR<T>& ray, const IndexSet& supp)
{
    Index n = cone.num_vars();
    Index m = cone.num_cons();

    VectorT<T> slacks(n+m);
    bool failed = false;
    cone.get_slacks(ray, slacks);
    for (Index j = 0; j < (Index) Base::supps_to_cons.size(); ++j) {
        if (Base::supp_types[j] == _4ti2_LB) {
            if ((slacks[Base::supps_to_cons[j]] > 0) != supp[j]) { 
                *out << "Check LB Supp failed: " << supp[j] << " " << Base::supps_to_cons[j] << " " << slacks[Base::supps_to_cons[j]] << "\n";
                failed = true;
            }
            if (cone.get_constraint_type(Base::supps_to_cons[j]) == _4ti2_LB && slacks[Base::supps_to_cons[j]] < 0) { 
                *out << "Check LB failed.\n"; failed = true;
            }
        } 
        else if (Base::supp_types[j] == _4ti2_UB) {
            if ((slacks[Base::supps_to_cons[j]] < 0) != supp[j]) { 
                *out << "Check UB Supp failed: " << supp[j] << " " << Base::supps_to_cons[j] << " " << slacks[Base::supps_to_cons[j]] << "\n";
                failed = true;
            }
            if (cone.get_constraint_type(Base::supps_to_cons[j]) == _4ti2_UB && slacks[Base::supps_to_cons[j]] > 0) {
                * out << "Check UB failed.\n"; failed = true;
            }
        }
        else if (Base::supp_types[j] == _4ti2_EQ) {
            if (slacks[Base::supps_to_cons[j]] != 0) { *out << "Check EQ failed.\n"; failed = true; }
            if (supp[j]) { *out << "Check EQ Supps failed.\n"; failed = true; }
        }
        else if (Base::supp_types[j] == _4ti2_DB) {
            if (slacks[Base::supps_to_cons[j]] == 0 && supp[j]) { *out << "Support Check DB failed.\n"; failed = true; }
        }
        else if (Base::supp_types[j] == _4ti2_FR) {
            if (supp[j]) { *out << "Check FR Supps failed.\n"; failed = true; failed = true; }
        }
    }
    if (failed) {
        *out << "Ray:    " << ray << "\nSupp:   " << supp << "\nSlacks: " << slacks << "\n\n"; 
    }
    return failed;
}

/////////////////
// RaySubState //
/////////////////

template <class T, class IndexSet>
RaySubState<T,IndexSet>::RaySubState(RayState<T,IndexSet>& _state, const ConeT<T>& _cone, VectorArrayT<T>& _rays, std::vector<IndexSet>& _supps, const Index& _next)
        : RaySubStateAPI<IndexSet>(), supps(_supps), state(_state), cone(_cone), rays(_rays), next(_next),
          new_rays(0,rays.get_size()), temp(rays.get_size())
{
}

template <class T, class IndexSet>
RaySubState<T,IndexSet>::~RaySubState()
{
}

template <class T, class IndexSet>
inline
void
RaySubState<T,IndexSet>::transfer()
{
    rays.transfer(new_rays, 0, new_rays.get_number(), rays.get_number());
    supps.insert(supps.end(), new_supps.begin(), new_supps.end());
    new_supps.clear();
}

template <class T, class IndexSet>
inline  
void    
RaySubState<T,IndexSet>::set_r1_index(Index r1)
{
    _r1 = r1;
    cone.get_slack(rays[r1], next, s1);
}

template <class T, class IndexSet>
void
RaySubState<T,IndexSet>::create_ray(Index r2)
{
    IndexSet temp_supp(supps[_r1]);
    temp_supp.set_union(supps[r2]);
    new_supps.push_back(temp_supp);

    //T s1; cone.get_slack(rays[r1], next, s1);
    T s2; cone.get_slack(rays[r2], next, s2);
    if (_r1 < r2) { temp.add(rays[r2], s1, rays[_r1], -s2); }
    else { temp.add(rays[_r1], s2, rays[r2], -s1); }
    temp.normalise();
    new_rays.insert(temp);

#if DEBUG_4ti2 > 0
    *out << "\nADDING NEW RAY.\n";
    *out << "R1: " << _r1 << "\n";
    *out << rays[_r1] << "\n";
    *out << "R2: " << r2 << "\n";
    *out << rays[r2] << "\n";
    *out << "NEW:\n";
    *out << temp << "\n";
#endif
}

template <class T, class IndexSet>
void
RaySubState<T,IndexSet>::create_circuit(Index i2)
{
    IndexSet tmp_union(supps[i2]);
    tmp_union.set_union(supps[_r1]);
    new_supps.push_back(tmp_union);

    const VectorR<T>& r1 = rays[_r1];
    const VectorR<T>& r2 = rays[i2];

    //T s1; cone.get_slack(r1, next, s1); 
    T s2; cone.get_slack(r2, next, s2); 

    if (s1 > 0) {
        if (s2 > 0) { temp.add(r1,s2,r2,-s1); } // r1 - r2.
        else { temp.add(r1,-s2,r2,s1); } // r1 + r2.
    }
    else { 
        if (s2 > 0) { temp.add(r1,s2,r2,-s1); } // r1 + r2.
        else { temp.add(r1,s2,r2,-s1); } // -r1 + r2.
    } 

    temp.normalise();
    new_rays.insert(temp);

#if DEBUG_4ti2 > 0
    if (state.check(temp, tmp_union)) {
        *out << "Ray1 " << _r1 << " " << s1 << " : " << r1 << "\n";
        *out << "Ray2 " << i2 << " " << s2 << " : " << r2 << "\n";
        *out << "Ray0 " << temp << "\n";
        *out << "Cir1 " << supps[_r1] << "\n";
        *out << "Cir2 " << supps[i2] << "\n";
        *out << "Cir0 " << tmp_union << "\n";
    }
#endif
}

#if 0
template <class T, class IndexSet>
inline
void
RaySubState<T,IndexSet>::project_cone(
                IndexSet& zero_supp,
                std::vector<Index>& con_map,
                IndexSet& zeros)
{
    zero_supp.set_complement();
    zero_supp.set_difference(state.rel);
    sub_cone.project_cone(cone, zero_supp, con_map);

    const MatrixT<T>& trans = sub_cone.get_matrix();
    zeros.zero();
    Index n = trans.get_m();
    Index m = trans.get_n();
    for (Index i = 0; i < n; ++i) {
        if (con_map[i] < 0) { continue; }
        zeros.set(con_map[i]);
        for (Index j = 0; j < m; ++j) {
            // TODO: trans row.
            if (trans(i,j) != 0 && con_map[j+n] >= 0) { zeros.unset(con_map[i]); break; }
        }
    }
}
#endif

#if 1
template <class T, class IndexSet>
inline
void
RaySubState<T,IndexSet>::project_cone(
                IndexSet& temp_supp,
                std::vector<Index>& con_map,
                IndexSet& zeros)
{
    IndexSet zero_supp(cone.num_vars()+cone.num_cons(),0);
    for (typename IndexSet::Iter it = temp_supp.begin(); it != temp_supp.end(); ++it) {
        zero_supp.set(state.supps_to_cons[*it]);
    }

    zero_supp.set_complement();
    zero_supp.set_difference(state.rel);
    sub_cone.project_cone(cone, zero_supp, con_map);
    for (Index i = 0; i < (Index) con_map.size(); ++i) {
        if (con_map[i] != -1) { con_map[i] = state.cons_to_supps[con_map[i]]; }
    }

    const MatrixT<T>& trans = sub_cone.get_matrix();
    zeros.zero();
    Index n = trans.get_m();
    Index m = trans.get_n();
    for (Index i = 0; i < n; ++i) {
        if (con_map[i] < 0) { continue; }
        zeros.set(con_map[i]);
        for (Index j = 0; j < m; ++j) {
            // TODO: trans row.
            if (trans(i,j) != 0 && con_map[j+n] >= 0) { zeros.unset(con_map[i]); break; }
        }
    }
}
#endif

// Checks whether the given support determines a two dimensional face of the cone.
template <class T, class IndexSet>
inline
bool
RaySubState<T,IndexSet>::is_two_dimensional_face(
            const std::vector<Index>& con_map,
            const IndexSet& diff)
{
    Index n = sub_cone.num_vars();
    Index m = sub_cone.num_cons();

    IndexSet vars(n,0);
    for (Index i = 0; i < n; ++i) { 
        if (con_map[i] != -1 && diff[con_map[i]]) { vars.set(i); }
    }
    IndexSet cons(m,0);
    for (Index i = n; i < m+n; ++i) {
        if (con_map[i] != -1 && !diff[con_map[i]]) { cons.set(i-n); }
    }
#if DEBUG_4ti2 > 0
    *out << "\nis_two_dimensional_face\n";
    *out << "Vars:\n" << vars << "\nCons:\n" << cons << "\n";
#endif

    return sub_cone.is_one_dimensional_face(vars, cons);
}

#undef DEBUG_4ti2
