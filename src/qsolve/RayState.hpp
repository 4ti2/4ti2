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
#define DEBUG_4ti2(X)

using namespace _4ti2_;

template <class T, class IndexSet>
RayState<T,IndexSet>::RayState(const ConeT<T>& _cone, VectorArrayT<T>& _rays, std::vector<IndexSet>& _supps, const IndexSet& _rem, const IndexSet& _ray_mask, const Index& _next)
        : cone(_cone), rays(_rays), supps(_supps), rem(_rem), next(_next), ray_mask(_ray_mask),
          new_rays(0,rays.get_size()), temp(rays.get_size())
{
}

template <class T, class IndexSet>
RayState<T,IndexSet>::~RayState()
{
}

template <class T, class IndexSet>
RayState<T,IndexSet>*
RayState<T,IndexSet>::clone()
{
    return new RayState(cone, rays, supps, rem, ray_mask, next);
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
RayState<T,IndexSet>::sort_nonzeros(Index start, Index end, Index next_col, Index& middle)
{
    assert(start >= 0 && start <= end && end <= rays.get_number());
    T slack;
    Index index = start;
    for (Index i = start; i < end; ++i) {
        cone.get_slack(rays[i], next_col, slack);
        if (slack != 0) {
            rays.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes positives to the beginning.
template <class T, class IndexSet>
void
RayState<T,IndexSet>::sort_positives(Index start, Index end, Index next_col, Index& middle)
{
    assert(start >= 0 && start <= end && end <= rays.get_number());
    T slack;
    Index index = start;
    for (Index i = start; i < end; ++i) {
        cone.get_slack(rays[i], next_col, slack);
        if (slack > 0) {
            rays.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes positives to the beginning.
template <class T, class IndexSet>
void
RayState<T,IndexSet>::sort_negatives(Index start, Index end, Index next_col, Index& middle)
{
    assert(start >= 0 && start <= end && end <= rays.get_number());
    T slack;
    Index index = start;
    for (Index i = start; i < end; ++i) {
        cone.get_slack(rays[i], next_col, slack);
        if (slack < 0) {
            rays.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            ++index;
        }
    }
    middle = index;
}

// Pushes rays to the beginning.
template <class T, class IndexSet>
void
RayState<T,IndexSet>::sort_filter(const IndexSet& filter, Index start, Index end, Index& middle)
{
    Index index = start;
    for (Index i = start; i < end; ++i) {
        if (!filter.set_disjoint(supps[i])) {
            rays.swap_vectors(i,index);
            IndexSet::swap(supps[i], supps[index]);
            ++index;
        }
    }
    middle = index;
}

template <class T, class IndexSet>
void
RayState<T,IndexSet>::update(Index index, Index start, Index end)
{
    // Update supports
    for (Index i = start; i < end; ++i) { supps[i].set(index); }
}

template <class T, class IndexSet>
void
RayState<T,IndexSet>::remove(Index start, Index end)
{
    supps.erase(supps.begin()+start, supps.begin()+end);
    rays.remove(start, end);
}

template <class T, class IndexSet>
void
RayState<T,IndexSet>::flip(Index start, Index end)
{
    for (Index i = start; i < end; ++i) { supps[i].swap_odd_n_even(); }
}

template <class T, class IndexSet>
void
RayState<T,IndexSet>::resize(Size size)
{
    for (Index i = 0; i < (Index) supps.size(); ++i) { supps[i].resize(size); }
}

template <class T, class IndexSet>
Index
RayState<T,IndexSet>::next_constraint(const ConsOrder& order, Index& pos_start, Index& pos_end, Index& neg_start, Index& neg_end)
{
    assert(!rem.empty());
    // First, we choose the next constraint to add.
    Index next_con = next_constraint(order, rem);

    // TODO: Should we use a vector of slacks?
    Index start = 0; Index end = rays.get_number(); Index middle;

    // We sort the vectors into nonzeros and then zeros.
    sort_nonzeros(start, end, next_con, middle);
    Index nonzero_start = start, nonzero_end = middle;

    // We sort the rays into positives and then negatives.
    sort_positives(nonzero_start, nonzero_end, next_con, middle);
    pos_start = nonzero_start; pos_end = middle;
    neg_start = middle; neg_end = nonzero_end;

    return next_con;
}

template <class T, class IndexSet>
Index
RayState<T,IndexSet>::next_constraint(
                const ConsOrder& order,
                Index& pos_ray_start, Index& pos_ray_end, Index& neg_ray_start, Index& neg_ray_end,
                Index& pos_cir_start, Index& pos_cir_end, Index& neg_cir_start, Index& neg_cir_end)
{
    assert(!rem.empty());
    // First, we choose the next constraint to add.
    //TODO Index next_col = next_circuit_constraint(order, rem);
    Index next_col = next_constraint(order, rem);
    
    Index start = 0; Index end = rays.get_number(); Index middle;

    // We sort the vectors into nonzeros and then zeros.
    sort_nonzeros(start, end, next_col, middle);
    Index nonzero_start = start, nonzero_end = middle;

    // We sort the nonzeros into positives and then negatives.
    sort_positives(nonzero_start, nonzero_end, next_col, middle);
    Index pos_start = nonzero_start; Index pos_end = middle;
    Index neg_start = middle; Index neg_end = nonzero_end;

    // We sort the positives into rays and then circuits.
    sort_filter(ray_mask, pos_start, pos_end, middle);
    pos_ray_start = pos_start; pos_ray_end = middle;
    pos_cir_start = middle; pos_cir_end = pos_end;

    // We sort the negatives into circuits then rays.
    IndexSet cir_mask(ray_mask); cir_mask.set_complement();
    sort_filter(cir_mask, neg_start, neg_end, middle);
    neg_cir_start = neg_start; neg_cir_end = middle;
    neg_ray_start = middle; neg_ray_end = neg_end;

    return next_col;
}

template <class T, class IndexSet>
Index
RayState<T,IndexSet>::sort_count(Size count, Index start, Index end)
{
    // We sort the r2's into vectors where r2_supp.count()==cons_added+1.
    Index middle = start;
    for (Index i = start; i < end; ++i) {
        if (supps[i].count() == count) {
            rays.swap_rows(i, middle);
            IndexSet::swap(supps[i], supps[middle]);
            ++middle;
        }
    }
    return middle;
}

template <class T, class IndexSet>
inline
void
RayState<T,IndexSet>::transfer()
{
    rays.transfer(new_rays, 0, new_rays.get_number(), rays.get_number());
    supps.insert(supps.end(), new_supps.begin(), new_supps.end());
    new_supps.clear();
}

template <class T, class IndexSet>
inline  
void    
RayState<T,IndexSet>::set_r1_index(Index r1)
{
    _r1 = r1;
    cone.get_slack(rays[r1], next, s1);
}

template <class T, class IndexSet>
void
RayState<T,IndexSet>::create_ray(Index r2)
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

    DEBUG_4ti2(
    *out << "\nADDING VECTOR.\n";
    *out << "R1: " << _r1 << "\n";
    *out << rays[_r1] << "\n";
    *out << "R2: " << r2 << "\n";
    *out << rays[r2] << "\n";
    *out << "NEW:\n";
    *out << temp << "\n";
    )
}

template <class T, class IndexSet>
void
RayState<T,IndexSet>::create_circuit(Index i2)
{
    DEBUG_4ti2(*out << "Creating new circuit.\n";)

    IndexSet tmp_union(supps[i2]);
    tmp_union.set_union(supps[_r1]);
    new_supps.push_back(tmp_union);

    DEBUG_4ti2(
        *out << "Cir1 " << supps[i1] << "\n";
        *out << "Cir2 " << supps[i2] << "\n";
        *out << "Cir0 " << tmp_union << "\n";
    )

    DEBUG_4ti2(*out << "Creating new circuit.\n";)
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

    DEBUG_4ti2(
        *out << "Ray1 " << i1 << " " << s1 << " : " << r1 << "\n";
        *out << "Ray2 " << i2 << " " << s2 << " : " << r2 << "\n";
        *out << "Ray0 " << temp << "\n";
    )
}


template <class T, class IndexSet>
Index
RayState<T,IndexSet>::next_constraint(const ConsOrder& order, const IndexSet& rem)
{
    // Sanity Check
    assert(rays.get_size() == rem.get_size());

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
    DEBUG_4ti2(*out << "Next Constraint is " << next_con << "\n";)
    return next_con;
}

template <class T, class IndexSet>
inline
void
RayState<T,IndexSet>::project_cone(
                const IndexSet& zero_supp,
                std::vector<Index>& con_map,
                IndexSet& zeros)
{
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
    DEBUG_4ti2(*out << "Zeros:\n" << zeros << "\n";)
}

// Checks whether the given support determines a two dimensional face of the cone.
template <class T, class IndexSet>
inline
bool
RayState<T,IndexSet>::is_two_dimensional_face(
            const std::vector<Index>& con_map,
            const IndexSet& diff)
{
    DEBUG_4ti2(*out << "\nis_two_dimensional_face\n";)
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
    DEBUG_4ti2(*out << "Vars:\n" << vars << "\nCons:\n" << cons << "\n";)

    return sub_cone.is_one_dimensional_face(vars, cons);
}

#if 0
template <class T>
MatrixHelper<T>::MatrixHelper(
                const ConeT<T>& _cone, VectorArrayT<T>& _rays, const Index& _next)
        : cone(_cone), rays(_rays), next(_next), new_rays(0,_rays.get_size()), temp(_cone.num_vars())
{
}

template <class T>
MatrixHelper<T>*
MatrixHelper<T>::clone()
{
    return new MatrixHelper<T>(cone, rays, next);
}

template <class T>
MatrixHelper<T>::~MatrixHelper()
{
}

template <class T>
inline
void
MatrixHelper<T>::project_cone(
                const IndexSetD& zero_supp,
                std::vector<Index>& con_map,
                IndexSetD& zeros)
{
    project_coneT(zero_supp, con_map, zeros);
}

template <class T>
inline
void
MatrixHelper<T>::project_cone(
                const IndexSetDS& zero_supp,
                std::vector<Index>& con_map,
                IndexSetDS& zeros)
{
    project_coneT(zero_supp, con_map, zeros);
}

template <class T> template <class IndexSet>
inline
void
MatrixHelper<T>::project_coneT(
                const IndexSet& zero_supp,
                std::vector<Index>& con_map,
                IndexSet& zeros)
{
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
    DEBUG_4ti2(*out << "Zeros:\n" << zeros << "\n";)
}

template <class T>
bool
MatrixHelper<T>::is_two_dimensional_face(
            const std::vector<Index>& con_map,
            const IndexSetD& diff)
{
    return is_two_dimensional_faceT(con_map, diff);
}

template <class T>
bool
MatrixHelper<T>::is_two_dimensional_face(
            const std::vector<Index>& con_map,
            const IndexSetDS& diff)
{
    return is_two_dimensional_faceT(con_map, diff);
}

// Checks whether the given support determines a two dimensional face of the cone.
template <class T> template <class IndexSet>
inline
bool
MatrixHelper<T>::is_two_dimensional_faceT(
            const std::vector<Index>& con_map,
            const IndexSet& diff)
{
    DEBUG_4ti2(*out << "\nis_two_dimensional_face\n";)
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
    DEBUG_4ti2(*out << "Vars:\n" << vars << "\nCons:\n" << cons << "\n";)

    return sub_cone.is_one_dimensional_face(vars, cons);
}

template <class T>
inline
void
MatrixHelper<T>::transfer()
{
    rays.transfer(new_rays, 0, new_rays.get_number(), rays.get_number());
}

template <class T>
inline  
void    
MatrixHelper<T>::set_r1_index(Index r1)
{
    _r1 = r1;
    cone.get_slack(rays[r1], next, s1);
}

template <class T>
void
MatrixHelper<T>::create_rays(Index i1, std::vector<Index>& r2s)
{
    VectorR<T>& r1 = rays[i1];
    T s1, s2; 
    cone.get_slack(r1, next, s1);
    for (std::vector<Index>::iterator i2 = r2s.begin(); i2 != r2s.end(); ++i2) {
        VectorR<T>& r2 = rays[*i2];
        cone.get_slack(r2, next, s2);
        if (i1 < *i2) { temp.add(r2, s1, r1, -s2); }
        else { temp.add(r1, s2, r2, -s1); }

        temp.normalise();
        new_rays.insert(temp);
    }
}

template <class T>
void
MatrixHelper<T>::create_ray(Index r2)
{
    //T s1; cone.get_slack(rays[r1], next, s1);
    T s2; cone.get_slack(rays[r2], next, s2);
    if (_r1 < r2) { temp.add(rays[r2], s1, rays[_r1], -s2); }
    else { temp.add(rays[_r1], s2, rays[r2], -s1); }

    temp.normalise();
    new_rays.insert(temp);

    DEBUG_4ti2(
    *out << "\nADDING VECTOR.\n";
    *out << "R1: " << _r1 << "\n";
    *out << rays[_r1] << "\n";
    *out << "R2: " << r2 << "\n";
    *out << rays[r2] << "\n";
    *out << "NEW:\n";
    *out << temp << "\n";
    )
}

template <class T>
bool
MatrixHelper<T>::create_circuit(Index i2)
{
    DEBUG_4ti2(*out << "Creating new circuit.\n";)
    const Index i1 = _r1;
    const VectorR<T>& r1 = rays[i1];
    const VectorR<T>& r2 = rays[i2];

    //T s1; cone.get_slack(r1, next, s1); 
    T s2; cone.get_slack(r2, next, s2); 
    if (s1 > 0) { temp.add(r1,-s2,r2,-s1); }
    else { temp.add(r1,-s2,r2,s1); }
    temp.normalise();
    new_rays.insert(temp);

    if (s1 > 0) { return true; }
    return false;

    DEBUG_4ti2(
        *out << "Ray1 " << i1 << " " << s1 << " : " << r1 << "\n";
        *out << "Ray2 " << i2 << " " << s2 << " : " << r2 << "\n";
        *out << "Ray0 " << temp << "\n";
    )
}

template <class T>
SupportHelper<T>::SupportHelper(const ConeT<T>& _cone, VectorArrayT<T>& _rays, const Index& _next)
        : cone(_cone), rays(_rays), new_rays(0,_rays.get_size()), next(_next), temp(_cone.num_vars())
{
}

template <class T>
SupportHelper<T>::~SupportHelper()
{
}

template <class T>
inline
void
SupportHelper<T>::transfer()
{
    rays.transfer(new_rays, 0, new_rays.get_number(), rays.get_number());
}

template <class T>
inline
void
SupportHelper<T>::set_r1_index(Index _r1)
{
    r1 = _r1;
    cone.get_slack(rays[r1], next, s1);
}

template <class T>
inline
void
SupportHelper<T>::create_ray(Index r2)
{
    T s2; cone.get_slack(rays[r2], next, s2);
    if (r1 < r2) { temp.add(rays[r2], s1, rays[r1], -s2); }
    else { temp.add(rays[r1], s2, rays[r2], -s1); }
    temp.normalise();
    new_rays.insert(temp);

    DEBUG_4ti2(
    *out << rays[r1] << "\n";
    *out << rays[r2] << "\n";
    )
}

template <class T>
bool
SupportHelper<T>::create_circuit(Index i1, Index i2)
{
    DEBUG_4ti2(*out << "Creating new circuit.\n";)
    const VectorR<T>& r1 = rays[i1];
    const VectorR<T>& r2 = rays[i2];

    //T s1; cone.get_slack(r1, next, s1); 
    T s2; cone.get_slack(r2, next, s2); 

    if (s2 > 0) { temp.add(r1,s2,r2,-s1); }
    else { temp.add(r1,-s2,r2,s1); }
    temp.normalise();
    new_rays.insert(temp);

    if (s1 > 0) { return true; }
    return false;

    DEBUG_4ti2(
        *out << "Ray1 " << i1 << " " << s1 << " : " << r1 << "\n";
        *out << "Ray2 " << i2 << " " << s2 << " : " << r2 << "\n";
        *out << "Ray0 " << temp << "\n";
    )
}
#endif

#undef DEBUG_4ti2

