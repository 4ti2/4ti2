/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Peter Malkin.

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

#ifndef _4ti2_qsolve__RayState_
#define _4ti2_qsolve__RayState_

#include "qsolve/Vector.h"
#include "qsolve/VectorArray.h"
#include "qsolve/QSolveConsOrder.h"
#include "qsolve/Cone.h"
#include "qsolve/ConeC.h"

#include "qsolve/BinaryTree.h"
#include "qsolve/FullTree.h"
#include "qsolve/MultiTree.h"

//#define SUPPORTTREE FullTree  // TODO: DOES NOT WORK!
//#define SUPPORTTREE BinaryTree
#define SUPPORTTREE MultiTree

namespace _4ti2_
{

template <class IndexSet> class RaySubStateAPI;
template <class T, class IndexSet> class RaySubState;

template <class IndexSet>
class RayStateAPI
{
public:
    RayStateAPI(const ConeAPI& cone);
    virtual ~RayStateAPI();
    virtual RaySubStateAPI<IndexSet>* clone() = 0;

    virtual Size num_vars() = 0;
    virtual Size num_cons() = 0;
    virtual Size num_gens() = 0;

    virtual Index next_constraint(const ConsOrder& order, const IndexSet& rem, Index& pos_start, Index& pos_end, Index& neg_start, Index& neg_end) = 0;
    virtual Index next_constraint(const ConsOrder& order, const IndexSet& rem, const IndexSet& ray_mask,
                Index& pos_ray_start, Index& pos_ray_end, Index& neg_ray_start, Index& neg_ray_end,
                Index& pos_cir_start, Index& pos_cir_end, Index& neg_cir_start, Index& neg_cir_end) = 0;
    virtual Index sort_count(Size count, Index start, Index end) = 0;
    void update(Index index, Index start, Index end);
    void flip(Index start, Index end);
    void resize(Size s);
    virtual void remove(Index start, Index end) = 0;
    virtual void check() = 0;

public:
    std::vector<IndexSet> supps;
    SUPPORTTREE<IndexSet> tree;
    const ConeAPI& cone_api;
    Index next;
    Index cons_added;
    IndexSet rem;
    IndexSet rel;
    IndexSet ray_mask;

    std::vector<Index> supps_to_cons;
    std::vector<_4ti2_constraint> supp_types;
    std::vector<Index> cons_to_supps;
};

template <class T, class IndexSet>
class RayState : public RayStateAPI<IndexSet>
{
public:
    RayState(const ConeT<T>& _cone, VectorArrayT<T>& _rays);
    virtual ~RayState();
    virtual RaySubState<T,IndexSet>* clone();

    Size num_vars();
    Size num_cons();
    Size num_gens();

    Index next_constraint(const ConsOrder& order, const IndexSet& rem, Index& pos_start, Index& pos_end, Index& neg_start, Index& neg_end);
    Index next_constraint(const ConsOrder& order, const IndexSet& rem, const IndexSet& ray_mask,
                Index& pos_ray_start, Index& pos_ray_end, Index& neg_ray_start, Index& neg_ray_end,
                Index& pos_cir_start, Index& pos_cir_end, Index& neg_cir_start, Index& neg_cir_end);
    Index next_constraint(const ConsOrder& order, const IndexSet& rem);
    void check();
    bool check(const VectorR<T>& ray, const IndexSet& supp, Size d);

    Index sort_count(Size count, Index start, Index end);
    void remove(Index start, Index end);

public:
    const ConeT<T>& cone;
    VectorArrayT<T>& rays;

protected:
    typedef RayStateAPI<IndexSet> Base;

    void sort_nonzeros(Index start, Index end, Index& middle);
    void sort_positives(Index start, Index end, Index& middle);
    void sort_negatives(Index start, Index end, Index& middle);
    void sort_filter(const IndexSet& filter, bool pos, Index start, Index end, Index& middle);
};

template <class IndexSet>
class RaySubStateAPI
{
public:
    RaySubStateAPI() : r1_supp(0) {}
    virtual ~RaySubStateAPI() {}

    virtual void set_r1_index(Index r1) = 0;
    virtual void create_ray(Index r2) = 0;
    virtual void create_circuit(Index r2) = 0;
    virtual void transfer() = 0;

    virtual void project_cone(IndexSet& zero_supp, std::vector<Index>& con_map, IndexSet& zeros) = 0;
    virtual bool is_two_dimensional_face(const std::vector<Index>& con_map, const IndexSet& diff) = 0;

    IndexSet r1_supp;
};

template <class T, class IndexSet>
class RaySubState : public RaySubStateAPI<IndexSet>
{
public:
    RaySubState(RayState<T,IndexSet>& _state);
    virtual ~RaySubState();

    void set_r1_index(Index r1);
    void create_ray(Index r2);
    void create_circuit(Index r2);
    void transfer();

    void project_cone(IndexSet& zero_supp, std::vector<Index>& con_map, IndexSet& zeros);
    bool is_two_dimensional_face(const std::vector<Index>& con_map, const IndexSet& diff);
protected:
    typedef RaySubStateAPI<IndexSet> Base;

    RayState<T,IndexSet>& state;
    VectorArrayT<T> new_rays;
    std::vector<IndexSet> new_supps;

    // Temporary variables.
    ConeC<T> sub_cone;
    VectorT<T> temp;
    Index i1;
    T s1;
};

} // namespace _4ti2_

// Definitions of template class functions.
#include "qsolve/RayState.hpp"

#endif
