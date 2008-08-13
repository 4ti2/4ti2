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

#ifndef _4ti2_groebner__CircuitSupportAlgorithm_
#define _4ti2_groebner__CircuitSupportAlgorithm_

#include "groebner/CircuitImplementation.h"
#include "groebner/VectorArray.h"
#include "groebner/SupportTree.h"

namespace _4ti2_
{

template <class IndexSet>
class CircuitSupportAlgorithm : public CircuitImplementation<IndexSet>
{
public:
    CircuitSupportAlgorithm();
    virtual ~CircuitSupportAlgorithm();

    virtual void compute(
                    VectorArray& matrix,
                    VectorArray& vs,
                    VectorArray& circuits,
                    const IndexSet& rs,
                    const IndexSet& cirs);

protected:
    void compute1(  VectorArray& matrix,
                    VectorArray& vs,
                    VectorArray& circuits,
                    const IndexSet& rs,
                    const IndexSet& cirs);

    void create(VectorArray& vs,
                int next_col,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps,
                int r1, int r2,
                Vector& temp, IndexSet& temp_supp, IndexSet& full_temp_supp);

    void compute(
                SupportTree<IndexSet>& tree,
                VectorArray& vs,
                int next_col,
                int full_num_cols,
                int num_remaining,
                int num_relaxed,
                int codim,
                int r1_start, int r1_end,
                int r2_start, int r2_end,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps);

    Timer t;
};

} // namespace _4ti2_

// Include template definitions.
#include "groebner/CircuitSupportAlgorithm.tpp"

#endif
