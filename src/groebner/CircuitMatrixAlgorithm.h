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

#ifndef _4ti2_groebner__CircuitMatrixAlgorithm_
#define _4ti2_groebner__CircuitMatrixAlgorithm_

#include "groebner/CircuitImplementation.h"
#include "groebner/VectorArray.h"
#include "groebner/Timer.h"

namespace _4ti2_
{

template <class IndexSet>
class CircuitMatrixAlgorithm : public CircuitImplementation<IndexSet>
{
public:
    CircuitMatrixAlgorithm();
    virtual ~CircuitMatrixAlgorithm();

    virtual void compute(
                    VectorArray& matrix,
                    VectorArray& vs,
                    VectorArray& circuits,
                    const IndexSet& rs,
                    const IndexSet& cirs);

protected:
    void compute0(  VectorArray& matrix,
                    VectorArray& vs,
                    VectorArray& circuits,
                    const IndexSet& rs,
                    const IndexSet& cirs);
    void compute1(  VectorArray& matrix,
                    VectorArray& vs,
                    VectorArray& circuits,
                    const IndexSet& rs,
                    const IndexSet& cirs);

    bool rank_check(VectorArray& matrix,
                    VectorArray& temp_matrix,
                    IndexSet& temp_diff,
                    int r1_rows);

    void zero_cols( VectorArray& matrix,
                    IndexSet& r1_supp,
                    IndexSet& temp_zero_cols,
                    int r1_rows);

    void compute(
                VectorArray& matrix,
                VectorArray& vs,
                int next_col,
                int num_remaining,
                IndexSet remaining,
                int remaining_row,
                int r1_start, int r1_end,
                int r2_start, int r2_end,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps);

    void create(VectorArray& vs,
                int next_col,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps,
                int r1, int r2,
                Vector& temp, IndexSet& temp_supp);

    Timer t;
};

} // namespace _4ti2_

// Include the template definitions.
#include "groebner/CircuitMatrixAlgorithm.tpp"

#endif
