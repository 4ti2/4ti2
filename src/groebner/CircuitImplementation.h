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

#ifndef _4ti2__CircuitImplementation_
#define _4ti2__CircuitImplementation_

#include "VectorArray.h"
#include "Timer.h"

namespace _4ti2_
{

template <class IndexSet>
class CircuitImplementation
{
public:
    CircuitImplementation() {}
    virtual ~CircuitImplementation() {}

    virtual void compute(
                    VectorArray& matrix,
                    VectorArray& vs,
                    const IndexSet& rs,
                    const IndexSet& cirs);

protected:
    void column_count(
                    const VectorArray& vs,
                    Index c,
                    int& pos_count,
                    int& neg_count,
                    int& zero_count);

    int next_column(const VectorArray& vs, const IndexSet& remaining);

    void sort_nonzeros(
                VectorArray& vs,
                int start, int end,
                std::vector<bool>& rays,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps,
                int next_col,
                int& middle);
    void sort_positives(
                VectorArray& vs,
                int start, int end,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps,
                int next_col,
                int& middle);
    void sort_negatives(
                VectorArray& vs,
                int start, int end,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps,
                int next_col,
                int& middle);
    void sort_rays(
                VectorArray& vs,
                int start, int end,
                std::vector<bool>& rays,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps,
                int& middle);

    void split_rays(
                VectorArray& vs,
                std::vector<bool>& rays,
                VectorArray& circuits);

    void update_supports(
                    std::vector<IndexSet>& supps,
                    int index, int start, int end);

    void switch_supports(
                int start, int end,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps);

    void check(
                const Vector& v,
                const IndexSet& supp,
                const IndexSet& pos_supp,
                const IndexSet& neg_supp,
                const IndexSet& remaining);
};

} // namespace _4ti2_

// Include template definitions.
#include "CircuitImplementation.tpp"

#endif
