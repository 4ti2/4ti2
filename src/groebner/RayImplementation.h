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

#ifndef _4ti2_groebner__RayImplementation_
#define _4ti2_groebner__RayImplementation_

#include "groebner/VectorArray.h"
#include "groebner/QSolveConsOrder.h"

namespace _4ti2_
{

typedef std::vector<int> Fathers;
typedef std::vector<int> History;

template <class IndexSet>
class RayImplementation
{
public:
    RayImplementation();
    RayImplementation(QSolveConsOrder o);
    virtual ~RayImplementation();

    void set_constraint_order(QSolveConsOrder o);

    virtual IndexSet compute(
                    const VectorArray& matrix,
                    VectorArray& vs,
                    const IndexSet& rs) = 0;
    virtual IndexSet compute(
                    const VectorArray& matrix,
                    VectorArray& vs,
                    std::vector<IndexSet>& supports,
                    const IndexSet& rs) = 0;

protected:
    void column_count(
                    const VectorArray& vs,
                    Index c,
                    int& positive_count,
                    int& negative_count,
                    int& zero_count);

    int next_column(const VectorArray& vs,
                    const IndexSet& remaining,
                    int& positive_count,
                    int& negative_count,
                    int& zero_count);

    bool (*compare)(int next_positive_count,
                    int next_negative_count,
                    int next_zero_count,
                    int positive_count,
                    int negative_count,
                    int zero_count);
    static bool maxinter(
                    int next_positive_count,
                    int next_negative_count,
                    int next_zero_count,
                    int positive_count,
                    int negative_count,
                    int zero_count);
    static bool maxcutoff(
                    int next_positive_count,
                    int next_negative_count,
                    int next_zero_count,
                    int positive_count,
                    int negative_count,
                    int zero_count);
    static bool mincutoff(
                    int next_positive_count,
                    int next_negative_count,
                    int next_zero_count,
                    int positive_count,
                    int negative_count,
                    int zero_count);
    static bool minindex(
                    int next_positive_count,
                    int next_negative_count,
                    int next_zero_count,
                    int positive_count,
                    int negative_count,
                    int zero_count);

    void sort(      VectorArray& vs,
                    std::vector<IndexSet>& supports,
                    int next_col,
                    int next_zero_count,
                    int next_positive_count, 
                    int next_negative_count);
    void sort(      VectorArray& vs,
                    std::vector<IndexSet>& supports,
                    Fathers& fathers, 
                    std::vector<IndexSet>& zeros, 
                    int next_col,
                    int next_zero_count,
                    int next_positive_count, 
                    int next_negative_count);

    void create_new_vector(
                    VectorArray& vs,
                    std::vector<IndexSet>& supports,
                    int r1, int r2, int next_col,
                    int next_positive_count, int next_negative_count,
                    Vector& temp, IndexSet& temp_supp);
    void create_new_vector(
                    VectorArray& vs,
                    std::vector<IndexSet>& supports,
                    Fathers& fathers,
                    int r1, int r2, int next_col,
                    int next_positive_count, int next_negative_count,
                    Vector& temp, IndexSet& temp_supp);

    QSolveConsOrder order;
};

} // namespace _4ti2_

// Definitions of template class functions.
#include "groebner/RayImplementation.tpp"

#endif
