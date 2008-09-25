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

#include "groebner/CircuitSupportAlgorithm.h"
#include "groebner/RaySupportAlgorithm.h"
#include "groebner/DiagonalAlgorithm.h"
#include "groebner/HermiteAlgorithm.h"
#include "groebner/Euclidean.h"
#include "groebner/SupportTree.h"
#include "groebner/OnesTree.h"
#include "groebner/Globals.h"

#include "groebner/VectorArrayStream.h"
#include "groebner/VectorStream.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/BitSetStream.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

//#define DEBUG_4ti2(X) X
#include "groebner/Debug.h"

#define TREE SupportTree

using namespace _4ti2_;

template <class IndexSet>
CircuitSupportAlgorithm<IndexSet>::CircuitSupportAlgorithm()
{
}

template <class IndexSet>
CircuitSupportAlgorithm<IndexSet>::~CircuitSupportAlgorithm()
{
}

template <class IndexSet>
void
CircuitSupportAlgorithm<IndexSet>::compute(
                const VectorArray& matrix,
                VectorArray& vs,
                VectorArray& circuits,
                const IndexSet& rs,
                const IndexSet& cirs)
{
    assert(IndexSet::set_disjoint(rs, cirs));
    compute1(matrix, vs, circuits, rs, cirs);
}

template <class IndexSet>
void
CircuitSupportAlgorithm<IndexSet>::compute1(
                const VectorArray& matrix,
                VectorArray& vs,
                VectorArray& circuits,
                const IndexSet& rs,
                const IndexSet& cirs)
{
    assert(matrix.get_size() == vs.get_size());
    assert(cirs.get_size() == vs.get_size());

    t.reset();
    Index num_cols = vs.get_size();
    Index full_num_cols = vs.get_size()+cirs.count();
    IndexSet urs(rs); // The variables that are unrestricted in sign.
    urs.set_union(cirs);
    urs.set_complement();

    // We duplicate the circuit components for the support sets.
    // The col_map maps from a circuit component to its duplicate.
    std::vector<int> col_map(num_cols, -1);
    Index index = num_cols;
    DEBUG_4ti2(*out << "Column Map:\n";)
    for (int i = 0; i < num_cols; ++i)
    {
        if (cirs[i]) { col_map[i] = index; ++index; }
        DEBUG_4ti2(*out << " " << col_map[i];)
    }
    DEBUG_4ti2(*out << "\n";)
    assert(index == full_num_cols);

    DEBUG_4ti2(*out << "Dimension = " << vs.get_number() << "\n";)

    Index ray_rows = diagonal(vs, rs); // Compute ray diagonal normal form.
    Index cir_rows = diagonal(vs, cirs, ray_rows); // Compute circuit diagonal normal form.
    vs.remove(cir_rows, vs.get_number()); // Remove unwanted rows.

    int codim = vs.get_size() - vs.get_number();
    DEBUG_4ti2(*out << "Codimension = " << codim << "\n";)

    circuits.renumber(0);
    VectorArray::transfer(vs, ray_rows, vs.get_number(), circuits, 0);

    VectorArray ray_matrix(matrix);
    // We find the entries on the diagonal for the circuits.
    // Add unit vectors to the ray matrix for the circuit diagonals.
    IndexSet diagonals(num_cols);
    Index col = 0;
    for (Index r = 0; r < circuits.get_number(); ++r)
    {
        while (circuits[r][col] == 0) { ++col; }
        diagonals.set(col);
        Vector unit(num_cols, 0);
        unit[col] = 1;
        ray_matrix.insert(unit);
        ++col;
    }

    DEBUG_4ti2(*out << "Ray dimension is " << vs.get_number() << "\n";)
    DEBUG_4ti2(*out << "Circuit dimension is " << circuits.get_number() << "\n";)

    // Compute the rays.
    std::vector<IndexSet> supps;
    RaySupportAlgorithm<IndexSet> ray_algorithm;
    ray_algorithm.compute(ray_matrix, vs, supps, rs);
    if (cirs.empty()) { return; }
    *out << "Circuit Support Algorithm.\n";
    std::vector<IndexSet> pos_supps(supps.size(), IndexSet(full_num_cols, false));
    std::vector<IndexSet> neg_supps(supps.size(), IndexSet(full_num_cols, false));
    std::vector<bool> rays(supps.size(), true);
    for (int i = 0; i < (int) supps.size(); ++i)
    {
        for (int j = 0; j < num_cols; ++j)
        {
            if (supps[i][j])
            {
                pos_supps[i].set(j);
                if (cirs[j]) { neg_supps[i].set(col_map[j]); }
            }
        }
    }

    // Add the circuits.
    col = 0;
    for (Index r = 0; r < circuits.get_number(); ++r)
    {
        while (!diagonals[col]) { ++col; }

        IndexSet supp(num_cols, false);
        supp.set(col);
        supps.push_back(supp);

        IndexSet neg_supp(full_num_cols, false);
        neg_supp.set(col_map[col]);
        neg_supps.push_back(neg_supp);

        IndexSet pos_supp(full_num_cols, false);
        pos_supp.set(col);
        pos_supps.push_back(pos_supp);

        rays.push_back(false);
        ++col;
    }
    VectorArray::transfer(circuits, vs);

    // The remaining columns to process.
    IndexSet remaining(cirs);
    remaining.set_difference(diagonals);
    int num_remaining = remaining.count();

    // The columns with relaxed non-negativity constraints.
    IndexSet relaxed(remaining);
    relaxed.set_union(urs);
    int num_relaxed = relaxed.count();

    while (vs.get_number() > 0 && num_remaining > 0)
    {
        DEBUG_4ti2(
            for (int i = 0; i < vs.get_number(); ++i)
            {
                *out << "VS[" << i << "] " << rays[i] << "\n";
                *out << vs[i] << "\n" << supps[i] << "\n";
                *out << pos_supps[i] << "\n" << neg_supps[i] << "\n";
                check(vs[i], supps[i], pos_supps[i], neg_supps[i], remaining);
            }
        )
 
        // Find the next column.
        Index next_col = next_column(vs, remaining);

        int start = 0; int end = vs.get_number(); int middle;
        // We sort the vectors into nonzeros and then zeros.
        sort_nonzeros(vs, start, end, rays, supps, pos_supps, neg_supps, next_col, middle);
        int nonzero_start = start, nonzero_end = middle;
        //int zero_start = middle, zero_end = end;
        // We sort the nonzeros into rays and circuits.
        sort_rays(vs, nonzero_start, nonzero_end, rays, supps, pos_supps, neg_supps, middle);
        int ray_start = nonzero_start, ray_end = middle;
        int cir_start = middle, cir_end = nonzero_end;
        // We sort the rays into positives and then negatives.
        sort_positives(vs, ray_start, ray_end, supps, pos_supps, neg_supps, next_col, middle);
        int pos_ray_start = ray_start, pos_ray_end = middle;
        int neg_ray_start = middle, neg_ray_end = ray_end;
        // We sort the circuits into positives and the negatives.
        sort_positives(vs, cir_start, cir_end, supps, pos_supps, neg_supps, next_col, middle);
        int pos_cir_start = cir_start, pos_cir_end = middle;
        int neg_cir_start = middle, neg_cir_end = cir_end;

        // Ouput statistics.
        *out << "\r";
        *out << "  Left = " << std::setw(3) << num_remaining;
        *out << ",  Col = " << std::setw(3) << next_col;
        *out << ",  Size = " << std::setw(8) << vs.get_number();
        DEBUG_4ti2(
            *out << "Rays (+,-) = (" << pos_ray_end-pos_ray_start;
            *out << "," << neg_ray_end-neg_ray_start << ") ";
            *out << "Circuits (+,-) = (" << pos_cir_end-pos_cir_start;
            *out << "," << neg_cir_end-neg_cir_start << ")\n";
        )

        // Note that the tree needs the ordering of the current vectors to be
        // constant.
        DEBUG_4ti2(*out << "Building Tree ... " << std::flush;)
        TREE<IndexSet> tree(pos_supps, pos_supps.size());
        for (int i = 0; i < vs.get_number(); ++i)
        {
            if (!rays[i]) { tree.insert(neg_supps[i], i); }
        }
        DEBUG_4ti2(*out << "done." << std::endl;)
        DEBUG_4ti2(tree.dump();)

        // Switch the positive and negative supports, so that it is as if all
        // vectors have a positive entry in the next column.
        switch_supports(neg_ray_start, neg_ray_end, pos_supps, neg_supps);
        switch_supports(neg_cir_start, neg_cir_end, pos_supps, neg_supps);

        //DEBUG_4ti2(*out << "Remaining row " << remaining_row << "\n";)
        int previous_size = vs.get_number();
        // Positive ray combinations.
        DEBUG_4ti2(*out << "+r\n";)
        compute(tree, vs, next_col, full_num_cols, num_remaining, num_relaxed, codim,
                        pos_ray_start, pos_ray_end, neg_ray_start, cir_end,
                        supps, pos_supps, neg_supps);
        // Negative ray combinations;
        DEBUG_4ti2(*out << "-r\n";)
        compute(tree, vs, next_col, full_num_cols, num_remaining, num_relaxed, codim,
                        neg_ray_start, neg_ray_end, cir_start, cir_end,
                        supps, pos_supps, neg_supps);
        rays.insert(rays.end(), vs.get_number()-previous_size, true);
        previous_size = vs.get_number();
        // Circuit combinations.
        DEBUG_4ti2(*out << "c\n";)
        compute(tree, vs, next_col, full_num_cols, num_remaining, num_relaxed, codim,
                        cir_start, cir_end, cir_start, cir_end,
                        supps, pos_supps, neg_supps);
        rays.insert(rays.end(), vs.get_number()-previous_size, false);

        // Switch back the positive and negative supports.
        switch_supports(neg_ray_start, neg_ray_end, pos_supps, neg_supps);
        switch_supports(neg_cir_start, neg_cir_end, pos_supps, neg_supps);

        // Update the supp vectors for the next_col.
        update_supports(supps, next_col, nonzero_start, nonzero_end);
        update_supports(pos_supps, next_col, pos_ray_start, pos_ray_end);
        update_supports(pos_supps, next_col, pos_cir_start, pos_cir_end);
        update_supports(pos_supps, col_map[next_col], neg_ray_start, neg_ray_end);
        update_supports(pos_supps, col_map[next_col], neg_cir_start, neg_cir_end);
        update_supports(neg_supps, next_col, neg_ray_start, neg_ray_end);
        update_supports(neg_supps, next_col, neg_cir_start, neg_cir_end);
        update_supports(neg_supps, col_map[next_col], pos_ray_start, pos_ray_end);
        update_supports(neg_supps, col_map[next_col], pos_cir_start, pos_cir_end);

        *out << "\r";
        *out << "  Left = " << std::setw(3) << num_remaining;
        *out << ",  Col = " << std::setw(3) << next_col;
        *out << ",  Size = " << std::setw(8) << vs.get_number() << ",  Time: " << t;
        *out << "                \n";

        remaining.unset(next_col);
        --num_remaining;
        relaxed.unset(next_col);
        --num_relaxed;
    }

    CircuitImplementation<IndexSet>::split_rays(vs, rays, circuits);

    Vector zero(vs.get_size(), 0);
    for (int i = 0; i < circuits.get_number(); ++i)
    {
        if (circuits[i] <= zero) { vs[i].mul(-1); }
    }
}

#if 1
template <class IndexSet>
void
CircuitSupportAlgorithm<IndexSet>::compute(
                TREE<IndexSet>& tree,
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
                std::vector<IndexSet>& neg_supps)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }
    DEBUG_4ti2(*out << "R1 [" << r1_start << "..." << r1_end << "]\n";)
    DEBUG_4ti2(*out << "R2 [" << r2_start << "..." << r2_end << "]\n";)
    int num_cols = vs.get_size();

    char buffer[256];
    sprintf(buffer, "  Left = %3d,  Col = %3d,", num_remaining, next_col);

    IndexSet temp_supp(num_cols);
    IndexSet full_temp_supp(full_num_cols);
    IndexSet r1_supp(num_cols);
    IndexSet r1_pos_supp(full_num_cols);
    IndexSet r1_neg_supp(full_num_cols);

    Vector temp(num_cols);

    int index_count = 0;
    for (int r1 = r1_start; r1 < r1_end; ++r1)
    {
        r1_supp = supps[r1];
        r1_pos_supp = pos_supps[r1];
        r1_neg_supp = neg_supps[r1];
        if (r2_start == r1) { ++r2_start; }

        if (!r1_supp.less_than_equal(codim-num_relaxed))
        {
            for (Index r2 = r2_start; r2 < r2_end; ++r2)
            {
                if (!IndexSet::set_disjoint(r1_pos_supp, pos_supps[r2])) { continue; }
                IndexSet::set_difference(supps[r2], r1_supp, temp_supp);
                if (temp_supp.power_of_2())
                {
                    create(vs, next_col, supps, pos_supps, neg_supps,
                                    r1, r2, temp, temp_supp, full_temp_supp);
                }
            }
        }
        else
        {
            int r1_count = r1_supp.count();
            for (Index r2 = r2_start; r2 < r2_end; ++r2)
            {
                if (!IndexSet::set_disjoint(r1_pos_supp, pos_supps[r2])) { continue; }
                IndexSet::set_difference(supps[r2], r1_supp, temp_supp);
                if (!temp_supp.less_than_equal(codim-num_relaxed-r1_count+2)) { continue; }
                IndexSet::set_union(r1_pos_supp,neg_supps[r2],full_temp_supp);
                if (!tree.dominated(full_temp_supp, r1, r2))
                {
                    create(vs, next_col, supps, pos_supps, neg_supps,
                                    r1, r2, temp, temp_supp, full_temp_supp);
                }
            }
        }
        if (index_count % Globals::output_freq == 0)
        {
           *out << "\r" << buffer;
           *out << "  Size = " << std::setw(8) << vs.get_number();
           *out << ",  Index = " << r1 << "/" << r2_end << std::flush;
        }
        ++index_count;
    }
    *out << "\r" << buffer;
    *out << "  Size = " << std::setw(8) << vs.get_number();
    *out << ",  Index = " << r1_end << "/" << r2_end << std::flush;
}
#else
template <class IndexSet>
void
CircuitSupportAlgorithm<IndexSet>::compute(
                TREE<IndexSet>& tree,
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
                std::vector<IndexSet>& neg_supps)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }
    DEBUG_4ti2(*out << "R1 [" << r1_start << "..." << r1_end << "]\n";)
    DEBUG_4ti2(*out << "R2 [" << r2_start << "..." << r2_end << "]\n";)
    int num_cols = vs.get_size();

    char buffer[256];
    sprintf(buffer, "  Left = %3d  Col = %3d", num_remaining, next_col);

    IndexSet temp_supp(num_cols);
    IndexSet zero_supp(full_num_cols);
    IndexSet full_temp_supp(full_num_cols);
    IndexSet r1_supp(num_cols);
    IndexSet r1_pos_supp(full_num_cols);
    IndexSet r1_neg_supp(full_num_cols);

    Vector temp(num_cols);
    std::vector<int> indices;

    int index_count = 0;
    for (int r1 = r1_start; r1 < r1_end; ++r1)
    {
        r1_supp = supps[r1];
        r1_pos_supp = pos_supps[r1];
        r1_neg_supp = neg_supps[r1];
        if (r2_start == r1) { ++r2_start; }

        //if (r1_supp.count() == codim-num_relaxed+1)
        if (!r1_supp.less_than_equal(codim-num_relaxed))
        {
            for (Index r2 = r2_start; r2 < r2_end; ++r2)
            {
                if (!IndexSet::set_disjoint(r1_pos_supp, pos_supps[r2])) { continue; }
                IndexSet::set_difference(supps[r2], r1_supp, temp_supp);
                if (temp_supp.power_of_2())
                {
                    create(vs, next_col, supps, pos_supps, neg_supps,
                                    r1, r2, temp, temp_supp, full_temp_supp);
                }
            }
        }
        else
        {
            //*out << "R1 " << r1;
            indices.clear();
            zero_supp.zero();
            IndexSet::set_complement(r1_pos_supp, full_temp_supp);
            tree.find_diff(indices, full_temp_supp, 1);
            std::sort(indices.begin(), indices.end());
            indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
            for (unsigned int i = 0; i < indices.size(); ++i)
            {
                int r2 = indices[i];
                if (IndexSet::set_disjoint(r1_pos_supp, neg_supps[r2]))
                {
                    zero_supp.set_union(pos_supps[r2]);
                    if (r2 >= r2_start && r2 < r2_end && is_r2_positive != is_r1_positive)
                    {
                        create(vs, next_col, supps, pos_supps, neg_supps,
                                    r1, r2, temp, temp_supp, full_temp_supp);
                    }
                }
                if (IndexSet::set_disjoint(r1_pos_supp, pos_supps[r2]))
                {
                    zero_supp.set_union(neg_supps[r2]);
                    if (r2 >= r2_start && r2 < r2_end && is_r2_positive == is_r1_positive)
                    {
                        create(vs, next_col, supps, pos_supps, neg_supps,
                                    r1, r2, temp, temp_supp, full_temp_supp);
                    }
                }
            }
            zero_supp.set_difference(r1_pos_supp);
            zero_supp.set_union(r1_neg_supp);
            //*out << "\nR1 " << r1 << "\n";
            //*out << vs[r1] << "\n";
            //*out << supps[r1] << "\n";
            //*out << "Zero Supp:\n" << zero_supp << "\n";
            int r1_count = r1_supp.count();
            for (Index r2 = r2_start; r2 < r2_end; ++r2)
            {
                if (!IndexSet::set_disjoint(zero_supp, neg_supps[r2])) { continue; }
                //if (!IndexSet::set_disjoint(r1_pos_supp, pos_supps[r2])) { continue; }
                IndexSet::set_difference(supps[r2], r1_supp, temp_supp);
                if (!temp_supp.less_than_equal(codim-num_relaxed-r1_count+2)) { continue; }
                IndexSet::set_union(r1_pos_supp,neg_supps[r2],full_temp_supp);
                if (!tree.dominated(full_temp_supp, r1, r2))
                {
                    create(vs, next_col, supps, pos_supps, neg_supps,
                                    r1, r2, temp, temp_supp, full_temp_supp);
                }
            }
        }
        if (index_count % Globals::output_freq == 0)
        {
           *out << "\r" << buffer;
           *out << "  Size = " << std::setw(8) << vs.get_number();
           *out << ",  Index = " << r1 << "/" << r2_end << std::flush;
        }
        ++index_count;
    }
    *out << "\r" << buffer;
    *out << "  Size = " << std::setw(8) << vs.get_number();
    *out << ",  Index = " << r1_end << "/" << r2_end << std::flush;
}
#endif

template <class IndexSet>
void
CircuitSupportAlgorithm<IndexSet>::create(
                VectorArray& vs,
                int next_col,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps,
                int r1, int r2,
                Vector& temp, IndexSet& temp_supp, IndexSet& full_temp_supp)
{
    if (vs[r2][next_col] > 0)
    {
        Vector::sub(vs[r1],vs[r2][next_col],vs[r2],vs[r1][next_col],temp);
    }
    else
    {
        Vector::sub(vs[r2],vs[r1][next_col],vs[r1],vs[r2][next_col],temp);
    }
    temp.normalise();
    vs.insert(temp);
    IndexSet::set_union(supps[r1],supps[r2],temp_supp);
    supps.push_back(temp_supp);

    if (vs[r1][next_col] > 0)
    {
        IndexSet::set_union(pos_supps[r1], neg_supps[r2], full_temp_supp);
        pos_supps.push_back(full_temp_supp);
        IndexSet::set_union(neg_supps[r1], pos_supps[r2], full_temp_supp);
        neg_supps.push_back(full_temp_supp);
    }
    else
    {
        IndexSet::set_union(neg_supps[r1], pos_supps[r2], full_temp_supp);
        pos_supps.push_back(full_temp_supp);
        IndexSet::set_union(pos_supps[r1], neg_supps[r2], full_temp_supp);
        neg_supps.push_back(full_temp_supp);
    }

    DEBUG_4ti2(
        *out << "\nNEW VECTOR:\n";
        *out << "R1 " << r1 << "\n";
        *out << vs[r1] << "\n";
        *out << supps[r1] << "\n";
        *out << pos_supps[r1] << "\n";
        *out << neg_supps[r1] << "\n";
        *out << "R2 " << r2 << "\n";
        *out << vs[r2] << "\n";
        *out << supps[r2] << "\n";
        *out << pos_supps[r2] << "\n";
        *out << neg_supps[r2] << "\n";
        int i = vs.get_number()-1;
        *out << "V:\n" << vs[i] << "\n";
        *out << supps[i] << "\n";
        *out << pos_supps[i] << "\n";
        *out << neg_supps[i] << "\n";
    )
}
