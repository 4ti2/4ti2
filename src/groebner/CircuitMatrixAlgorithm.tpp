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

#include "CircuitMatrixAlgorithm.h"
#include "RayMatrixAlgorithm.h"
#include "DiagonalAlgorithm.h"
#include "HermiteAlgorithm.h"
#include "Debug.h"
#include "Globals.h"
#include "Timer.h"

#include "VectorArrayStream.h"
#include "VectorStream.h"
#include "VectorArrayStream.h"
#include <iostream>
#include <iomanip>

using namespace _4ti2_;

template <class IndexSet>
CircuitMatrixAlgorithm<IndexSet>::CircuitMatrixAlgorithm()
{
}

template <class IndexSet>
CircuitMatrixAlgorithm<IndexSet>::~CircuitMatrixAlgorithm()
{
}

template <class IndexSet>
void
CircuitMatrixAlgorithm<IndexSet>::compute(
                VectorArray& matrix,
                VectorArray& vs,
                VectorArray& circuits,
                const IndexSet& rs,
                const IndexSet& cirs)
{
    assert(IndexSet::set_disjoint(rs, cirs));
    compute1(matrix, vs, circuits, rs, cirs);
}

// Assumes cone is pointed.
template <class IndexSet>
void
CircuitMatrixAlgorithm<IndexSet>::compute1(
                VectorArray& matrix,
                VectorArray& vs,
                VectorArray& circuits,
                const IndexSet& rs,
                const IndexSet& cirs)
{
    assert(matrix.get_size() == vs.get_size());
    assert(cirs.get_size() == vs.get_size());

    t.reset();
    Index num_cols = vs.get_size();

    DEBUG_4ti2(*out << "Dimension = " << vs.get_number() << "\n";)

    int codim = upper_triangle(matrix);
    matrix.remove(codim, matrix.get_number());
    DEBUG_4ti2(*out << "Codimension = " << codim << "\n";)
    VectorArray orig_matrix(matrix);

    Index ray_rows = diagonal(vs, rs); // Compute ray diagonal normal form.
    Index cir_rows = diagonal(vs, cirs, ray_rows); // Compute circuit diagonal normal form.
    vs.remove(cir_rows, vs.get_number()); // Remove unwanted rows.

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
    RayMatrixAlgorithm<IndexSet> ray_algorithm;
    ray_algorithm.compute(ray_matrix, vs, supps, rs);

    // If there are no circuit components, then there is nothing to do.
    if (cirs.empty()) { return; }
    *out << "Circuit Matrix Algorithm.\n";

    std::vector<IndexSet> pos_supps(supps);
    std::vector<IndexSet> neg_supps(supps.size(), IndexSet(num_cols, false));
    std::vector<bool> rays(supps.size(), true);

    // Add the circuits.
    col = 0;
    for (Index r = 0; r < circuits.get_number(); ++r)
    {
        while (!diagonals[col]) { ++col; }
        IndexSet supp(num_cols, false);
        neg_supps.push_back(supp);
        supp.set(col);
        supps.push_back(supp);
        pos_supps.push_back(supp);
        rays.push_back(false);
        ++col;
    }
    VectorArray::transfer(circuits, vs);

    // The remaining columns to process.
    IndexSet remaining(cirs);
    remaining.set_difference(diagonals);
    int num_remaining = remaining.count();

    while (vs.get_number() > 0 && num_remaining > 0)
    {
        DEBUG_4ti2(
            for (int i = 0; i < vs.get_number(); ++i)
            {
                *out << "VS[" << i << "]" << "\n";
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
        // We sort the circuits into positives and then negatives.
        sort_positives(vs, cir_start, cir_end, supps, pos_supps, neg_supps, next_col, middle);
        int pos_cir_start = cir_start, pos_cir_end = middle;
        int neg_cir_start = middle, neg_cir_end = cir_end;

        *out << "\r";
        *out << "  Left = " << std::setw(3) << num_remaining;
        *out << "  Col = " << std::setw(3) << next_col;
        *out << "  Size = " << std::setw(8) << vs.get_number();
        DEBUG_4ti2(
            *out << "Rays (+,-) = (" << pos_ray_end-pos_ray_start;
            *out << "," << neg_ray_end-neg_ray_start << ") ";
            *out << "Circuits (+,-) = (" << pos_cir_end-pos_cir_start;
            *out << "," << neg_cir_end-neg_cir_start << ")\n";
        )

        // Switch the positive and negative supports, so that it is as if all
        // vectors have a positive entry in the next column.
        switch_supports(neg_ray_start, neg_ray_end, pos_supps, neg_supps);
        switch_supports(neg_cir_start, neg_cir_end, pos_supps, neg_supps);

        matrix = orig_matrix;
        int remaining_row = upper_triangle(matrix, remaining, 0);
        DEBUG_4ti2(*out << "Remaining row " << remaining_row << "\n";)
        int previous_size = vs.get_number();
        // Positive ray combinations.
        DEBUG_4ti2(*out << "+r\n";)
        compute(matrix, vs, next_col, num_remaining, remaining, remaining_row,
                        pos_ray_start, pos_ray_end, neg_ray_start, cir_end,
                        supps, pos_supps, neg_supps);
        // Negative ray combinations.
        DEBUG_4ti2(*out << "-r\n";)
        compute(matrix, vs, next_col, num_remaining, remaining, remaining_row,
                        neg_ray_start, neg_ray_end, cir_start, cir_end,
                        supps, pos_supps, neg_supps);
        rays.insert(rays.end(), vs.get_number()-previous_size, true);
        previous_size = vs.get_number();
        // circuit combinations.
        DEBUG_4ti2(*out << "c\n";)
        compute(matrix, vs, next_col, num_remaining, remaining, remaining_row,
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
        update_supports(neg_supps, next_col, neg_ray_start, neg_ray_end);
        update_supports(neg_supps, next_col, neg_cir_start, neg_cir_end);

        *out << "\r";
        *out << "  Left = " << std::setw(3) << num_remaining;
        *out << "  Col = " << std::setw(3) << next_col;
        *out << "  Size = " << std::setw(8) << vs.get_number();
        *out << "  Time: " << t << "                \n";

        remaining.unset(next_col);
        --num_remaining;
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
CircuitMatrixAlgorithm<IndexSet>::compute(
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
                std::vector<IndexSet>& neg_supps)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }

    char buffer[256];
    sprintf(buffer, "  Left = %3d  Col = %3d", num_remaining, next_col);

    DEBUG_4ti2(*out << "R1 [" << r1_start << "..." << r1_end << "]\n";)
    DEBUG_4ti2(*out << "R2 [" << r2_start << "..." << r2_end << "]\n";)

    int num_cols = matrix.get_size();
    int codim = matrix.get_number();

    IndexSet temp_supp(num_cols);
    IndexSet r1_supp(num_cols);
    IndexSet r1_pos_supp(num_cols);
    IndexSet r1_neg_supp(num_cols);
    IndexSet zero_cols_supp(num_cols);

    Vector temp(num_cols);
    VectorArray orig_matrix(matrix);
    VectorArray temp_matrix(matrix.get_number(), matrix.get_size(), 0);

    int index_count = 0;
    for (int r1 = r1_start; r1 < r1_end; ++r1)
    {
        r1_supp = supps[r1];
        r1_pos_supp = pos_supps[r1];
        r1_neg_supp = neg_supps[r1];
        if (r2_start == r1) { ++r2_start; }

        //if (r1_supp.count() == codim-num_remaining+1)
        if (!r1_supp.less_than_equal(codim-num_remaining))
        {
            for (Index r2 = r2_start; r2 < r2_end; ++r2)
            {
                IndexSet::set_difference(supps[r2], r1_supp, temp_supp);
                if (!temp_supp.power_of_2()) { continue; }
                if (IndexSet::set_disjoint(r1_pos_supp, pos_supps[r2]) &&
                    IndexSet::set_disjoint(r1_neg_supp, neg_supps[r2]))
                {
                    create(vs, next_col, supps, pos_supps, neg_supps,
                                    r1, r2, temp, temp_supp);
                }
            }
        }
        else
        {
            // TODO: Avoid unnecessary copying of rows.
            matrix = orig_matrix;
            int r1_rows = upper_triangle(matrix, r1_supp, remaining_row);
            // Find the columns in the matrix which are zero.
            zero_cols(matrix, r1_supp, zero_cols_supp, r1_rows);
            // Positive and Positive combinations.
            for (Index r2 = r2_start; r2 < r2_end; ++r2)
            {
                IndexSet::set_intersection(supps[r2], zero_cols_supp, temp_supp);
                //if (temp_supp.count() <= 1)
                if (temp_supp.power_of_2())
                {
                    IndexSet::set_difference(supps[r2], r1_supp, temp_supp);
                    //if (temp_supp.count() > codim-r1_rows+1)
                    if (temp_supp.less_than_equal(codim-r1_rows+1))
                    {
                        if (IndexSet::set_disjoint(r1_pos_supp, pos_supps[r2]) &&
                            IndexSet::set_disjoint(r1_neg_supp, neg_supps[r2]))
                        {
                            if (rank_check(matrix, temp_matrix, temp_supp, r1_rows))
                            {
                                create(vs, next_col, supps, pos_supps, neg_supps,
                                    r1, r2, temp, temp_supp);
                            }
                        }
                    }
                }
            }
        }
        if (index_count % Globals::output_freq == 0)
        {
           *out << "\r" << buffer;
           *out << "  Size = " << std::setw(8) << vs.get_number();
           *out << "  Index = " << r1 << "/" << r2_end << std::flush;
        }
        ++index_count;
    }
    *out << "\r" << buffer;
    *out << "  Size = " << std::setw(8) << vs.get_number();
    *out << "  Index = " << r1_end << "/" << r2_end << std::flush;
}
#else
template <class IndexSet>
void
CircuitMatrixAlgorithm<IndexSet>::compute(
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
                std::vector<IndexSet>& neg_supps)
{
    if (r1_start == r1_end || r2_start == r2_end) { return; }

    char buffer[256];
    sprintf(buffer, "  Left = %3d  Col = %3d", num_remaining, next_col);

    DEBUG_4ti2(*out << "R1 [" << r1_start << "..." << r1_end << "]\n";)
    DEBUG_4ti2(*out << "R2 [" << r2_start << "..." << r2_end << "]\n";)

    int num_cols = matrix.get_size();
    int codim = matrix.get_number();

    IndexSet temp_supp(num_cols);
    IndexSet r1_supp(num_cols);
    IndexSet r1_pos_supp(num_cols);
    IndexSet r1_neg_supp(num_cols);
    IndexSet zero_cols_supp(num_cols);

    Vector temp(num_cols);
    VectorArray orig_matrix(matrix);
    VectorArray temp_matrix(matrix.get_number(), matrix.get_size(), 0);

    int index_count = 0;
    for (int r1 = r1_start; r1 < r1_end; ++r1)
    {
        r1_supp = supps[r1];
        r1_pos_supp = pos_supps[r1];
        r1_neg_supp = neg_supps[r1];
        if (r2_start == r1) { ++r2_start; }

        //if (r1_supp.count() == codim-num_remaining+1)
        if (!r1_supp.less_than_equal(codim-num_remaining))
        {
            for (Index r2 = r2_start; r2 < r2_end; ++r2)
            {
                IndexSet::set_difference(supps[r2], r1_supp, temp_supp);
                if (!temp_supp.power_of_2()) { continue; }
                if (IndexSet::set_disjoint(r1_pos_supp, pos_supps[r2]) &&
                    IndexSet::set_disjoint(r1_neg_supp, neg_supps[r2]))
                {
                    create(vs, next_col, supps, pos_supps, neg_supps,
                                    r1, r2, temp, temp_supp);
                }
            }
        }
        else
        {
            // TODO: Avoid unnecessary copying of rows.
            matrix = orig_matrix;
            int r1_rows = upper_triangle(matrix, r1_supp, remaining_row);
            // Find the columns in the matrix which are zero.
            zero_cols(matrix, r1_supp, zero_cols_supp, r1_rows);
            // Positive and Positive combinations.
            for (Index r2 = r2_start; r2 < r2_end; ++r2)
            {
                if (IndexSet::set_disjoint(supps[r2], zero_cols_supp))
                {
                    IndexSet::set_difference(supps[r2], r1_supp, temp_supp);
                    //if (temp_supp.count() > codim-r1_rows+1)
                    if (temp_supp.less_than_equal(codim-r1_rows+1))
                    {
                        if (IndexSet::set_disjoint(r1_pos_supp, pos_supps[r2]) &&
                            IndexSet::set_disjoint(r1_neg_supp, neg_supps[r2]))
                        {
                            if (rank_check(matrix, temp_matrix, temp_supp, r1_rows))
                            {
                                create(vs, next_col, supps, pos_supps, neg_supps,
                                    r1, r2, temp, temp_supp);
                            }
                        }
                    }
                }
                else
                {
                    IndexSet::set_difference(supps[r2], r1_supp, temp_supp);
                    if (temp_supp.power_of_2())
                    {
                        if (IndexSet::set_disjoint(r1_pos_supp, pos_supps[r2]) &&
                            IndexSet::set_disjoint(r1_neg_supp, neg_supps[r2]))
                        {
                            create(vs, next_col, supps, pos_supps, neg_supps,
                                    r1, r2, temp, temp_supp);
                        }
                    }
                }
            }
        }
        if (index_count % Globals::output_freq == 0)
        {
           *out << "\r" << buffer;
           *out << "  Size = " << std::setw(8) << vs.get_number();
           *out << "  Index = " << r1 << "/" << r2_end << std::flush;
        }
        ++index_count;
    }
    *out << "\r" << buffer;
    *out << "  Size = " << std::setw(8) << vs.get_number();
    *out << "  Index = " << r1_end << "/" << r2_end << std::flush;
}
#endif


template <class IndexSet>
bool
CircuitMatrixAlgorithm<IndexSet>::rank_check(
                VectorArray& matrix,
                VectorArray& _temp_matrix,
                IndexSet& temp_diff,
                int r1_rows)
{
    int m = matrix.get_number()-r1_rows;
    int n = temp_diff.count();
    int col_index = 0;
    VectorArray temp_matrix(m,n); // TODO: Handle memory management better.
    for (int c = 0; c < matrix.get_size(); ++c)
    {
        if (temp_diff[c])
        {
            for (int r = 0; r < m; ++r)
            {
                temp_matrix[r][col_index] = matrix[r+r1_rows][c];
            }
            ++col_index;
        }
    }
    int rank = upper_triangle(temp_matrix, m, n);
    return (rank == n-1);
}

template <class IndexSet>
void
CircuitMatrixAlgorithm<IndexSet>::zero_cols(
                VectorArray& matrix,
                IndexSet& r1_supp,
                IndexSet& temp_zero_cols,
                int r1_rows)
{
    temp_zero_cols.zero();
    for (int i = 0; i < temp_zero_cols.get_size(); ++i)
    {
        if (!r1_supp[i])
        {
            int r = r1_rows;
            while (r < matrix.get_number() && matrix[r][i] == 0) { ++r; }
            if (r == matrix.get_number()) { temp_zero_cols.set(i); }
        }
    }
}

template <class IndexSet>
void
CircuitMatrixAlgorithm<IndexSet>::create(
                VectorArray& vs,
                int next_col,
                std::vector<IndexSet>& supps,
                std::vector<IndexSet>& pos_supps,
                std::vector<IndexSet>& neg_supps,
                int r1, int r2,
                Vector& temp, IndexSet& temp_supp)
{
    Vector& v1 = vs[r1];
    Vector& v2 = vs[r2];
    if (v2[next_col] > 0)
    {
        Vector::sub(v1,v2[next_col],v2,v1[next_col],temp);
    }
    else
    {
        Vector::sub(v2,v1[next_col],v1,v2[next_col],temp);
    }
    temp.normalise();
    vs.insert(temp);
    IndexSet::set_union(supps[r1],supps[r2],temp_supp);
    supps.push_back(temp_supp);

    if (v1[next_col] > 0)
    {
        IndexSet::set_union(pos_supps[r1], neg_supps[r2], temp_supp);
        pos_supps.push_back(temp_supp);
        IndexSet::set_union(neg_supps[r1], pos_supps[r2], temp_supp);
        neg_supps.push_back(temp_supp);
    }
    else
    {
        IndexSet::set_union(neg_supps[r1], pos_supps[r2], temp_supp);
        pos_supps.push_back(temp_supp);
        IndexSet::set_union(pos_supps[r1], neg_supps[r2], temp_supp);
        neg_supps.push_back(temp_supp);
    }

    DEBUG_4ti2(
        *out << "\nNEW VECTOR:\n";
        *out << "R1 " << r1 << "\n";
        *out << v1 << "\n";
        *out << supps[r1] << "\n";
        *out << pos_supps[r1] << "\n";
        *out << neg_supps[r1] << "\n";
        *out << "R2 " << r2 << "\n";
        *out << v2 << "\n";
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

