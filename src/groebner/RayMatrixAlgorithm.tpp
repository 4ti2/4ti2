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

#include "RayMatrixAlgorithm.h"
#include "DiagonalAlgorithm.h"
#include "HermiteAlgorithm.h"
#include "Euclidean.h"
#include "Globals.h"
#include "Timer.h"

#include "VectorArrayStream.h"
#include "VectorStream.h"
#include "VectorArrayStream.h"
#include "BitSetStream.h"
#include "SupportTree.h"
#include "OnesTree.h"
#include <iostream>
#include <iomanip>

#include "Debug.h"

#define TREE SupportTree

using namespace _4ti2_;

template <class IndexSet>
RayMatrixAlgorithm<IndexSet>::RayMatrixAlgorithm()
{
}

template <class IndexSet>
RayMatrixAlgorithm<IndexSet>::~RayMatrixAlgorithm()
{
}

template <class IndexSet>
IndexSet
RayMatrixAlgorithm<IndexSet>::compute(
                VectorArray& matrix, 
                VectorArray& vs, 
                std::vector<IndexSet>& supps,  
                const IndexSet& rs)
{
    return compute3(matrix, vs, supps, rs);
}

template <class IndexSet>
IndexSet
RayMatrixAlgorithm<IndexSet>::compute(
                VectorArray& matrix, 
                VectorArray& vs, 
                const IndexSet& rs)
{
    std::vector<IndexSet> supports;
    return compute(matrix, vs, supports, rs);
}

// It is assumed that the cone is pointed.
// This is a slightly more complicated version of the ray computation algorithm.
template <class IndexSet>
IndexSet
RayMatrixAlgorithm<IndexSet>::compute1(
                VectorArray& matrix, 
                VectorArray& vs, 
                std::vector<IndexSet>& supports, 
                const IndexSet& rs)
{
    // Sanity Checks.
    assert(matrix.get_size() == vs.get_size());
    assert(rs.get_size() == vs.get_size());

    *out << "Ray Matrix Algorithm.\n";
    Timer t;

    int num_cols = vs.get_size();
    IndexSet urs(rs); // The variables that are unrestricted in sign.
    urs.set_complement();

    DEBUG_4ti2(*out << "RS\n" << rs << "\n";)
    DEBUG_4ti2(*out << "URS\n" << urs << "\n";)
    DEBUG_4ti2(*out << "The dimension is " << vs.get_number() << ".\n";)

    Index rows = diagonal(vs, rs); // Compute diagonal normal form.
    vs.remove(rows, vs.get_number());

    // We find the entries on the diagonal.
    supports.clear();
    Index col = 0;
    IndexSet diagonals(num_cols);
    for (Index r = 0; r < vs.get_number(); ++r)
    {
        while (vs[r][col] == 0 || urs[col]) { ++col; }
        diagonals.set(col);
        IndexSet support(num_cols, false);
        support.set(col);
        supports.push_back(support);
        ++col;
    }

    int codim = upper_triangle(matrix);
    matrix.remove(codim, matrix.get_number());
    DEBUG_4ti2(*out << "The codimension is " << codim << ".\n";)
    VectorArray orig_matrix(matrix);

    // The remaining columns to compute rays for.
    IndexSet remaining(rs);
    remaining.set_difference(diagonals);
    int num_remaining = remaining.count();

    // Temporary variables.
    IndexSet temp_supp(num_cols);
    IndexSet temp_diff(num_cols);
    IndexSet temp_diff2(num_cols);
    IndexSet temp_zero_cols(num_cols);
    IndexSet r1_supp(num_cols);
    VectorArray temp_matrix(matrix.get_number(), matrix.get_size(), 0);
    Vector temp(vs.get_size());

    long long unsigned int num_dominated = 0;
    long long unsigned int num_added = 0;
    long long unsigned int num_high_dim = 0;
    long long unsigned int num_adjacency_checks = 0;
    while (vs.get_number() > 0 && num_remaining > 0)
    {
        // Statistics.
        //long int num_zero_rows = 0;
        //long int num_rank_checks = 0;

        // Find the next column.
        int next_positive_count, next_negative_count, next_zero_count;
        Index next_col = next_column(vs, remaining,
                                        next_positive_count,
                                        next_negative_count,
                                        next_zero_count);

        char buffer[256];
        sprintf(buffer, "  Left = %3d  Col = %3d", num_remaining, next_col);
        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << "  Time: " << t;
        DEBUG_4ti2(
            *out << "(+,0,-) = (" << next_positive_count << ",";
            *out << next_zero_count << "," << next_negative_count << ")\n";
        )

        // We sort the vectors into zeros, positives, then negatives.
        sort(vs, supports, next_col, next_zero_count, next_positive_count,
                        next_negative_count);

        matrix = orig_matrix;
        temp_supp = remaining;
        temp_supp.set_union(urs);
        int remaining_row = upper_triangle(matrix, temp_supp, 0);
        VectorArray test_matrix(matrix);

        int original_size = vs.get_number();
        int positive_start = next_zero_count;
        int negative_start = next_zero_count+next_positive_count;

        // We wish to reduce the number of matrices we triangularise so we
        // choose the smaller out of the positives and negatives.
        int r1_start;
        int r1_finish;
        int r2_start;
        int r2_finish;
        int index_max;
        if (next_positive_count <= next_negative_count)
        {
            //*out << "Using positive vectors.\n";
            r1_start = positive_start;
            r1_finish = negative_start;
            r2_start = negative_start;
            r2_finish = original_size;
            index_max = next_positive_count;
        }
        else
        {
            //*out << "Using negative vectors.\n";
            r1_start = negative_start;
            r1_finish = original_size;
            r2_start = positive_start;
            r2_finish = negative_start;
            index_max = next_negative_count;
        }

        int index_count = 0;
        for (int r1 = r1_start; r1 < r1_finish; ++r1)
        {
            r1_supp = supports[r1];
            // TODO: Avoid unnecessary copying of rows.
            matrix = test_matrix;
            int r1_rows = upper_triangle(matrix, r1_supp, remaining_row);
            // Find the columns in the matrix which are zero.
            for (Index r2 = r2_start; r2 < r2_finish; ++r2)
            {
                ++num_adjacency_checks;
                IndexSet::set_difference(supports[r2], r1_supp, temp_diff);
                //if (temp_diff.count() <= codim-r1_rows+1)
                if (temp_diff.less_than_equal(codim-r1_rows+1))
                {
                    IndexSet::set_difference(supports[r2], r1_supp, temp_diff2);
                    if (temp_diff2.power_of_2())
                    {
                        create_new_vector(vs, supports, r1, r2, next_col,
                                            next_positive_count, next_negative_count,
                                            temp, temp_supp);
                        continue;
                    }
                    IndexSet::set_difference(r1_supp, supports[r2], temp_diff2);
                    if (temp_diff2.power_of_2())
                    {
                        create_new_vector(vs, supports, r1, r2, next_col,
                                            next_positive_count, next_negative_count,
                                            temp, temp_supp);
                        continue;
                    }
                    if (rank_check(matrix, temp_matrix, temp_diff, r1_rows))
                    {
                        create_new_vector(vs, supports, r1, r2, next_col,
                                        next_positive_count, next_negative_count,
                                        temp, temp_supp);
                        ++num_added;
                        continue;
                    }
                    ++num_dominated;
                }
                else
                {
                    ++num_high_dim;
                }
            }
            if (index_count % Globals::output_freq == 0)
            {
                *out << "\r" << buffer;
                *out << "  Size = " << std::setw(8) << vs.get_number()-next_negative_count << ", ";
                *out << "  Index = " << index_count << "/" << index_max << std::flush;
            }
            ++index_count;
        }
        // Update the support vectors for the next_col.
        for (int r1 = positive_start; r1 < negative_start; ++r1)
        {
            supports[r1].set(next_col);
        }

        // Delete all the vectors with a negative entry in the next_col.
        vs.remove(negative_start, original_size);
        supports.erase(supports.begin()+negative_start,
                        supports.begin()+original_size);

        //DEBUG_4ti2(*out << "Rank Checks " << num_rank_checks << "\n";)
        //*out << "Size = " << vs.get_number() << "\n" << std::endl;

        remaining.unset(next_col);
        --num_remaining;

        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << ", ";
        *out << "  Time: " << t << "                \n";
    }
    *out << "Added        " << num_added << "\n";
    *out << "Dominated    " << num_dominated << "\n";
    *out << "Num high dim " << num_high_dim << "\n";
    *out << "Total        " << num_adjacency_checks << "\n";
    return remaining;
}

// It is assumed that the cone is pointed.
// This is a straight-forward version of the ray computation algorithm.
template <class IndexSet>
IndexSet
RayMatrixAlgorithm<IndexSet>::compute0(
                VectorArray& matrix, 
                VectorArray& vs, 
                std::vector<IndexSet>& supports, 
                const IndexSet& rs)
{
    // Sanity Checks.
    assert(matrix.get_size() == vs.get_size());
    assert(rs.get_size() == vs.get_size());

    *out << "Ray Matrix Algorithm.\n";
    Timer t;

    int num_cols = vs.get_size();
    IndexSet urs(rs); // The variables that are unrestricted in sign.
    urs.set_complement();

    DEBUG_4ti2(*out << "RS\n" << rs << "\n";)
    DEBUG_4ti2(*out << "URS\n" << urs << "\n";)
    DEBUG_4ti2(*out << "The dimension is " << vs.get_number() << ".\n";)

    Index rows = diagonal(vs, rs); // Compute diagonal normal form.
    vs.remove(rows, vs.get_number());

    // We find the entries on the diagonal.
    supports.clear();
    Index col = 0;
    IndexSet diagonals(num_cols);
    for (Index r = 0; r < vs.get_number(); ++r)
    {
        while (vs[r][col] == 0 || urs[col]) { ++col; }
        diagonals.set(col);
        IndexSet support(num_cols, false);
        support.set(col);
        supports.push_back(support);
        ++col;
    }

    int codim = upper_triangle(matrix);
    matrix.remove(codim, matrix.get_number());
    DEBUG_4ti2(*out << "The codimension is " << codim << ".\n";)
    VectorArray orig_matrix(matrix);

    // The remaining columns to compute rays for.
    IndexSet remaining(rs);
    remaining.set_difference(diagonals);
    int num_remaining = remaining.count();

    // Temporary variables.
    IndexSet temp_supp(num_cols);
    IndexSet temp_diff(num_cols);
    IndexSet temp_diff2(num_cols);
    IndexSet temp_zero_cols(num_cols);
    IndexSet r1_supp(num_cols);
    VectorArray temp_matrix(matrix.get_number(), matrix.get_size(), 0);
    Vector temp(vs.get_size());

    unsigned long long int num_dominated = 0;
    unsigned long long int num_added = 0;
    unsigned long long int num_checks = 0;
    unsigned long long int num_one_diff_added = 0;
    while (vs.get_number() > 0 && num_remaining > 0)
    {
        // Find the next column.
        int next_positive_count, next_negative_count, next_zero_count;
        Index next_col = next_column(vs, remaining,
                                        next_positive_count,
                                        next_negative_count,
                                        next_zero_count);

        char buffer[256];
        sprintf(buffer, "  Left = %3d  Col = %3d", num_remaining, next_col);
        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << "  Time: " << t;
        DEBUG_4ti2(
            *out << "(+,0,-) = (" << next_positive_count << ",";
            *out << next_zero_count << "," << next_negative_count << ")\n";
        )

        // We sort the vectors into zeros, positives, then negatives.
        sort(vs, supports, next_col, next_zero_count, next_positive_count,
                        next_negative_count);

        matrix = orig_matrix;
        temp_supp = remaining;
        temp_supp.set_union(urs);
        int remaining_row = upper_triangle(matrix, temp_supp, 0);
        VectorArray test_matrix(matrix);

        int original_size = vs.get_number();
        int positive_start = next_zero_count;
        int negative_start = next_zero_count+next_positive_count;

        // We wish to reduce the number of matrices we triangularise so we
        // choose the smaller out of the positives and negatives.
        int r1_start;
        int r1_finish;
        int r2_start;
        int r2_finish;
        int index_max;
        if (next_positive_count <= next_negative_count)
        {
            //*out << "Using positive vectors.\n";
            r1_start = positive_start;
            r1_finish = negative_start;
            r2_start = negative_start;
            r2_finish = original_size;
            index_max = next_positive_count;
        }
        else
        {
            //*out << "Using negative vectors.\n";
            r1_start = negative_start;
            r1_finish = original_size;
            r2_start = positive_start;
            r2_finish = negative_start;
            index_max = next_negative_count;
        }

        int index_count = 0;
        for (int r1 = r1_start; r1 < r1_finish; ++r1)
        {
            r1_supp = supports[r1];
            if (r1_supp.count() == codim-num_remaining+1)
            {
                for (Index r2 = r2_start; r2 < r2_finish; ++r2)
                {
                    IndexSet::set_difference(supports[r2], r1_supp, temp_diff);
                    if (temp_diff.power_of_2())
                    {
                        create_new_vector(vs, supports, r1, r2, next_col,
                                            next_positive_count, next_negative_count,
                                            temp, temp_supp);
                        ++num_added;
                    }
                }
            }
            else
            {
                // TODO: Avoid unnecessary copying of rows.
                matrix = test_matrix;
                int r1_rows = upper_triangle(matrix, r1_supp, remaining_row);
                // Find the columns in the matrix which are zero.
                zero_cols(matrix, r1_supp, temp_zero_cols, r1_rows);
                for (Index r2 = r2_start; r2 < r2_finish; ++r2)
                {
                    if (IndexSet::set_disjoint(supports[r2], temp_zero_cols))
                    {
                        IndexSet::set_difference(supports[r2], r1_supp, temp_diff);
                        //if (temp_diff.count() <= codim-r1_rows+1)
                        if (temp_diff.less_than_equal(codim-r1_rows+1))
                        {
#if 1
                            IndexSet::set_difference(r1_supp, supports[r2], temp_diff2);
                            if (temp_diff2.power_of_2())
                            {
                                create_new_vector(vs, supports, r1, r2, next_col,
                                    next_positive_count, next_negative_count,
                                    temp, temp_supp);
                                ++num_one_diff_added;
                                continue;
                            }
#endif

                            if (rank_check(matrix, temp_matrix, temp_diff, r1_rows))
                            {
                                create_new_vector(vs, supports, r1, r2, next_col,
                                                next_positive_count, next_negative_count,
                                                temp, temp_supp);
                                ++num_added;
                            }
                            else
                            {
                                ++num_dominated;
                            }
                        }
                    }
                    else
                    {
                        IndexSet::set_difference(supports[r2], r1_supp, temp_diff);
                        if (temp_diff.power_of_2())
                        {
                            create_new_vector(vs, supports, r1, r2, next_col,
                                next_positive_count, next_negative_count,
                                temp, temp_supp);
                            ++num_added;
                        }
                        //else
                        //{
                        //    ++num_checks;
                        //}
                    }
                }
            }
            if (index_count % Globals::output_freq == 0)
            {
                *out << "\r" << buffer;
                *out << "  Size = " << std::setw(8) << vs.get_number()-next_negative_count;
                *out << ",   Index = " << index_count << "/" << index_max << std::flush;
            }
            ++index_count;
        }
        // Update the support vectors for the next_col.
        for (int r1 = positive_start; r1 < negative_start; ++r1)
        {
            supports[r1].set(next_col);
        }

        // Delete all the vectors with a negative entry in the next_col.
        vs.remove(negative_start, original_size);
        supports.erase(supports.begin()+negative_start,
                        supports.begin()+original_size);

        remaining.unset(next_col);
        --num_remaining;

        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << ", ";
        *out << "  Time: " << t << "                \n";

        *out << "Added      " << num_added << "\n";
        *out << "Dominated  " << num_dominated << "\n";
        *out << "Num Checks " << num_checks << "\n";
        *out << "One diff   " << num_one_diff_added << "\n";
    }
    return remaining;
}



// It is assumed that the cone is pointed.
// This version has adjacency propagation.
template <class IndexSet>
IndexSet
RayMatrixAlgorithm<IndexSet>::compute2(
                VectorArray& matrix, 
                VectorArray& vs, 
                std::vector<IndexSet>& supports, 
                const IndexSet& rs)
{
    // Sanity Checks.
    assert(matrix.get_size() == vs.get_size());
    assert(rs.get_size() == vs.get_size());

    *out << "Ray Matrix Algorithm.\n";
    Timer t;

    int num_cols = vs.get_size();
    IndexSet urs(rs); // The variables that are unrestricted in sign.
    urs.set_complement();

    DEBUG_4ti2(*out << "RS\n" << rs << "\n";)
    DEBUG_4ti2(*out << "URS\n" << urs << "\n";)
    DEBUG_4ti2(*out << "The dimension is " << vs.get_number() << ".\n";)

    Index rows = diagonal(vs, rs); // Compute diagonal normal form.
    vs.remove(rows, vs.get_number());

    // We find the entries on the diagonal.
    Fathers fathers;
    std::vector<IndexSet> zeros;
    supports.clear();
    Index col = 0;
    IndexSet diagonals(num_cols);
    for (Index r = 0; r < vs.get_number(); ++r)
    {
        while (vs[r][col] == 0 || urs[col]) { ++col; }
        diagonals.set(col);
        IndexSet support(num_cols, false);
        support.set(col);
        supports.push_back(support);
        fathers.push_back(-1);
        ++col;
    }
    for (Index r = 0; r < vs.get_number(); ++r)
    {
        IndexSet zero(diagonals);
        zero.set_difference(supports[r]);
        zeros.push_back(zero);
        ++col;
    }

    int codim = upper_triangle(matrix);
    matrix.remove(codim, matrix.get_number());
    *out << "The codimension is " << codim << ".\n";
    VectorArray orig_matrix(matrix);

    // The remaining columns to compute rays for.
    IndexSet remaining(rs);
    remaining.set_difference(diagonals);
    int num_remaining = remaining.count();

    // Temporary variables.
    IndexSet temp_supp(num_cols);
    IndexSet temp_diff(num_cols);
    IndexSet temp_diff2(num_cols);
    IndexSet temp_zero_cols(num_cols);
    IndexSet r1_supp(num_cols);
    IndexSet r1_zero(num_cols);
    VectorArray temp_matrix(matrix.get_number(), matrix.get_size(), 0);
    Vector temp(vs.get_size());

    while (vs.get_number() > 0 && num_remaining > 0)
    {
        // Statistics.
        long int num_dominated = 0;
        long int num_added = 0;
        long int num_one_diff_added = 0;
        //long int num_high_dim = 0;
        //long int num_zero_rows = 0;
        //long int num_rank_checks = 0;

        // Find the next column.
        int next_positive_count, next_negative_count, next_zero_count;
        Index next_col = next_column(vs, remaining,
                                        next_positive_count,
                                        next_negative_count,
                                        next_zero_count);

        char buffer[256];
        sprintf(buffer, "  Left = %3d  Col = %3d", num_remaining, next_col);
        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << "  Time: " << t;
        DEBUG_4ti2(
            *out << "(+,0,-) = (" << next_positive_count << ",";
            *out << next_zero_count << "," << next_negative_count << ")\n";
        )

        // We sort the vectors into zeros, positives, then negatives.
        sort(vs, supports, fathers, zeros, next_col, next_zero_count, next_positive_count,
                        next_negative_count);

        matrix = orig_matrix;
        temp_supp = remaining;
        temp_supp.set_union(urs);
        int remaining_row = upper_triangle(matrix, temp_supp, 0);
        VectorArray test_matrix(matrix);

        int original_size = vs.get_number();
        int positive_start = next_zero_count;
        int negative_start = next_zero_count+next_positive_count;

        // We wish to reduce the number of matrices we triangularise so we
        // choose the smaller out of the positives and negatives.
        int r1_start;
        int r1_finish;
        int r2_start;
        int r2_finish;
        int index_max;
        if (next_positive_count <= next_negative_count)
        {
            //*out << "Using positive vectors.\n";
            r1_start = positive_start;
            r1_finish = negative_start;
            r2_start = negative_start;
            r2_finish = original_size;
            index_max = next_positive_count;
        }
        else
        {
            //*out << "Using negative vectors.\n";
            r1_start = negative_start;
            r1_finish = original_size;
            r2_start = positive_start;
            r2_finish = negative_start;
            index_max = next_negative_count;
        }

        int index_count = 0;
        for (int r1 = r1_start; r1 < r1_finish; ++r1)
        {
            r1_supp = supports[r1];
            if (r1_supp.count() == codim-num_remaining+1)
            {
                for (Index r2 = r2_start; r2 < r2_finish; ++r2)
                {
                    IndexSet::set_difference(supports[r2], r1_supp, temp_diff);
                    if (temp_diff.power_of_2())
                    {
                        create_new_vector(vs, supports, fathers, r1, r2, next_col,
                                            next_positive_count, next_negative_count,
                                            temp, temp_supp);
                        ++num_added;
                    }
                }
            }
            else
            {
                // TODO: Avoid unnecessary copying of rows.
                r1_zero = zeros[r1];
                matrix = test_matrix;
                int r1_rows = upper_triangle(matrix, r1_supp, remaining_row);
                // Find the columns in the matrix which are zero.
                zero_cols(matrix, r1_supp, temp_zero_cols, r1_rows);
                for (Index r2 = r2_start; r2 < r2_finish; ++r2)
                {
                    if (IndexSet::set_disjoint(supports[r2], temp_zero_cols))
                    {
#if 1
                        if (IndexSet::set_disjoint(r1_zero, zeros[r2]))
                        {
                            if (fathers[r1] == r2 || fathers[r2] == r1)
                            {
                                create_new_vector(vs, supports, fathers, r1, r2, next_col,
                                    next_positive_count, next_negative_count,
                                    temp, temp_supp);
                                //++num_fathers;
                            }
                            continue;
                        }
#endif
                        IndexSet::set_difference(supports[r2], r1_supp, temp_diff);
                        //if (temp_diff.count() <= codim-r1_rows+1)
                        if (temp_diff.less_than_equal(codim-r1_rows+1))
                        {
#if 1
                            IndexSet::set_difference(r1_supp, supports[r2], temp_diff2);
                            if (temp_diff2.power_of_2())
                            {
                                create_new_vector(vs, supports, fathers, r1, r2, next_col,
                                    next_positive_count, next_negative_count,
                                    temp, temp_supp);
                                //++num_one_diff_added;
                                continue;
                            }
#endif
                            if (rank_check(matrix, temp_matrix, temp_diff, r1_rows))
                            {
                                create_new_vector(vs, supports, fathers, r1, r2, next_col,
                                                next_positive_count, next_negative_count,
                                                temp, temp_supp);
                                ++num_added;
                            }
                            else
                            {
                                ++num_dominated;
                            }
                        }
                    }
                    else
                    {
                        IndexSet::set_difference(supports[r2], r1_supp, temp_diff);
                        if (temp_diff.power_of_2())
                        {
                            create_new_vector(vs, supports, fathers, r1, r2, next_col,
                                next_positive_count, next_negative_count,
                                temp, temp_supp);
                            ++num_added;
                        }
                    }
                }
            }
            if (index_count % Globals::output_freq == 0)
            {
                *out << "\r" << buffer;
                *out << "  Size = " << std::setw(8) << vs.get_number()-next_negative_count << ", ";
                *out << "  Index = " << index_count << "/" << index_max << std::flush;
            }
            ++index_count;
        }
        // Update the support vectors for the next_col.
        for (int r1 = positive_start; r1 < negative_start; ++r1)
        {
            supports[r1].set(next_col);
        }
        for (int r1 = 0; r1 < positive_start; ++r1)
        {
            zeros[r1].set(next_col);
        }

        // Delete all the vectors with a negative entry in the next_col.
        int num_new_vectors = vs.get_number() - original_size;
        vs.remove(negative_start, original_size);
        supports.erase(supports.begin()+negative_start,
                        supports.begin()+original_size);
        fathers.erase(fathers.begin()+negative_start,
                        fathers.begin()+original_size);
        for (int i = 0; i < negative_start; ++i)
        {
            if (fathers[i] >= negative_start) { fathers[i] = -1; }
        }
        zeros.erase(zeros.begin()+negative_start,
                        zeros.begin()+original_size);
        IndexSet zero(num_cols);
        zero.set(next_col);
        zeros.insert(zeros.end(), num_new_vectors, zero);
        //for (int i = negative_start; i < fathers.size(); ++i)
        //{
        //    if (fathers[i] > negative_start || fathers[i] < positive_start)
        //    { std::cout << "Error\n"; }
        //}

        remaining.unset(next_col);
        --num_remaining;

        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << ", ";
        *out << "  Time: " << t << "                \n";

        //*out << "Added     " << num_added << "\n";
        //*out << "One diff  " << num_one_diff_added << "\n";
        //*out << "Dominated " << num_dominated << "\n";
        //DEBUG_4ti2(*out << "Rank Checks " << num_rank_checks << "\n";)
        //*out << "Num high dim " << num_high_dim << "\n";
        //*out << "Size = " << vs.get_number() << "\n" << std::endl;
    }
    return remaining;
}

// It is assumed that the cone is pointed.
template <class IndexSet>
IndexSet
RayMatrixAlgorithm<IndexSet>::compute3(
                VectorArray& matrix, 
                VectorArray& vs, 
                std::vector<IndexSet>& supports, 
                const IndexSet& rs)
{
    // Sanity Checks.
    assert(matrix.get_size() == vs.get_size());
    assert(rs.get_size() == vs.get_size());

    *out << "Ray Matrix Algorithm.\n";
    Timer t;

    int num_cols = vs.get_size();
    IndexSet urs(rs); // The variables that are unrestricted in sign.
    urs.set_complement();

    DEBUG_4ti2(*out << "RS\n" << rs << "\n";)
    DEBUG_4ti2(*out << "URS\n" << urs << "\n";)
    DEBUG_4ti2(*out << "The dimension is " << vs.get_number() << ".\n";)

    Index rows = diagonal(vs, rs); // Compute diagonal normal form.
    vs.remove(rows, vs.get_number());

    // We find the entries on the diagonal.
    supports.clear();
    Index col = 0;
    IndexSet diagonals(num_cols);
    for (Index r = 0; r < vs.get_number(); ++r)
    {
        while (vs[r][col] == 0 || urs[col]) { ++col; }
        diagonals.set(col);
        IndexSet support(num_cols, false);
        support.set(col);
        supports.push_back(support);
        ++col;
    }

    int codim = upper_triangle(matrix);
    matrix.remove(codim, matrix.get_number());
    DEBUG_4ti2(*out << "The codimension is " << codim << ".\n";)
    VectorArray orig_matrix(matrix);

    // The remaining columns to compute rays for.
    IndexSet remaining(rs);
    remaining.set_difference(diagonals);
    int num_remaining = remaining.count();

    // Temporary variables.
    IndexSet temp_supp(num_cols);
    IndexSet temp_diff(num_cols);
    IndexSet temp_diff2(num_cols);
    IndexSet temp_zero_cols(num_cols);
    IndexSet r1_supp(num_cols);
    VectorArray temp_matrix(matrix.get_number(), matrix.get_size(), 0);
    Vector temp(vs.get_size());

    DEBUG_4ti2(unsigned long long int num_dominated = 0;)
    DEBUG_4ti2(unsigned long long int num_added = 0;)
    DEBUG_4ti2(unsigned long long int num_one_diff_checks = 0;)
    DEBUG_4ti2(unsigned long long int num_one_diff_added = 0;)
    DEBUG_4ti2(unsigned long long int num_total = 0;)
    while (vs.get_number() > 0 && num_remaining > 0)
    {
        DEBUG_4ti2(unsigned long long int num_checks = 0;)
        // Find the next column.
        int next_positive_count, next_negative_count, next_zero_count;
        Index next_col = next_column(vs, remaining,
                                        next_positive_count,
                                        next_negative_count,
                                        next_zero_count);

        char buffer[256];
        sprintf(buffer, "  Left = %3d  Col = %3d", num_remaining, next_col);
        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << "  Time: " << t;
        DEBUG_4ti2(
            *out << "(+,0,-) = (" << next_positive_count << ",";
            *out << next_zero_count << "," << next_negative_count << ")\n";
        )

        // We sort the vectors into zeros, positives, then negatives.
        sort(vs, supports, next_col, next_zero_count, next_positive_count,
                        next_negative_count);

        DEBUG_4ti2(*out << "Rays:\n" << vs << "\n";)

        matrix = orig_matrix;
        temp_supp = remaining;
        temp_supp.set_union(urs);
        int remaining_row = upper_triangle(matrix, temp_supp, 0);
        VectorArray test_matrix(matrix);

        int original_size = vs.get_number();
        int positive_start = next_zero_count;
        int negative_start = next_zero_count+next_positive_count;

        // We wish to reduce the number of matrices we triangularise so we
        // choose the smaller out of the positives and negatives.
        int r1_start;
        int r1_finish;
        int r2_start;
        int r2_finish;
        int index_max;
        if (next_positive_count <= next_negative_count)
        {
            //*out << "Using positive vectors.\n";
            r1_start = positive_start;
            r1_finish = negative_start;
            r2_start = negative_start;
            r2_finish = original_size;
            index_max = next_positive_count;
        }
        else
        {
            //*out << "Using negative vectors.\n";
            r1_start = negative_start;
            r1_finish = original_size;
            r2_start = positive_start;
            r2_finish = negative_start;
            index_max = next_negative_count;
        }
        // We sort the r2's into vectors where r2_supp.count()==codim-num_remaining+1.
        int r2_index = r2_start;
        for (int r2 = r2_start; r2 < r2_finish; ++r2)
        {
            if (supports[r2].count() == codim-num_remaining+1)
            {
                vs.swap_vectors(r2, r2_index);
                IndexSet::swap(supports[r2], supports[r2_index]);
                ++r2_index;
            }
        }
        //std::cout << "\nR2: " << r2_start << " " << r2_index << " " << r2_finish << "\n";

        int index_count = 0;
        for (int r1 = r1_start; r1 < r1_finish; ++r1)
        {
            r1_supp = supports[r1];
            if (r1_supp.count() == codim-num_remaining+1)
            {
                for (Index r2 = r2_start; r2 < r2_finish; ++r2)
                {
                    IndexSet::set_difference(supports[r2], r1_supp, temp_diff);
                    if (temp_diff.power_of_2())
                    {
                        create_new_vector(vs, supports, r1, r2, next_col,
                                            next_positive_count, next_negative_count,
                                            temp, temp_supp);
                        DEBUG_4ti2(++num_one_diff_added;)
                        DEBUG_4ti2(++num_added;)
                    }
                }
            }
            else
            {
                for (Index r2 = r2_start; r2 < r2_index; ++r2)
                {
                    IndexSet::set_difference(r1_supp, supports[r2], temp_diff);
                    if (temp_diff.power_of_2())
                    {
                        create_new_vector(vs, supports, r1, r2, next_col,
                                            next_positive_count, next_negative_count,
                                            temp, temp_supp);
                        DEBUG_4ti2(++num_one_diff_added;)
                        DEBUG_4ti2(++num_added;)
                    }
                }

                // TODO: Avoid unnecessary copying of rows.
                matrix = test_matrix;
                int r1_rows = upper_triangle(matrix, r1_supp, remaining_row);
                // Find the columns in the matrix which are zero.
                zero_cols(matrix, r1_supp, temp_zero_cols, r1_rows);
                for (Index r2 = r2_index; r2 < r2_finish; ++r2)
                {
                    if (IndexSet::set_disjoint(supports[r2], temp_zero_cols))
                    {
                        IndexSet::set_difference(supports[r2], r1_supp, temp_diff);
                        //if (temp_diff.count() <= codim-r1_rows+1)
                        if (temp_diff.less_than_equal(codim-r1_rows+1))
                        {
#if 1
                            IndexSet::set_difference(r1_supp, supports[r2], temp_diff2);
                            if (temp_diff2.power_of_2())
                            {
                                create_new_vector(vs, supports, r1, r2, next_col,
                                    next_positive_count, next_negative_count,
                                    temp, temp_supp);
                                DEBUG_4ti2(++num_one_diff_added;)
                                DEBUG_4ti2(++num_added;)
                                continue;
                            }
#endif

                            DEBUG_4ti2(++num_checks;)
                            if (rank_check(matrix, temp_matrix, temp_diff, r1_rows))
                            {
                                create_new_vector(vs, supports, r1, r2, next_col,
                                                next_positive_count, next_negative_count,
                                                temp, temp_supp);
                                DEBUG_4ti2(++num_added;)
                            }
                            else
                            {
                                DEBUG_4ti2(++num_dominated;)
                            }
                        }
                    }
                    else
                    {
                        IndexSet::set_difference(supports[r2], r1_supp, temp_diff);
                        if (temp_diff.power_of_2())
                        {
                            create_new_vector(vs, supports, r1, r2, next_col,
                                next_positive_count, next_negative_count,
                                temp, temp_supp);
                            DEBUG_4ti2(++num_one_diff_added;)
                            DEBUG_4ti2(++num_added;)
                        }
                    }
                }
            }
            if (index_count % Globals::output_freq == 0)
            {
                *out << "\r" << buffer;
                *out << "  Size = " << std::setw(8) << vs.get_number()-next_negative_count;
                *out << ",   Index = " << index_count << "/" << index_max << std::flush;
            }
            ++index_count;
        }
        // Update the support vectors for the next_col.
        for (int r1 = positive_start; r1 < negative_start; ++r1)
        {
            supports[r1].set(next_col);
        }

        // Delete all the vectors with a negative entry in the next_col.
        vs.remove(negative_start, original_size);
        supports.erase(supports.begin()+negative_start,
                        supports.begin()+original_size);

        remaining.unset(next_col);
        --num_remaining;

        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << ", ";
        *out << "  Time: " << t << "                \n";

        DEBUG_4ti2(*out << "Added      " << num_added << "\n";)
        DEBUG_4ti2(*out << "Dominated  " << num_dominated << "\n";)
        DEBUG_4ti2(*out << "Num Checks " << num_checks << "\n";)
        DEBUG_4ti2(*out << "One diff   " << num_one_diff_added << "\n";)
        DEBUG_4ti2(num_one_diff_checks += next_positive_count + next_negative_count;)
       DEBUG_4ti2( *out << "One Checks " << num_one_diff_checks << "\n";)
    }
    return remaining;
}

template <class IndexSet>
bool
RayMatrixAlgorithm<IndexSet>::rank_check(
                VectorArray& matrix,
                VectorArray& _temp_matrix,
                IndexSet& temp_diff,
                int r1_rows)
{
    int m = matrix.get_number()-r1_rows;
    int n = temp_diff.count();
    VectorArray temp_matrix(m,n); // TODO: Handle memory management better.
    int col_index = 0;
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
    DEBUG_4ti2(*out << "\nTemp Matrix:\n" << temp_matrix << "\n";)
    int rank = upper_triangle(temp_matrix);
    return (rank == n-1);
}

template <class IndexSet>
void
RayMatrixAlgorithm<IndexSet>::zero_cols(
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
