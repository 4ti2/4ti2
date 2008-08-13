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

#include "groebner/RaySupportAlgorithm.h"
#include "groebner/DiagonalAlgorithm.h"
#include "groebner/HermiteAlgorithm.h"
#include "groebner/Euclidean.h"
#include "groebner/SupportTree.h"
#include "groebner/OnesTree.h"
#include "groebner/Globals.h"
#include "groebner/Timer.h"
#include "groebner/Debug.h"

#include "groebner/VectorArrayStream.h"
#include "groebner/VectorStream.h"
#include "groebner/VectorArrayStream.h"
#include <iostream>
#include <iomanip>

#define TREE SupportTree

using namespace _4ti2_;

template <class IndexSet>
RaySupportAlgorithm<IndexSet>::RaySupportAlgorithm()
{
}

template <class IndexSet>
RaySupportAlgorithm<IndexSet>::~RaySupportAlgorithm()
{
}

template <class IndexSet>
IndexSet
RaySupportAlgorithm<IndexSet>::compute(
                VectorArray& matrix,
                VectorArray& vs,
                std::vector<IndexSet>& supps,
                const IndexSet& rs)
{
    return compute3(matrix, vs, supps, rs);
}

template <class IndexSet>
IndexSet
RaySupportAlgorithm<IndexSet>::compute(
                VectorArray& matrix,
                VectorArray& vs,
                const IndexSet& rs)
{
    std::vector<IndexSet> supports;
    return compute(matrix, vs, supports, rs);
}

// It is assumed that the cone is pointed.
template <class IndexSet>
IndexSet
RaySupportAlgorithm<IndexSet>::compute0(
                VectorArray& matrix,
                VectorArray& vs,
                std::vector<IndexSet>& supports,
                const IndexSet& rs)
{
    // Sanity Checks.
    assert(matrix.get_size() == vs.get_size());
    assert(rs.get_size() == vs.get_size());

    *out << "Ray Support Algorithm.\n";
    Timer t;

    int num_cols = vs.get_size();
    IndexSet urs(rs); // The variables that are unrestricted in sign.
    urs.set_complement();

    DEBUG_4ti2(*out << "RS\n" << rs << "\n";)
    DEBUG_4ti2(*out << "URS\n" << urs << "\n";)
    DEBUG_4ti2(*out << "The dimension is " << vs.get_number() << "\n";)

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
    DEBUG_4ti2(*out << "The codimension is " << codim << "\n";)
    VectorArray orig_matrix(matrix);

    // The remaining columns to compute rays for.
    IndexSet remaining(rs);
    remaining.set_difference(diagonals);
    int num_remaining = remaining.count();

    // The columns with relaxed non-negativity constraints.
    IndexSet relaxed(remaining);
    relaxed.set_union(urs);
    int num_relaxed = relaxed.count();

    // Temporary variables.
    IndexSet temp_supp(num_cols);
    IndexSet temp_diff(num_cols);
    IndexSet r1_supp(num_cols);
    Vector temp(vs.get_size());

    while (vs.get_number() > 0 && num_remaining > 0)
    {
        // Statistics.
        long int num_added = 0;
        long int num_support_checks = 0;

        // Find the next column.
        int next_positive_count, next_negative_count, next_zero_count;
        Index next_col = next_column(vs, remaining,
                                        next_positive_count,
                                        next_negative_count,
                                        next_zero_count);

        char buffer[256];
        sprintf(buffer, "  Left = %3d,  Col = %3d,", num_remaining, next_col);
        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << "  Time: " << t;
        DEBUG_4ti2(
            *out << "(+,0,-) = (" << next_positive_count << ",";
            *out << next_zero_count << "," << next_negative_count << ")\n";
        )

        // We sort the vectors into zeros, positives, then negatives.
        sort(vs, supports, next_col, next_zero_count, next_positive_count,
                        next_negative_count);

        // Note that the tree needs the ordering of the current vectors to be
        // constant.
        DEBUG_4ti2(*out << "Building Tree ... " << std::flush;)
        TREE<IndexSet> tree(supports, supports.size());
        DEBUG_4ti2(*out << "done." << std::endl;)
        DEBUG_4ti2(*out << "VS:\n" << vs << "\n";)
        DEBUG_4ti2(tree.dump();)

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
            int r1_count = r1_supp.count();
            if (r1_count == codim-num_relaxed+1)
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
                for (Index r2 = r2_start; r2 < r2_finish; ++r2)
                {
                    IndexSet::set_difference(supports[r2],supports[r1],temp_supp);
                    //if (temp_supp.count() <= codim-num_relaxed-r1_count+2) 
                    if (temp_supp.less_than_equal(codim-num_relaxed-r1_count+2))
                    {
                        IndexSet::set_union(r1_supp,supports[r2],temp_supp);
                        ++num_support_checks;
                        if (!tree.dominated(temp_supp, r1, r2))
                        {
                            create_new_vector(vs, supports, r1, r2, next_col,
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

        // Delete all the vectors with a negative entry in the next_col.
        vs.remove(negative_start, original_size);
        supports.erase(supports.begin()+negative_start,
                        supports.begin()+original_size);

        DEBUG_4ti2(*out << "Added  " << num_added << "\n";)
        DEBUG_4ti2(*out << "Support Checks " << num_support_checks << "\n";)

        remaining.unset(next_col);
        --num_remaining;
        relaxed.unset(next_col);
        --num_relaxed;

        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << ",";
        *out << "  Time: " << t << "                \n";
    }
    return remaining;
}

// It is assumed that the cone is pointed.
// Finds the rays whose support differs by one from the current ray's support.
template <class IndexSet>
IndexSet
RaySupportAlgorithm<IndexSet>::compute1(
                VectorArray& matrix,
                VectorArray& vs,
                std::vector<IndexSet>& supports,
                const IndexSet& rs)
{
    // Sanity Checks.
    assert(matrix.get_size() == vs.get_size());
    assert(rs.get_size() == vs.get_size());

    *out << "Ray Support Algorithm.\n";
    Timer t;

    int num_cols = vs.get_size();
    IndexSet urs(rs); // The variables that are unrestricted in sign.
    urs.set_complement();

    DEBUG_4ti2(*out << "The dimension is " << vs.get_number() << "\n";)

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
    DEBUG_4ti2(*out << "The codimension is " << codim << "\n";)
    VectorArray orig_matrix(matrix);

    // The remaining columns to compute rays for.
    IndexSet remaining(rs);
    remaining.set_difference(diagonals);
    int num_remaining = remaining.count();

    // The columns with relaxed non-negativity constraints.
    IndexSet relaxed(remaining);
    relaxed.set_union(urs);
    int num_relaxed = relaxed.count();

    // Temporary variables.
    IndexSet temp_supp(num_cols);
    IndexSet temp_diff(num_cols);
    IndexSet zero_supp(num_cols);
    IndexSet r1_supp(num_cols);
    Vector temp(vs.get_size());
    std::vector<int> indices;

    //*out << "Matrix:\n" << matrix << "\n";
    //*out << "VS:\n" << vs << "\n";

    while (vs.get_number() > 0 && num_remaining > 0)
    {
        // Statistics.
        long int num_added = 0;
        long int num_support_checks = 0;

        // Find the next column.
        int next_positive_count, next_negative_count, next_zero_count;
        Index next_col = next_column(vs, remaining,
                                        next_positive_count,
                                        next_negative_count,
                                        next_zero_count);

        char buffer[256];
        sprintf(buffer, "  Left = %3d,  Col = %3d,", num_remaining, next_col);
        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << "  Time: " << t;
        DEBUG_4ti2(
            *out << "(+,0,-) = (" << next_positive_count << ",";
            *out << next_zero_count << "," << next_negative_count << ")\n";
        )

        // We sort the vectors into zeros, positives, then negatives.
        sort(vs, supports, next_col, next_zero_count, next_positive_count,
                        next_negative_count);

        // Note that the tree needs the ordering of the current vectors to be
        // constant.
        DEBUG_4ti2(*out << "Building Tree ... " << std::flush;)
        TREE<IndexSet> tree(supports, supports.size());
        DEBUG_4ti2(*out << "done." << std::endl;)
        DEBUG_4ti2(tree.dump();)

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
            int r1_count = r1_supp.count();
            if (r1_count == codim-num_relaxed+1)
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
                zero_supp.zero();
                indices.clear();
                IndexSet::set_complement(r1_supp, temp_supp);
                // Find the rays whose support differs by one from the current ray's support.
                tree.find_diff(indices, temp_supp, 1);
                for (unsigned int i = 0; i < indices.size(); ++i)
                {
                    int index = indices[i];
                    zero_supp.set_union(supports[index]);
                    if (index >= r2_start && index < r2_finish)
                    {
                        create_new_vector(vs, supports, r1, index, next_col,
                            next_positive_count, next_negative_count,
                            temp, temp_supp);
                        ++num_added;
                    }
                }
                zero_supp.set_difference(r1_supp);
                for (Index r2 = r2_start; r2 < r2_finish; ++r2)
                {
                    if (!IndexSet::set_disjoint(zero_supp, supports[r2])) { continue; }
                    IndexSet::set_difference(supports[r2],supports[r1],temp_supp);
                    //if (temp_supp.count() <= codim-num_relaxed-r1_count+2) 
                    if (temp_supp.less_than_equal(codim-num_relaxed-r1_count+2))
                    {
                        IndexSet::set_union(r1_supp,supports[r2],temp_supp);
                        ++num_support_checks;
                        if (!tree.dominated(temp_supp, r1, r2))
                        {
                            create_new_vector(vs, supports, r1, r2, next_col,
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

        // Delete all the vectors with a negative entry in the next_col.
        vs.remove(negative_start, original_size);
        supports.erase(supports.begin()+negative_start,
                        supports.begin()+original_size);

        DEBUG_4ti2(*out << "Added  " << num_added << "\n";)
        DEBUG_4ti2(*out << "Support Checks " << num_support_checks << "\n";)

        remaining.unset(next_col);
        --num_remaining;
        relaxed.unset(next_col);
        --num_relaxed;

        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << ",";
        *out << "  Time: " << t << "                \n";
    }
    return remaining;
}

// It is assumed that the cone is pointed.
// Uses look ahead.
template <class IndexSet>
IndexSet
RaySupportAlgorithm<IndexSet>::compute2(
                VectorArray& matrix,
                VectorArray& vs,
                std::vector<IndexSet>& supports,
                const IndexSet& rs)
{
    // Sanity Checks.
    assert(matrix.get_size() == vs.get_size());
    assert(rs.get_size() == vs.get_size());

    *out << "Ray Support Algorithm.\n";
    Timer t;

    int num_cols = vs.get_size();
    IndexSet urs(rs); // The variables that are unrestricted in sign.
    urs.set_complement();

    DEBUG_4ti2(*out << "The dimension is " << vs.get_number() << "\n";)

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
    DEBUG_4ti2(*out << "The codimension is " << codim << "\n";)
    VectorArray orig_matrix(matrix);

    // The remaining columns to compute rays for.
    IndexSet remaining(rs);
    remaining.set_difference(diagonals);
    int num_remaining = remaining.count();

    // The columns with relaxed non-negativity constraints.
    IndexSet relaxed(remaining);
    relaxed.set_union(urs);
    int num_relaxed = relaxed.count();

    // Temporary variables.
    IndexSet temp_supp(num_cols);
    IndexSet temp_diff(num_cols);
    IndexSet zero_supp(num_cols);
    IndexSet r1_supp(num_cols);
    Vector temp(vs.get_size());
    std::vector<int> indices;

    while (vs.get_number() > 0 && num_remaining > 0)
    {
        // Statistics.
        long int num_added = 0;
        long int num_support_checks = 0;

        // Find the next column.
        int next_positive_count, next_negative_count, next_zero_count;
        Index next_col = next_column(vs, remaining,
                                        next_positive_count,
                                        next_negative_count,
                                        next_zero_count);

        char buffer[256];
        sprintf(buffer, "  Left = %3d,  Col = %3d,", num_remaining, next_col);
        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << "  Time: " << t;
        DEBUG_4ti2(
            *out << "(+,0,-) = (" << next_positive_count << ",";
            *out << next_zero_count << "," << next_negative_count << ")\n";
        )

        // We sort the vectors into zeros, positives, then negatives.
        sort(vs, supports, next_col, next_zero_count, next_positive_count,
                        next_negative_count);

        // Note that the tree needs the ordering of the current vectors to be
        // constant.
        DEBUG_4ti2(*out << "Building Tree ... " << std::flush;)
        TREE<IndexSet> tree(supports, supports.size());
        DEBUG_4ti2(*out << "done." << std::endl;)
        //tree.dump();

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
        TREE<IndexSet> small_tree;
        for (int i = r2_start; i < r2_finish; ++i) { small_tree.insert(supports[i], i); }

        int index_count = 0;
        for (int r1 = r1_start; r1 < r1_finish; ++r1)
        {
            r1_supp = supports[r1];
            int r1_count = r1_supp.count();
            if (r1_count == codim-num_relaxed+1)
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
                zero_supp.zero();
                indices.clear();
                IndexSet::set_complement(r1_supp, temp_supp);
                tree.find_diff(indices, temp_supp, 1);
                for (unsigned int i = 0; i < indices.size(); ++i)
                {
                    int index = indices[i];
                    zero_supp.set_union(supports[index]);
                    if (index >= r2_start && index < r2_finish)
                    {
                        create_new_vector(vs, supports, r1, index, next_col,
                            next_positive_count, next_negative_count,
                            temp, temp_supp);
                        ++num_added;
                    }
                }
                zero_supp.set_difference(r1_supp);
                indices.clear();
                IndexSet::set_complement(r1_supp, temp_supp);
                small_tree.find_diff(indices, zero_supp, 0, temp_supp, codim-num_relaxed-r1_count+2);
                //*out << " " << indices.size();
                for (unsigned int i = 0; i < indices.size(); ++i)
                {
                    int r2 = indices[i];
                    //IndexSet::set_difference(supports[r2],supports[r1],temp_supp);
                    //if (temp_supp.count() <= codim-num_relaxed-r1_count+2) 
                    if (temp_supp.less_than_equal(codim-num_relaxed-r1_count+2))
                    {
                        IndexSet::set_union(r1_supp,supports[r2],temp_supp);
                        ++num_support_checks;
                        if (!tree.dominated(temp_supp, r1, r2))
                        {
                            create_new_vector(vs, supports, r1, r2, next_col,
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

        // Delete all the vectors with a negative entry in the next_col.
        vs.remove(negative_start, original_size);
        supports.erase(supports.begin()+negative_start,
                        supports.begin()+original_size);

        DEBUG_4ti2(*out << "Added  " << num_added << "\n";)
        DEBUG_4ti2(*out << "Support Checks " << num_support_checks << "\n";)

        remaining.unset(next_col);
        --num_remaining;
        relaxed.unset(next_col);
        --num_relaxed;

        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << ",";
        *out << "  Time: " << t << "                \n";
    }
    return remaining;
}

// It is assumed that the cone is pointed.
// Finds the rays whose support differs by one from the current ray's support.
template <class IndexSet>
IndexSet
RaySupportAlgorithm<IndexSet>::compute3(
                VectorArray& matrix,
                VectorArray& vs,
                std::vector<IndexSet>& supports,
                const IndexSet& rs)
{
    // Sanity Checks.
    assert(matrix.get_size() == vs.get_size());
    assert(rs.get_size() == vs.get_size());

    *out << "Ray Support Algorithm.\n";
    Timer t;

    int num_cols = vs.get_size();
    IndexSet urs(rs); // The variables that are unrestricted in sign.
    urs.set_complement();

    DEBUG_4ti2(*out << "The dimension is " << vs.get_number() << "\n";)

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
    DEBUG_4ti2(*out << "The codimension is " << codim << "\n";)
    VectorArray orig_matrix(matrix);

    // The remaining columns to compute rays for.
    IndexSet remaining(rs);
    remaining.set_difference(diagonals);
    int num_remaining = remaining.count();

    // The columns with relaxed non-negativity constraints.
    IndexSet relaxed(remaining);
    relaxed.set_union(urs);
    int num_relaxed = relaxed.count();

    // Temporary variables.
    IndexSet temp_supp(num_cols);
    IndexSet temp_diff(num_cols);
    IndexSet temp_diff2(num_cols);
    IndexSet zero_supp(num_cols);
    IndexSet r1_supp(num_cols);
    Vector temp(vs.get_size());
    std::vector<int> indices;

    //*out << "Matrix:\n" << matrix << "\n";
    //*out << "VS:\n" << vs << "\n";

    while (vs.get_number() > 0 && num_remaining > 0)
    {
        // Statistics.
        DEBUG_4ti2(unsigned long long int num_added = 0;)
        DEBUG_4ti2(unsigned long long int num_checks = 0;)

        // Find the next column.
        int next_positive_count, next_negative_count, next_zero_count;
        Index next_col = next_column(vs, remaining,
                                        next_positive_count,
                                        next_negative_count,
                                        next_zero_count);

        char buffer[256];
        sprintf(buffer, "  Left = %3d,  Col = %3d,", num_remaining, next_col);
        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << "  Time: " << t;
        DEBUG_4ti2(
            *out << "(+,0,-) = (" << next_positive_count << ",";
            *out << next_zero_count << "," << next_negative_count << ")\n";
        )

        // We sort the vectors into zeros, positives, then negatives.
        sort(vs, supports, next_col, next_zero_count, next_positive_count,
                        next_negative_count);

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
        // We sort the r2's into vectors where r2_supp.count()==codim-num_relaxed+1.
        int r2_index = r2_start;
        for (int r2 = r2_start; r2 < r2_finish; ++r2)
        {
            if (supports[r2].count() == codim-num_relaxed+1)
            {
                vs.swap_vectors(r2, r2_index);
                IndexSet::swap(supports[r2], supports[r2_index]);
                ++r2_index;
            }
        }
        //std::cout << "\nR2: " << r2_start << " " << r2_index << " " << r2_finish << "\n";

        // Note that the tree needs the ordering of the current vectors to be
        // constant.
        DEBUG_4ti2(*out << "Building Tree ... " << std::flush;)
        TREE<IndexSet> tree(supports, supports.size());
        DEBUG_4ti2(*out << "done." << std::endl;)
        DEBUG_4ti2(tree.dump();)

        int index_count = 0;
        for (int r1 = r1_start; r1 < r1_finish; ++r1)
        {
            r1_supp = supports[r1];
            int r1_count = r1_supp.count();
            if (r1_count == codim-num_relaxed+1)
            {
                for (Index r2 = r2_start; r2 < r2_finish; ++r2)
                {
                    IndexSet::set_difference(supports[r2], r1_supp, temp_diff);
                    if (temp_diff.power_of_2())
                    {
                        create_new_vector(vs, supports, r1, r2, next_col,
                                            next_positive_count, next_negative_count,
                                            temp, temp_supp);
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
                        DEBUG_4ti2(++num_added;)
                    }
                }

                zero_supp.zero();
                indices.clear();
                IndexSet::set_complement(r1_supp, temp_supp);
                // Find the rays whose support differs by one from the current ray's support.
                tree.find_diff(indices, temp_supp, 1);
                for (unsigned int i = 0; i < indices.size(); ++i)
                {
                    int index = indices[i];
                    zero_supp.set_union(supports[index]);
                    if (index >= r2_index && index < r2_finish)
                    {
                        create_new_vector(vs, supports, r1, index, next_col,
                            next_positive_count, next_negative_count,
                            temp, temp_supp);
                        DEBUG_4ti2(++num_added;)
                    }
                }
                zero_supp.set_difference(r1_supp);

                for (Index r2 = r2_index; r2 < r2_finish; ++r2)
                {
                    if (!IndexSet::set_disjoint(zero_supp, supports[r2])) { continue; }
                    IndexSet::set_difference(supports[r2],r1_supp,temp_supp);
                    //if (temp_supp.count() <= codim-num_relaxed-r1_count+2) 
                    if (temp_supp.less_than_equal(codim-num_relaxed-r1_count+2))
                    {
#if 1
                        IndexSet::set_difference(r1_supp, supports[r2], temp_diff2);
                        if (temp_diff2.power_of_2())
                        {
                            create_new_vector(vs, supports, r1, r2, next_col,
                                    next_positive_count, next_negative_count,
                                    temp, temp_supp);
                            DEBUG_4ti2(++num_added;)
                            continue;
                        }
#endif

                        IndexSet::set_union(r1_supp,supports[r2],temp_supp);
                        DEBUG_4ti2(++num_checks;)
                        if (!tree.dominated(temp_supp, r1, r2))
                        {
                            create_new_vector(vs, supports, r1, r2, next_col,
                                            next_positive_count, next_negative_count,
                                            temp, temp_supp);
                            DEBUG_4ti2(++num_added;)
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

        // Delete all the vectors with a negative entry in the next_col.
        vs.remove(negative_start, original_size);
        supports.erase(supports.begin()+negative_start,
                        supports.begin()+original_size);

        remaining.unset(next_col);
        --num_remaining;
        relaxed.unset(next_col);
        --num_relaxed;

        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << ",";
        *out << "  Time: " << t << "                \n";

        DEBUG_4ti2(*out << "Num Checks " << num_checks << "\n";)
        DEBUG_4ti2(*out << "Added  " << num_added << "\n";)
    }
    return remaining;
}


// It is assumed that the cone is pointed.
// Perform extreme ray support check as opposed to adjacency support check.
template <class IndexSet>
IndexSet
RaySupportAlgorithm<IndexSet>::compute4(
                VectorArray& matrix,
                VectorArray& vs,
                std::vector<IndexSet>& supports,
                const IndexSet& rs)
{
    // Sanity Checks.
    assert(matrix.get_size() == vs.get_size());
    assert(rs.get_size() == vs.get_size());

    *out << "Ray Support Algorithm.\n";
    Timer t;

    int num_cols = vs.get_size();
    IndexSet urs(rs); // The variables that are unrestricted in sign.
    urs.set_complement();

    DEBUG_4ti2(*out << "The dimension is " << vs.get_number() << "\n";)

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
    DEBUG_4ti2(*out << "The codimension is " << codim << "\n";)
    VectorArray orig_matrix(matrix);

    // The remaining columns to compute rays for.
    IndexSet remaining(rs);
    remaining.set_difference(diagonals);
    int num_remaining = remaining.count();

    // The columns with relaxed non-negativity constraints.
    IndexSet relaxed(remaining);
    relaxed.set_union(urs);
    int num_relaxed = relaxed.count();

    // Temporary variables.
    IndexSet temp_supp(num_cols);
    IndexSet temp_diff(num_cols);
    IndexSet temp_diff2(num_cols);
    IndexSet zero_supp(num_cols);
    IndexSet r1_supp(num_cols);
    Vector temp(vs.get_size());
    std::vector<int> indices;

    //*out << "Matrix:\n" << matrix << "\n";
    //*out << "VS:\n" << vs << "\n";

    while (vs.get_number() > 0 && num_remaining > 0)
    {
        // Statistics.
        long int num_added = 0;
        long int num_support_checks = 0;

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
#if 1
        // We sort the r2's into vectors where r2_supp.count()==codim-num_relaxed+1.
        int r2_index = r2_start;
        for (int r2 = r2_start; r2 < r2_finish; ++r2)
        {
            if (supports[r2].count() == codim-num_relaxed+1)
            {
                vs.swap_vectors(r2, r2_index);
                IndexSet::swap(supports[r2], supports[r2_index]);
                ++r2_index;
            }
        }
        //std::cout << "\nR2: " << r2_start << " " << r2_index << " " << r2_finish << "\n";
#endif

        // Note that the tree needs the ordering of the current vectors to be
        // constant.
        DEBUG_4ti2(*out << "Building Tree ... " << std::flush;)
        TREE<IndexSet> tree(supports, supports.size());
        DEBUG_4ti2(*out << "done." << std::endl;)
        DEBUG_4ti2(tree.dump();)

        std::vector<std::vector<IndexSet> > next_supports(num_cols);
        std::vector<std::vector<std::pair<int,int> > >next_indices(num_cols);

        int index_count = 0;
        for (int r1 = r1_start; r1 < r1_finish; ++r1)
        {
            r1_supp = supports[r1];
            int r1_count = r1_supp.count();
            if (r1_count == codim-num_relaxed+1)
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
                zero_supp.zero();
                indices.clear();
                IndexSet::set_complement(r1_supp, temp_supp);
                // Find the rays whose support differs by one from the current ray's support.
                tree.find_diff(indices, temp_supp, 1);
                for (unsigned int i = 0; i < indices.size(); ++i)
                {
                    int index = indices[i];
                    zero_supp.set_union(supports[index]);
                    if (index >= r2_start && index < r2_finish)
                    {
                        create_new_vector(vs, supports, r1, index, next_col,
                            next_positive_count, next_negative_count,
                            temp, temp_supp);
                        ++num_added;
                    }
                }
                zero_supp.set_difference(r1_supp);
#if 1
                for (Index r2 = r2_start; r2 < r2_index; ++r2)
                {
                    if (!IndexSet::set_disjoint(zero_supp, supports[r2])) { continue; }
                    IndexSet::set_difference(r1_supp, supports[r2], temp_diff);
                    if (temp_diff.power_of_2())
                    {
                        create_new_vector(vs, supports, r1, r2, next_col,
                                            next_positive_count, next_negative_count,
                                            temp, temp_supp);
                        DEBUG_4ti2(++num_added;)
                    }
                }
#endif
                for (Index r2 = r2_index; r2 < r2_finish; ++r2)
                {
                    if (!IndexSet::set_disjoint(zero_supp, supports[r2])) { continue; }
                    IndexSet::set_difference(supports[r2],supports[r1],temp_supp);
                    //if (temp_supp.count() <= codim-num_relaxed-r1_count+2) 
                    if (temp_supp.less_than_equal(codim-num_relaxed-r1_count+2))
                    {
#if 1
                        IndexSet::set_difference(r1_supp, supports[r2], temp_diff2);
                        if (temp_diff2.power_of_2())
                        {
                            create_new_vector(vs, supports, r1, r2, next_col,
                                    next_positive_count, next_negative_count,
                                    temp, temp_supp);
                            DEBUG_4ti2(++num_added;)
                            continue;
                        }
#endif
                        IndexSet::set_union(r1_supp,supports[r2],temp_supp);
                        ++num_support_checks;
                        //if (!tree.dominated(temp_supp, r1, r2))
                        {
                            int count = temp_supp.count();
                            next_supports[count].push_back(temp_supp);
                            next_indices[count].push_back(std::pair<int,int>(r1,r2));
                            //create_new_vector(vs, supports, r1, r2, next_col,
                            //                next_positive_count, next_negative_count,
                            //                temp, temp_supp);
                            //++num_added;
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

        std::cout << "\nNum checks: " << num_support_checks << "\n";
        TREE<IndexSet> next_tree;
        // Inserting the rays with zeros in the next_col component.
        for (int i = 0; i < next_zero_count; ++i)
        { next_tree.insert(supports[i], i); }
        // Inserting the new rays.
        for (int i = original_size; i < supports.size(); ++i)
        { next_tree.insert(supports[i], i); }
        for (int i = 0; i < num_cols; ++i)
        {
            for (int j = 0; j < next_supports[i].size(); ++j)
            {
                if (!next_tree.dominated(next_supports[i][j]))
                {
                    int r1 = next_indices[i][j].first;
                    int r2 = next_indices[i][j].second;
                    create_new_vector(vs, supports, r1, r2, next_col,
                                    next_positive_count, next_negative_count,
                                    temp, temp_supp);
                    next_tree.insert(next_supports[i][j], vs.get_number()-1);
                    //++num_added;
                }
            }
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

        DEBUG_4ti2(*out << "Added  " << num_added << "\n";)
        DEBUG_4ti2(*out << "Support Checks " << num_support_checks << "\n";)

        remaining.unset(next_col);
        --num_remaining;
        relaxed.unset(next_col);
        --num_relaxed;

        *out << "\r" << buffer;
        *out << "  Size = " << std::setw(8) << vs.get_number() << ",";
        *out << "  Time: " << t << "                \n";
    }
    return remaining;
}
