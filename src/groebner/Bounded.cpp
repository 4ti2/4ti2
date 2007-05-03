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

#include "Bounded.h"
#include "Globals.h"

#include "BitSetStream.h"
#include "VectorStream.h"
#include "VectorArrayStream.h"
#include "LatticeBasis.h"
#include "HermiteAlgorithm.h"
#include "RayAlgorithm.h"

//#define DEBUG_4ti2(X) X
#include "Debug.h"

extern "C" {
#include "glpk.h"
}

// TODO: Fix this up.
#ifdef _4ti2_GMP_
  #define DOUBLE(X) X.get_d()
#else
  #define DOUBLE(X) X
#endif

namespace _4ti2_ {
void matrix_bounded(
                const VectorArray& matrix,
                const BitSet& urs,
                BitSet& bounded,
                Vector& grading);
void lattice_unbounded(
                const VectorArray& lattice,
                const BitSet& urs,
                BitSet& unbounded,
                Vector& ray);
bool is_matrix_non_negative(
                const Vector& v,
                const BitSet& urs,
                const BitSet& bounded);
bool is_matrix_non_positive(
                const Vector& v,
                const BitSet& urs,
                const BitSet& bounded);
bool is_lattice_non_negative(
                const Vector& v,
                const BitSet& urs,
                const BitSet& bounded);
bool is_lattice_non_positive(
                const Vector& v,
                const BitSet& urs,
                const BitSet& bounded);
void add_positive_support(
                const Vector& v,
                const BitSet& urs,
                BitSet& bounded,
                Vector& grading);
void add_negative_support(
                const Vector& v,
                const BitSet& urs,
                BitSet& bounded,
                Vector& grading);
void lp_bounded(const VectorArray& matrix,
                const VectorArray& lattice,
                const BitSet& urs,
                BitSet& bounded,
                Vector& grading,
                BitSet& unbounded,
                Vector& ray);
void reconstruct_primal_integer_solution(
                const VectorArray& matrix,
                const BitSet& basic,
                const BitSet& ones,
                Vector& solution);
void reconstruct_dual_integer_solution(
                const VectorArray& matrix,
                const VectorArray& lattice,
                const BitSet& basic,
                const BitSet& ones,
                Vector& solution);
void reconstruct_primal_integer_solution(
                const VectorArray& matrix,
                const BitSet& basic,
                const Vector& rhs,
                Vector& solution);
void compute_ray(
                const VectorArray& lattice,
                const BitSet& bounded,
                const BitSet& unbounded,
                const BitSet& urs);

void load_matrix(LPX* lp, const VectorArray& matrix);
void load_matrix_transpose(LPX* lp, const VectorArray& matrix);

} // namespace _4ti2_;

using namespace _4ti2_;

void
_4ti2_::bounded(const VectorArray& matrix,
                const VectorArray& lattice,
                const BitSet& urs,
                BitSet& bounded,
                Vector& grading,
                BitSet& unbounded,
                Vector& ray)
{
    // If all components are labeled as either bounded, unbounded, or urs, then
    // we are finished.
    if (bounded.count()+unbounded.count()+urs.count() == matrix.get_size())
    { return; }

    // Compute the variables we can infer are bounded from the matrix.
    matrix_bounded(matrix, urs, bounded, grading);
    // If all components are labeled as either bounded, unbounded, or urs, then
    // we are finished.
    if (bounded.count()+unbounded.count()+urs.count() == matrix.get_size())
    { return; }

    // Compute the variables we can infer are unbounded from the lattice.
    lattice_unbounded(lattice, urs, unbounded, ray);
    // If all components are labeled as either bounded, unbounded, or urs, then
    // we are finished.
    if (bounded.count()+unbounded.count()+urs.count() == matrix.get_size())
    { return; }

    // Next we use linear programming to compute the variables that are
    // bounded and unbounded.
    lp_bounded(matrix, lattice, urs, bounded, grading, unbounded, ray);
}

// If a row of the matrix is non-negative or non-positive and all urs
// components are zero, then all non-zero components of the row must be
// bounded.
void
_4ti2_::matrix_bounded(
                const VectorArray& _matrix,
                const BitSet& urs,
                BitSet& bounded,
                Vector& grading)
{
    DEBUG_4ti2(*out << "Checking matrix for bounded components.\n";)
    VectorArray matrix(_matrix);
    // Eliminate the urs variables.
    int rows = upper_triangle(matrix, urs);
    matrix.remove(0, rows);
    DEBUG_4ti2(*out << "Matrix:\n" << matrix << "\n";)
    DEBUG_4ti2(*out << "URS:\n" << urs << "\n";)
    DEBUG_4ti2(*out << "Bounded:\n" << bounded << "\n";)
    DEBUG_4ti2(*out << "Grading:\n" << grading << "\n";)
    while (bounded.count()+urs.count() < bounded.get_size())
    {
        Size count = bounded.count();
        for (Index i = 0; i < matrix.get_number(); ++i)
        {
            if (is_matrix_non_negative(matrix[i], urs, bounded))
            {
                add_positive_support(matrix[i], urs, bounded, grading);
            }
            if (is_matrix_non_positive(matrix[i], urs, bounded))
            {
                add_negative_support(matrix[i], urs, bounded, grading);
            }
        }
        // If nothing changed, then stop.
        if (count == bounded.count()) { return; }
    }
}

// Returns true if the vector v is zero on the urs indices and non-negative on
// the non-bounded indices.
bool
_4ti2_::is_matrix_non_negative(
                const Vector& v,
                const BitSet& urs,
                const BitSet& bounded)
{
    bool result = false;
    for (Index i = 0; i < v.get_size(); ++i)
    {
        if (urs[i] && v[i] != 0) { return false; }
        if (!bounded[i])
        {
            if (v[i] < 0) { return false; }
            if (v[i] > 0) { result = true; }
        }
    }
    return result;
}

// Returns true if the vector v is zero on the urs indices and non-positive on
// the non-bounded indices.
bool
_4ti2_::is_matrix_non_positive(
                const Vector& v,
                const BitSet& urs,
                const BitSet& bounded)
{
    bool result = false;
    for (Index i = 0; i < v.get_size(); ++i)
    {
        if (urs[i] && v[i] != 0) { return false; }
        if (!bounded[i])
        {
            if (v[i] > 0) { return false; }
            if (v[i] < 0) { result = true; }
        }
    }
    return result;
}

// If a vector in the lattice basis is non-negative or non-positive
// excluding the urs components, then all non-zero components of the row
// must be unbounded.
void
_4ti2_::lattice_unbounded(
                const VectorArray& lattice,
                const BitSet& urs,
                BitSet& unbounded,
                Vector& ray)
{
    while (unbounded.count()+urs.count() < unbounded.get_size())
    {
        Size count = unbounded.count();
        for (Index i = 0; i < lattice.get_number(); ++i)
        {
            if (is_lattice_non_negative(lattice[i], urs, unbounded))
            {
                add_positive_support(lattice[i], urs, unbounded, ray);
            }
            if (is_lattice_non_positive(lattice[i], urs, unbounded))
            {
                add_negative_support(lattice[i], urs, unbounded, ray);
            }
        }
        // If nothing changed, then stop.
        if (count == unbounded.count()) { break; }
    }
}

bool
_4ti2_::is_lattice_non_negative(
                const Vector& v,
                const BitSet& urs,
                const BitSet& unbounded)
{
    bool result = false;
    for (Index i = 0; i < v.get_size(); ++i)
    {
        if (!urs[i] && !unbounded[i])
        {
            if (v[i] < 0) { return false; }
            if (v[i] > 0) { result = true; }
        }
    }
    return result;
}

bool
_4ti2_::is_lattice_non_positive(
                const Vector& v,
                const BitSet& urs,
                const BitSet& unbounded)
{
    bool result = false;
    for (Index i = 0; i < v.get_size(); ++i)
    {
        if (!urs[i] && !unbounded[i])
        {
            if (v[i] > 0) { return false; }
            if (v[i] < 0) { result = true; }
        }
    }
    return result;
}

void
_4ti2_::add_positive_support(
                const Vector& v,
                const BitSet& urs,
                BitSet& bounded,
                Vector& grading)
{
    IntegerType factor = 1;
    for (Index i = 0; i < v.get_size(); ++i)
    {
        if (!urs[i])
        {
            if (v[i] > 0) { bounded.set(i); }
            else if (v[i] < 0)
            {
                assert(bounded[i]);
                assert(grading[i]);
                IntegerType ratio = (-v[i])/grading[i] + 1;
                if (ratio > factor) { factor = ratio; }
            }
        }
    }
    Vector::add(grading, factor, v, 1, grading);
    grading.normalise();
}

void
_4ti2_::add_negative_support(
                const Vector& v,
                const BitSet& urs,
                BitSet& bounded,
                Vector& grading)
{
    IntegerType factor = 1;
    for (Index i = 0; i < v.get_size(); ++i)
    {
        if (!urs[i])
        {
            if (v[i] < 0) { bounded.set(i); }
            else if (v[i] > 0)
            {
                assert(bounded[i]);
                IntegerType ratio = v[i]/grading[i] + 1;
                if (ratio > factor) { factor = ratio; }
            }
        }
    }
    Vector::sub(grading, factor, v, 1, grading);
    grading.normalise();
}

bool
_4ti2_::lp_feasible(
                const VectorArray& lattice,
                const Vector& rhs)
{
    assert(lattice.get_size() == rhs.get_size());
    if (lattice.get_number() == 0) { return rhs.is_non_negative(); }

    int m = lattice.get_size();
    int n = lattice.get_number();
    DEBUG_4ti2(*out << "M = " << m << " N = " << n << "\n";)
    DEBUG_4ti2(*out << "lattice:\n" << lattice << "\n";)

    LPX *lp = lpx_create_prob();
    lpx_set_int_parm(lp,LPX_K_MSGLEV,0);
    DEBUG_4ti2(lpx_set_int_parm(lp,LPX_K_MSGLEV,2);)
    lpx_set_obj_dir(lp, LPX_MAX);

    lpx_add_rows(lp, m);
    for (int i = 1; i <= m; ++i)
    {
        lpx_set_row_bnds(lp, i, LPX_UP, 0.0, DOUBLE(rhs[i-1]));
    }

    lpx_add_cols(lp, n);
    for (int j = 1; j <= n; ++j)
    {
        lpx_set_col_bnds(lp, j, LPX_FR, 0.0, 0.0);
        lpx_set_obj_coef(lp, j, 0.0);
    }

    load_matrix_transpose(lp, lattice);
    //lpx_print_prob(lp, "model");

    lpx_simplex(lp);

    // TODO: Check for other exit statuses!
    bool feasible = true;
    switch(lpx_get_status(lp))
    {
        case LPX_NOFEAS:
        case LPX_INFEAS:
            feasible = false;
            break;
    };

    DEBUG_4ti2(
    if (feasible)
    {
        *out << "Primal Solution:\n";
        for (int j = 0; j < lattice.get_size(); ++j)
        {
            double acc = 0;
            for (int i = 0; i < lattice.get_number(); ++i)
            {
                acc += lpx_get_col_prim(lp, i+1)*DOUBLE(lattice[i][j]);
            }
            *out << " " << rhs[j] - acc;
        }
        *out << "\n";
    })

    lpx_delete_prob(lp);
    return feasible;
}

bool
_4ti2_::ip_feasible(
                const VectorArray& lattice,
                const Vector& rhs)
{
    assert(lattice.get_size() == rhs.get_size());
    if (lattice.get_number() == 0) { return rhs.is_non_negative(); }

    int m = lattice.get_size();
    int n = lattice.get_number();
    DEBUG_4ti2(*out << "M = " << m << " N = " << n << "\n";)
    DEBUG_4ti2(*out << "lattice:\n" << lattice << "\n";)

    LPX *lp = lpx_create_prob();
    lpx_set_int_parm(lp,LPX_K_MSGLEV,0);
    DEBUG_4ti2(lpx_set_int_parm(lp,LPX_K_MSGLEV,2);)
    lpx_set_obj_dir(lp, LPX_MAX);

    lpx_add_rows(lp, m);
    for (int i = 1; i <= m; ++i)
    {
        lpx_set_row_bnds(lp, i, LPX_UP, 0.0, DOUBLE(rhs[i-1]));
    }

    lpx_add_cols(lp, n);
    for (int j = 1; j <= n; ++j)
    {
        lpx_set_col_bnds(lp, j, LPX_FR, 0.0, 0.0);
        lpx_set_obj_coef(lp, j, 0.0);
    }

    load_matrix_transpose(lp, lattice);
    //lpx_print_prob(lp, "model");

    lpx_simplex(lp);

    // TODO: Check for other exit statuses!
    switch(lpx_get_status(lp))
    {
        case LPX_NOFEAS:
        case LPX_INFEAS:
            lpx_delete_prob(lp);
            return false;
            break;
    };

    DEBUG_4ti2(
    *out << "Linear Primal Solution:\n";
    for (int j = 0; j < lattice.get_size(); ++j)
    {
        double acc = 0;
        for (int i = 0; i < lattice.get_number(); ++i)
        {
            acc += lpx_get_col_prim(lp, i+1)*DOUBLE(lattice[i][j]);
        }
        *out << " " << rhs[j] - acc;
    }
    *out << "\n";
    )

    lpx_set_class(lp, LPX_MIP);
    for (int i = 1; i <= n; ++i)
    {
        lpx_set_col_kind(lp, i, LPX_IV);
    }

    lpx_integer(lp);

    // TODO: Check for other exit statuses!
    bool feasible = true;
    switch(lpx_mip_status(lp))
    {
        case LPX_I_NOFEAS:
            feasible = false;
            break;
    };

    if (feasible)
    {
        DEBUG_4ti2(
        *out << "Integer Primal Solution:\n";
        for (int j = 0; j < lattice.get_size(); ++j)
        {
            double acc = 0;
            for (int i = 0; i < lattice.get_number(); ++i)
            {
                acc += lpx_mip_col_val(lp, i+1)*DOUBLE(lattice[i][j]);
            }
            *out << " " << rhs[j] - acc;
        }
        *out << "\n";
        )
    }

    lpx_delete_prob(lp);
    return feasible;
}

void
_4ti2_::lp_bounded(
                const VectorArray& matrix,
                const VectorArray& tmp_lattice,
                const BitSet& urs,
                BitSet& bounded,
                Vector& grading,
                BitSet& unbounded,
                Vector& ray)
{
    DEBUG_4ti2(*out << "lp_bounded.\n";)
    VectorArray lattice(tmp_lattice);
    BitSet rs(urs);
    rs.set_complement();
    // TODO: We should not use the hermite normal form since sometimes this
    // introduces numerical problems.
    int rows = hermite(lattice, rs);
    lattice.remove(rows, lattice.get_number());

    LPX *lp = lpx_create_prob();
    lpx_set_int_parm(lp,LPX_K_MSGLEV,0);
    DEBUG_4ti2(lpx_set_int_parm(lp,LPX_K_MSGLEV,2);)
    lpx_set_obj_dir(lp, LPX_MAX);

    // Set RHS = 0.
    lpx_add_rows(lp, lattice.get_number());
    for (int i = 1; i <= lattice.get_number(); ++i)
    {
        lpx_set_row_bnds(lp, i, LPX_FX, 0.0, 0.0);
    }

    // Set 0 <= x <= 1.
    // Set c = 0
    lpx_add_cols(lp, lattice.get_size());
    for (int i = 1; i <= lattice.get_size(); ++i)
    {
        if (!urs[i-1]) { lpx_set_col_bnds(lp, i, LPX_DB, 0.0, 1.0); }
        else { lpx_set_col_bnds(lp, i, LPX_FX, 0.0, 0.0); }
        lpx_set_obj_coef(lp, i, 0.0);
    }

    // Input the lattice.
    int *ai = new int[lattice.get_size()*lattice.get_number()+1];
    int *aj = new int[lattice.get_size()*lattice.get_number()+1];
    double *ar = new double[lattice.get_size()*lattice.get_number()+1];
    int index = 1;
    for (int i = 0; i < lattice.get_number(); ++i)
    {
        for (int j = 0; j < lattice.get_size(); ++j)
        {
            if (!urs[j] && lattice[i][j] != 0)
            {
                ai[index] = i+1;
                aj[index] = j+1;
                ar[index] = DOUBLE(lattice[i][j]);
                ++index;
            }
        }
    }
    lpx_load_matrix(lp, index-1, ai, aj, ar);
    delete [] ai; delete [] aj; delete [] ar;

    while (urs.count()+bounded.count()+unbounded.count() < matrix.get_size())
    {
        // The set of components that we do not know whether they are urs,
        // bounded or unbounded.
        BitSet unknown(bounded.get_size(), true);
        unknown.set_difference(urs);
        unknown.set_difference(bounded);
        unknown.set_difference(unbounded);
        DEBUG_4ti2(*out << "Unknown:\n" << unknown << "\n";)
        DEBUG_4ti2(*out << "URS:\n" << urs << "\n";)
        DEBUG_4ti2(*out << "Bounded:\n" << bounded << "\n";)
        DEBUG_4ti2(*out << "Unbounded:\n" << unbounded << "\n";)

        // Run the simplex algorithm.
        for (int i = 1; i <= lattice.get_size(); ++i)
        {
            if (unknown[i-1])
            {
                lpx_set_obj_coef(lp, i, 1.0);
                lpx_set_col_bnds(lp, i, LPX_DB, 0.0, 1.0);
            }
            else
            {
                lpx_set_obj_coef(lp, i, 0.0);
                lpx_set_col_bnds(lp, i, LPX_LO, 0.0, 0.0);
            }
        }
        lpx_adv_basis(lp);
        lpx_simplex(lp);

        DEBUG_4ti2(
        *out << "Number of simplex iterations: " <<
               lpx_get_int_parm(lp, LPX_K_ITCNT) << "\n";
        *out << "Primal Solution:\n";
        for (int j = 1; j <= lattice.get_size(); ++j)
        {
            *out << " " << lpx_get_col_prim(lp, j);
        }
        *out << "\n";
        *out << "Dual Solution:\n";
        for (int i = 1; i <= lattice.get_number(); ++i)
        {
            *out << " " << lpx_get_row_dual(lp, i);
        }
        *out << "\n";
        *out << "Objective = " << lpx_get_obj_val(lp) << "\n";
        )

        // If the solution is 0, then all the unknown variables are unbounded.
        BitSet basic(lattice.get_size());
        BitSet ones(lattice.get_size());
        for (int j = 1; j <= lattice.get_size(); ++j)
        {
            switch(lpx_get_col_stat(lp, j))
            {
                case LPX_BS:
                    basic.set(j-1);
                    break;
                case LPX_NU:
                    ones.set(j-1);
                    break;
                case LPX_NL:
                case LPX_NS:
                    break;
                case LPX_NF:
                    std::cerr<<"Received LPX_NF for component "<<j-1<<".\n";
                default:
                    std::cerr<<"LP solver unexpected output error.\n";
                    exit(1);
                    break;
            };
        }
        DEBUG_4ti2(*out << "Basic variables:\n" << basic << "\n";)
        DEBUG_4ti2(*out << "Upper bound variables:\n" << ones << "\n";)

        Vector solution(lattice.get_size(),0);
        if (lpx_get_obj_val(lp) < 0.5) // It should be at least 1 if not zero.
        {
            DEBUG_4ti2(*out << "Unknown are all unbounded.\n";)
            reconstruct_dual_integer_solution(matrix,lattice,basic,unknown,solution);
            add_positive_support(solution, urs, unbounded, ray);
        }
        else
        {
            reconstruct_primal_integer_solution(lattice,basic,ones,solution);
            add_positive_support(solution, urs, bounded, grading);
        }
    }
    lpx_delete_prob(lp);
}

// TODO: This function is quite ugly and should be rewritten.
bool
_4ti2_::bounded(const VectorArray& matrix,
                const VectorArray& lattice,
                const BitSet& urs,
                const VectorArray& cost,
                const BitSet& bnd,
                const Vector& grading,
                const BitSet& unbnd,
                const Vector& ray,
                BitSet& cost_unbnd)
{
    DEBUG_4ti2(*out << "Determining if cost is bounded or unbounded...\n";)
    assert(matrix.get_size() == lattice.get_size());
    assert(matrix.get_size() == urs.get_size());
    assert(matrix.get_size() == cost.get_size());
    assert(matrix.get_size() == bnd.get_size());
    assert(matrix.get_size() == unbnd.get_size());

    // If every component is bnd, then there is nothing to do.
    if (unbnd.empty()) { cost_unbnd.zero(); return true; }
    // If there are no cost vectors, then there is nothing to do.
    if (cost.get_number() == 0) { cost_unbnd = unbnd; return true; }

    // Extend the matrix and the lattice to allow for the cost component.
    VectorArray full_matrix(matrix.get_number(), matrix.get_size()+1, 0);
    VectorArray::lift(matrix, 0, matrix.get_size(), full_matrix);
    DEBUG_4ti2(*out << "Extended Matrix:\n" << full_matrix << "\n";)
    VectorArray full_lattice(lattice.get_number(), lattice.get_size()+1, 0);
    VectorArray::lift(lattice, 0, lattice.get_size(), full_lattice);
    DEBUG_4ti2(*out << "Extended Lattice:\n" << full_lattice << "\n";)
    VectorArray full_cost(cost.get_number(), cost.get_size()+1, 0);
    VectorArray::lift(cost, 0, cost.get_size(), full_cost);
    DEBUG_4ti2(*out << "Extended Cost:\n" << full_cost << "\n";)

    // Extend the urs, bnd, and unbnd sets.
    BitSet full_urs(urs.get_size()+1);
    for (int i = 0; i < urs.get_size(); ++i)
    { if (urs[i]) { full_urs.set(i); } }
    DEBUG_4ti2(*out << "Extended URS:\n" << full_urs << "\n";)
    BitSet full_bnd(bnd.get_size()+1);
    for (int i = 0; i < bnd.get_size(); ++i)
    { if (bnd[i]) { full_bnd.set(i); } }
    DEBUG_4ti2(*out << "Extended BND:\n" << full_bnd << "\n";)
    BitSet full_unbnd(unbnd.get_size()+1);
    DEBUG_4ti2(*out << "Extended UNBND:\n" << full_unbnd << "\n";)

    // Extend the grading.
    Vector full_grading(grading.get_size()+1,0);
    Vector::lift(grading, 0, grading.get_size(), full_grading);
    DEBUG_4ti2(*out << "Extended Grading:\n" << full_grading << "\n";)

    int col = full_matrix.get_size()-1;
    BitSet col_set(full_matrix.get_size());
    col_set.set(col);
    int m = 0;
    while (m < cost.get_number())
    {
        // Insert the cost vector into the matrix.
        Vector  mtmp(full_matrix.get_size(),0);
        Vector::lift(cost[m], 0, cost.get_size(), mtmp);
        mtmp[col] = 1;
        full_matrix.insert(mtmp, 0);
        DEBUG_4ti2(*out << "Full Matrix:\n" << full_matrix << "\n";)

        // Add the cost component to the lattice.
        Vector ltmp(full_lattice.get_number(),0);
        VectorArray::dot(full_lattice, full_cost[m], ltmp);
        if (!ltmp.is_zero())
        {
            for (int i = 0; i < full_lattice.get_number(); ++i)
            {   full_lattice[i][col] = -ltmp[i]; }
            DEBUG_4ti2(*out << "Full Lattice:\n" << full_lattice << "\n";)
    
            // Check whether the cost component is bounded.
            full_unbnd.zero();
            Vector full_ray(ray.get_size()+1,0);
            full_grading[grading.get_size()]=0;
            bounded(full_matrix, full_lattice, full_urs,
                            full_bnd, full_grading,
                            full_unbnd, full_ray);
            DEBUG_4ti2(*out << "Grading:\n" << full_grading << "\n";)
            DEBUG_4ti2(if (!full_unbnd.empty()) {*out<<"Cost Ray:\n"<<ray<<"\n";})
            DEBUG_4ti2(*out << "EXT BND:\n" << full_bnd << "\n";)
            DEBUG_4ti2(*out << "EXT UNBND:\n" << full_unbnd << "\n";)
            if (!full_bnd[col]) { return false; }
    
            if (m+1 != cost.get_number())
            {
                // Zero the cost component.
                full_matrix[0][col] = 0;
                DEBUG_4ti2(*out << "Zeroed Full Matrix:\n" << full_matrix << "\n";)
                // Eliminate the cost component.
                upper_triangle(full_lattice, col_set);
                full_lattice.remove(0);
                DEBUG_4ti2(*out << "Eliminated Full Lattice:\n" << full_lattice << "\n";)
            }
        }
        ++m;
    }
    cost_unbnd.zero();
    for (int i = 0; i < cost_unbnd.get_size(); ++i)
    { if (full_unbnd[i]) { cost_unbnd.set(i); } }
    DEBUG_4ti2(*out << "Cost Unbounded:\n" << cost_unbnd << "\n";)
    return true;
}

// TODO: We should have a version of bounded_projection which uses linear
// programming instead of ray computation because of numerical stability.
void
_4ti2_::bounded_projection(
                const VectorArray& _matrix,
                const VectorArray& lattice,
                const BitSet& urs,
                const Vector& rhs,
                BitSet& proj)
{
#if 1
    VectorArray rays(lattice);
    VectorArray matrix(_matrix);
    VectorArray subspace(0,rays.get_size());
    BitSet rs(urs);
    rs.set_complement();
    // TODO: We suppress the output of the following ray calculation.
    // Is this what we really want?
    std::ostream* tmp_out = out;
    out = new std::ofstream;
    RayAlgorithm ray_algorithm;
    proj = ray_algorithm.compute(matrix, rays, subspace, rs);
    rays.clear();
    delete out; out = tmp_out;
    DEBUG_4ti2(*out << "Projection: " << proj << "\n";)
#endif
#if 0
    int dim = lattice.get_size();
    Vector tmp(dim);
    if (Globals::norm == 2)
    {
        lp_weight_l2(lattice, urs, rhs, tmp);
    }
    else
    {
        lp_weight_l1(lattice, urs, rhs, tmp);
    }
    proj.zero();
    for (int i = 0; i < dim; ++i) { if (tmp[i] == 0) { proj.set(i); } }
    std::cout << "FIN:\n" << proj << "\n";

    assert(BitSet::set_disjoint(proj, urs));
#endif
}

void
_4ti2_::lp_weight_l2(
                const VectorArray& tmp_lattice,
                const BitSet& urs,
                const Vector& rhs,
                Vector& weight)
{
    VectorArray basis(0, tmp_lattice.get_size());
    lattice_basis(tmp_lattice, basis);
    int rows = upper_triangle(basis, urs);
    DEBUG_4ti2(*out << "Basis\n" << basis << "\n";)
    basis.remove(0,rows);
    VectorArray matrix(0, tmp_lattice.get_size());
    lattice_basis(basis, matrix);
    DEBUG_4ti2(*out << "Basis\n" << basis << "\n";)
    DEBUG_4ti2(*out << "Matrix\n" << matrix << "\n";)

    BitSet rs(urs);
    rs.set_complement();
    VectorArray subspace(0,basis.get_size());
    RayAlgorithm algorithm;
    algorithm.compute(matrix, basis, subspace, rs);
    DEBUG_4ti2(*out << "Rays\n" << basis << "\n";)

    if (basis.get_number() == 0) { return; }
    int index = 0;
    RationalType dot = Vector::dot(rhs, basis[0]);
    RationalType ratio = 0;
    for (int j = 0; j < basis.get_size(); ++j)
    {
        ratio += basis[0][j]*(basis[0][j]/dot);
    }
    RationalType min = ratio;
    for (int i = 1; i < basis.get_number(); ++i)
    {
        ratio = 0;
        dot = Vector::dot(rhs, basis[i]);
        for (int j = 0; j < basis.get_size(); ++j)
        {
            ratio += basis[i][j]*(basis[i][j]/dot);
        }
        if (ratio > min) { index = i; min = ratio; }
    }
    DEBUG_4ti2(*out << "Index is " << index << " ratio is " << min << "\n";)
    weight = basis[index];
}

void
_4ti2_::lp_weight_l1(
                const VectorArray& tmp_lattice,
                const BitSet& urs,
                const Vector& rhs,
                Vector& weight)
{
    DEBUG_4ti2(*out << "Computing Weight vector...\n";)
    assert(tmp_lattice.get_size() == rhs.get_size());
    assert(tmp_lattice.get_size() == urs.get_size());
    assert(weight.get_size() == rhs.get_size());

    VectorArray lattice(tmp_lattice);
    //BitSet rs(urs);
    //rs.set_complement();
    //int rows = hermite(lattice, rs);
    //lattice.remove(rows, lattice.get_number());
    lattice.insert(Vector(lattice.get_size(), 1));

    LPX *lp = lpx_create_prob();
    lpx_set_int_parm(lp,LPX_K_MSGLEV,0);
    DEBUG_4ti2(lpx_set_int_parm(lp,LPX_K_MSGLEV,2);)
    lpx_set_obj_dir(lp, LPX_MIN);

    // Set RHS = 0 except for the last row which is 1 (the l1 norm).
    lpx_add_rows(lp, lattice.get_number());
    for (int i = 1; i <= lattice.get_number()-1; ++i)
    {
        lpx_set_row_bnds(lp, i, LPX_FX, 0.0, 0.0);
    }
    lpx_set_row_bnds(lp, lattice.get_number(), LPX_FX, 1.0, 1.0);

    // Set c = rhs
    lpx_add_cols(lp, lattice.get_size());
    for (int i = 1; i <= lattice.get_size(); ++i)
    {
        if (!urs[i-1]) { lpx_set_col_bnds(lp, i, LPX_LO, 0.0, 0.0); }
        else { lpx_set_col_bnds(lp, i, LPX_FX, 0.0, 0.0); }
        lpx_set_obj_coef(lp, i, DOUBLE(rhs[i-1]));
    }

    // Input the lattice.
    int *ai = new int[lattice.get_size()*lattice.get_number()+1];
    int *aj = new int[lattice.get_size()*lattice.get_number()+1];
    double *ar = new double[lattice.get_size()*lattice.get_number()+1];
    int index = 1;
    for (int i = 0; i < lattice.get_number(); ++i)
    {
        for (int j = 0; j < lattice.get_size(); ++j)
        {
            if (!urs[j] && lattice[i][j] != 0)
            {
                ai[index] = i+1;
                aj[index] = j+1;
                ar[index] = DOUBLE(lattice[i][j]);
                ++index;
            }
        }
    }
    lpx_load_matrix(lp, index-1, ai, aj, ar);
    delete [] ai; delete [] aj; delete [] ar;

    lpx_simplex(lp);

    switch(lpx_get_status(lp))
    {
        case LPX_NOFEAS:
        case LPX_INFEAS:
            return;
            break;
    };
    
    DEBUG_4ti2(
    *out << "Number of simplex iterations: " <<
               lpx_get_int_parm(lp, LPX_K_ITCNT) << "\n";
    *out << "Primal Solution:\n";
    for (int j = 1; j <= lattice.get_size(); ++j)
    {
        *out << " " << lpx_get_col_prim(lp, j);
    }
    *out << "\n";
    *out << "Dual Solution:\n";
    for (int i = 1; i <= lattice.get_number(); ++i)
    {
        *out << " " << lpx_get_row_dual(lp, i);
    }
    *out << "\n";
    *out << "Objective = " << lpx_get_obj_val(lp) << "\n";
    )

    BitSet basic(lattice.get_size());
    BitSet ones(lattice.get_size());
    for (int j = 1; j <= lattice.get_size(); ++j)
    {
         switch(lpx_get_col_stat(lp, j))
         {
             case LPX_BS:
                 basic.set(j-1);
                 break;
             case LPX_NU:
                 ones.set(j-1);
                 break;
             case LPX_NL:
             case LPX_NS:
                 break;
             case LPX_NF:
                 std::cerr<<"Received LPX_NF for component "<<j-1<<".\n";
             default:
                 std::cerr<<"LP solver unexpected output error.\n";
                 exit(1);
                 break;
         };
    }
    DEBUG_4ti2(
    *out << "Basic variables:\n" << basic << "\n";
    *out << "Upper bound variables:\n" << ones << "\n";
    )

    Vector sol(lattice.get_number(),0);
    sol[lattice.get_number()-1] = 1;
    reconstruct_primal_integer_solution(lattice,basic,sol,weight);
    DEBUG_4ti2(*out << "Weight:\n" << weight << "\n";)

    lpx_delete_prob(lp);
}

int
_4ti2_::lp_solve(
                const VectorArray& matrix,
                const Vector& rhs,
                const Vector& cost,
                const BitSet& urs,
                BitSet& basic,
                RationalType& objective)
{
    assert(matrix.get_number() == rhs.get_size());
    assert(matrix.get_size() == cost.get_size());
    assert(matrix.get_size() == urs.get_size());
    assert(matrix.get_size() == basic.get_size());

    LPX *lp = lpx_create_prob();
    lpx_set_int_parm(lp,LPX_K_MSGLEV,0);
    DEBUG_4ti2(lpx_set_int_parm(lp,LPX_K_MSGLEV,2);)
    lpx_set_obj_dir(lp, LPX_MIN);

    int m = matrix.get_number();
    int n = matrix.get_size();

    // Set rhs.
    lpx_add_rows(lp, m);
    for (int i = 1; i <= m; ++i)
    {
        lpx_set_row_bnds(lp, i, LPX_FX, DOUBLE(rhs[i-1]), 0.0);
    }

    // Set cost.
    lpx_add_cols(lp, n);
    for (int i = 1; i <= n; ++i)
    {
        lpx_set_obj_coef(lp, i, DOUBLE(cost[i-1]));
        if (urs[i-1]) { lpx_set_col_bnds(lp, i, LPX_FR, 0.0, 0.0); }
        else { lpx_set_col_bnds(lp, i, LPX_LO, 0.0, 0.0); }
    }

    // Input the matrix.
    load_matrix(lp, matrix);

    lpx_simplex(lp);

    switch(lpx_get_status(lp))
    {
        case LPX_OPT:
             DEBUG_4ti2(*out << "Optimal solution found.\n";)
             break;
        case LPX_NOFEAS:
        case LPX_INFEAS:
            DEBUG_4ti2(*out << "Problem is infeasible.\n";)
            return -1;
            break;
        case LPX_UNBND:
            DEBUG_4ti2(*out << "Problem is unbounded.\n";)
            return 1;
            break;
        case LPX_UNDEF:
        case LPX_FEAS:
        default:
            std::cerr << "Software Error: Received unexpected lp solver output.\n";
            exit(1);
    };

    DEBUG_4ti2(
    *out << "Number of simplex iterations: " <<
               lpx_get_int_parm(lp, LPX_K_ITCNT) << "\n";
    *out << "Primal Solution:\n";
    for (int j = 1; j <= matrix.get_size(); ++j)
    {
        *out << " " << lpx_get_col_prim(lp, j);
    }
    *out << "\n";
    *out << "Dual Solution:\n";
    for (int i = 1; i <= matrix.get_number(); ++i)
    {
        *out << " " << lpx_get_row_dual(lp, i);
    }
    *out << "\n";
    *out << "Objective = " << lpx_get_obj_val(lp) << "\n";
    )

    objective = lpx_get_obj_val(lp);

    // Find the basic variables.
    for (int j = 1; j <= n; ++j)
    {
         switch(lpx_get_col_stat(lp, j))
         {
             case LPX_BS:
                 basic.set(j-1);
                 break;
             case LPX_NU:
             case LPX_NL:
             case LPX_NS:
             case LPX_NF:
                 break;
             default:
                 std::cerr<<"LP solver unexpected output error.\n";
                 exit(1);
                 break;
         };
    }
    DEBUG_4ti2(*out << "Basic variables:\n" << basic << "\n";)

    lpx_delete_prob(lp);

    return 0;
}

void
_4ti2_::reconstruct_primal_integer_solution(
                const VectorArray& matrix,
                const BitSet& basic,
                const Vector& rhs,
                Vector& solution)
{
    assert(matrix.get_size() == basic.get_size());
    assert(matrix.get_number() == rhs.get_size());
    assert(matrix.get_size() == solution.get_size());
    VectorArray basic_matrix(matrix.get_number(), basic.count(), 0);
    VectorArray::project(matrix, basic, basic_matrix);
    DEBUG_4ti2(*out << "Matrix:\n" << basic_matrix << "\n";)

    Vector basic_solution(basic.count());
    if (solve(basic_matrix, rhs, basic_solution) == 0)
    {
        std::cerr << "Software Error: Unable to reconstruct primal solution.\n";
        exit(1);
    }

    solution.mul(0);
    Vector::lift(basic_solution, basic, solution);
    DEBUG_4ti2(*out << "Vector:\n" << solution << "\n";)
}

void
_4ti2_::reconstruct_primal_integer_solution(
                const VectorArray& matrix,
                const BitSet& basic,
                const BitSet& ones,
                Vector& solution)
{
    DEBUG_4ti2(*out << "Reconstructing primal integer solution.\n";)
    VectorArray basic_matrix(matrix.get_number(), basic.count(), 0);
    VectorArray::project(matrix, basic, basic_matrix);
    Vector rhs(matrix.get_number(),0);
    int index = 0;
    for (int j = 0; j < matrix.get_size(); ++j)
    {
        if (ones[j])
        {
            for (int i = 0; i < matrix.get_number(); ++i)
            {
                rhs[i] -= matrix[i][j];
            }
        }
    }
    DEBUG_4ti2(*out << "Matrix:\n" << basic_matrix << "\n";)
    DEBUG_4ti2(*out << "RHS:\n" << rhs << "\n";)

    Vector basic_solution(basic.count());
    IntegerType factor = solve(basic_matrix, rhs, basic_solution);
    DEBUG_4ti2(*out << "Basic Solution:\n" << basic_solution << "\n";)
    DEBUG_4ti2(*out << "Basic Solution Factor:\n" << factor << "\n";)

    if (factor == 0)
    {
        std::cerr << "Software Error: Unable to reconstruct primal solution.\n";
        exit(1);
    }

    Vector::lift(basic_solution, basic, solution);
    index = 0;
    for (int j = 0; j < solution.get_size(); ++j)
    {
        if (ones[j])
        {
            solution[j] = factor;
        }
    }
    DEBUG_4ti2(*out << "Vector:\n" << solution << "\n";)

    // Check vector is in null space of matrix.
    Vector check(matrix.get_number());
    VectorArray::dot(matrix, solution, check);
    Vector zero(matrix.get_number(),0);
    DEBUG_4ti2(*out << "Check Vector:\n" << check << "\n";)
    if (check != zero)
    {
        *out << "ERROR: Integer Solution not in matrix.\n";
        exit(1);
    }
    DEBUG_4ti2(*out << "Integer Solution:\n" << solution << "\n";)
}

void
_4ti2_::reconstruct_dual_integer_solution(
                const VectorArray& matrix,
                const VectorArray& lattice,
                const BitSet& basic,
                const BitSet& ones,
                Vector& solution)
{
    VectorArray temp(basic.count(), lattice.get_number()+1, 0);
    int index = 0;
    for (int j = 0; j < lattice.get_size(); ++j)
    {
        if (basic[j])
        {
            for (int i = 0; i < lattice.get_number(); ++i)
            {
                temp[index][i] = lattice[i][j];
            }
            if (ones[j])
            {
                temp[index][lattice.get_number()] = -1;
            }
            ++index;
        }
    }
    DEBUG_4ti2(*out << "Matrix:\n" << temp << "\n";)

    VectorArray basis(0, lattice.get_number()+1);
    lattice_basis(temp, basis);
    // TODO: There can be more than one vector in the basis.
    //assert(basis.get_number() == 1);
    //*out << "Basis:\n" << basis << "\n";
    index = 0;
    Vector dual(lattice.get_number());
    for (int j = 0; j < lattice.get_number(); ++j)
    {
        dual[j] = basis[0][j];
    }
    if (basis[0][lattice.get_number()] < 0) { dual.mul(-1); }
    DEBUG_4ti2(*out << "Integer Dual Solution:\n" << dual << "\n";)
    DEBUG_4ti2(*out << "Basis[0]:\n" << basis[0] << "\n";)
    VectorArray trans(lattice.get_size(), lattice.get_number());
    VectorArray::transpose(lattice, trans);
    VectorArray::dot(trans, dual, solution);
    DEBUG_4ti2(*out << "Dual Solution:\n" << solution << "\n";)
}

void
_4ti2_::compute_ray(
                const VectorArray& tmp_lattice,
                const BitSet& bounded,
                const BitSet& unbounded,
                const BitSet& urs)
{
    *out << "Compute Rays.\n";
    *out << "Unbounded:\n" << unbounded << "\n";
    VectorArray lattice(tmp_lattice);

    int rows = upper_triangle(lattice, bounded);
    lattice.remove(0,rows);
    if (lattice.get_number() == 0) { return; }

    int m = lattice.get_size();
    int n = lattice.get_number();
    DEBUG_4ti2(*out << "M = " << m << " N = " << n << "\n";)
    DEBUG_4ti2(*out << "lattice:\n" << lattice << "\n";)

    LPX *lp = lpx_create_prob();
    lpx_set_int_parm(lp,LPX_K_MSGLEV,0);
    DEBUG_4ti2(lpx_set_int_parm(lp,LPX_K_MSGLEV,2);)
    lpx_set_obj_dir(lp, LPX_MAX);

    lpx_add_rows(lp, m);
    for (int i = 1; i <= m; ++i)
    {
        if (unbounded[i-1])
        {
           lpx_set_row_bnds(lp, i, LPX_LO, 1.0, 0.0);
        }
        else
        {
           lpx_set_row_bnds(lp, i, LPX_FR, 0.0, 0.0);
        }
    }

    lpx_add_cols(lp, n);
    for (int j = 1; j <= n; ++j)
    {
        lpx_set_col_bnds(lp, j, LPX_FR, 0.0, 0.0);
        lpx_set_obj_coef(lp, j, 0.0);
    }

    // Input the (transpose) lattice.
    load_matrix_transpose(lp, lattice);

    //lpx_print_prob(lp, "model");

    lpx_simplex(lp);

    // TODO: Check for other exit statuses!
    bool feasible = true;
    switch(lpx_get_status(lp))
    {
        case LPX_NOFEAS:
        case LPX_INFEAS:
            feasible = false;
            lpx_delete_prob(lp);
            *out << "Not feasible.\n";
            return;
            break;
    };

    DEBUG_4ti2(
    if (feasible)
    {
        *out << "Primal Solution:\n";
        for (int j = 0; j < lattice.get_size(); ++j)
        {
            double acc = 0;
            for (int i = 0; i < lattice.get_number(); ++i)
            {
                acc += lpx_get_col_prim(lp, i+1)*DOUBLE(lattice[i][j]);
            }
            *out << " " << acc;
        }
        *out << "\n";
    })

    lpx_set_class(lp, LPX_MIP);
    for (int i = 1; i <= n; ++i)
    {
        lpx_set_col_kind(lp, i, LPX_IV);
    }

    lpx_integer(lp);

    // TODO: Check for other exit statuses!
    feasible = true;
    switch(lpx_mip_status(lp))
    {
        case LPX_I_NOFEAS:
            feasible = false;
            break;
    };

    if (feasible)
    {
        DEBUG_4ti2(
        *out << "Integer Primal Solution:\n";
        for (int j = 0; j < lattice.get_size(); ++j)
        {
            double acc = 0;
            for (int i = 0; i < lattice.get_number(); ++i)
            {
                acc += lpx_mip_col_val(lp, i+1)*DOUBLE(lattice[i][j]);
            }
            *out << " " << acc;
        }
        *out << "\n";
        )
    }

    lpx_delete_prob(lp);
    return;
}

void
_4ti2_::load_matrix(LPX* lp, const VectorArray& matrix)
{
    // Input the matrix.
    int *ai = new int[matrix.get_size()*matrix.get_number()+1];
    int *aj = new int[matrix.get_size()*matrix.get_number()+1];
    double *ar = new double[matrix.get_size()*matrix.get_number()+1];
    int index = 1;
    for (int i = 0; i < matrix.get_number(); ++i)
    {
        for (int j = 0; j < matrix.get_size(); ++j)
        {
            if (matrix[i][j] != 0)
            {
                ai[index] = i+1;
                aj[index] = j+1;
                ar[index] = DOUBLE(matrix[i][j]);
                ++index;
            }
        }
    }
    lpx_load_matrix(lp, index-1, ai, aj, ar);
    delete [] ai; delete [] aj; delete [] ar;
}

void
_4ti2_::load_matrix_transpose(LPX* lp, const VectorArray& matrix)
{
    int m = matrix.get_size();
    int n = matrix.get_number();
    // Input the matrix.
    int *ai = new int[m*n+1];
    int *aj = new int[m*n+1];
    double *ar = new double[m*n+1];
    int index = 1;
    for (int i = 1; i <= m; ++i)
    {
        for (int j = 1; j <= n; ++j)
        {
            if (matrix[j-1][i-1] != 0)
            {
                ai[index] = i;
                aj[index] = j;
                ar[index] = DOUBLE(matrix[j-1][i-1]);
                //*out << "(" << r << "," << c << "," << lattice[j-1][i-1] << ")\n";
                ++index;
            }
        }
    }
    lpx_load_matrix(lp, index-1, ai, aj, ar);
    delete [] ai; delete [] aj; delete [] ar;
}
