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

#include "Optimise.h"
#include "Statistics.h"
#include "Bounded.h"
#include "Globals.h"

#include "BasicCompletion.h"
#include "SyzygyCompletion.h"
#include "HermiteAlgorithm.h"
#include "Completion.h"
#include "Minimize.h"
#include "Timer.h"
#include "Markov.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include "BitSetStream.h"
#include "VectorArrayStream.h"
#include "FeasibleStream.h"

//#define DEBUG_4ti2(X) X
#include "Debug.h"

using namespace _4ti2_;

Optimise::Optimise()
{
}

Optimise::~Optimise()
{
}

int
Optimise::compute(
                Feasible& feasible,
                Vector& cost,
                Vector& sol)
{
    BitSet rs(feasible.get_urs());
    rs.set_complement();
    // We use different algorithms depending on whether the given solution is
    // feasible or not.
    if (sol.is_non_negative(rs))
    {
        return compute_feasible(feasible, cost, sol);
    }
    else
    {
        return compute_infeasible(feasible, cost, sol);
    }
}

// TODO: Deal with unbounded components even after cost introduced.
// We can project out these components.
int
Optimise::compute_infeasible(
                Feasible& feasible,
                Vector& cost,
                Vector& sol)
{
    assert(feasible.get_dimension() == cost.get_size());
    assert(feasible.get_dimension() == sol.get_size());

    Timer t;

    *out << "Optimizing.\n";

    int n = feasible.get_dimension();
    BitSet rs(feasible.get_urs());
    rs.set_complement();

    Vector rhs(feasible.get_matrix().get_number());
    VectorArray::dot(feasible.get_matrix(), sol, rhs);
    BitSet basic(n);
    RationalType objective;
    int status = lp_solve(feasible.get_matrix(), rhs, cost, feasible.get_urs(),
                    basic, objective);
    switch (status)
    {
        case 0:
            *out << "Objective value = " << objective << "\n";
            DEBUG_4ti2(*out << "Basic variables:\n" << basic << "\n";)
            break;
        case -1:
            *out << "Problem is infeasible.\n";
            return -1;
        case 1:
            *out << "Problem is unbounded.\n";
            return 1;
        default:
            std::cerr << "Software Error: Unexpected LP solver output.\n";
            exit(1);
    }

    VectorArray costs(0, n);
    costs.insert(cost);
    VectorArray sols(0, n);
    sols.insert(sol);

    // We project onto the non-basic variables.
    BitSet proj(feasible.get_urs());
    proj.set_union(basic);
    Feasible* projection = new Feasible(feasible, proj);
    Feasible* next_projection = 0;

    DEBUG_4ti2(*out << *projection << "\n";)
    VectorArray gens(feasible.get_basis());

    assert(projection->get_bnd().empty());
    DEBUG_4ti2(*out << "Bounded:\n" << projection->get_bnd() << "\n";)
    DEBUG_4ti2(*out << "Unbounded:\n" << projection->get_unbnd() << "\n";)

    // TODO: We don't need to do the following step. We can just compute an
    // upper triangle form of the basis.
    gens.insert(projection->get_ray());
    make_feasible(sols, projection->get_ray()); 
    DEBUG_4ti2(*out << "Solution:\n" << sols[0] << "\n";)

    *out << "Computing Groebner basis for the group relaxation...\n";
    Completion algorithm;
    algorithm.compute(*projection, costs, gens, sols);
    *out << "Optimal Solution:\n" << sols[0] << "\n";
    *out << "Objective = " << Vector::dot(sols[0], cost) << "\n";
    if (sols[0].is_non_negative(rs))
    {
        *out << "Solution is optimal.\n";
        sol = sols[0];
        DEBUG_4ti2(*out << sols[0] << "\n";)
        *out << "Done. " << " Time: " << t << " / ";
        *out << Timer::global << " secs.\n";
        delete projection;
        return 0;
    }

    *out << "Computing Groebner bases for the extended group relaxations...\n";
    BitSet remaining(basic);
    while (!remaining.empty())
    {
        DEBUG_4ti2(*out << "NEXT ITERATION\n";)
        DEBUG_4ti2(*out << "Left = " << remaining.count() << "\n";)
        DEBUG_4ti2(*out << "Remaining:\n" << remaining << "\n";)
        DEBUG_4ti2(*out << "Projection:\n" << proj << "\n";)

        int column = next_support(gens, remaining, sols[0]);
        DEBUG_4ti2(*out << "Next column is " << column << "\n";)
        proj.unset(column);
        remaining.unset(column);

        next_projection = new Feasible(feasible, proj);
        if (next_projection->get_bnd()[column])
        {
            DEBUG_4ti2(*out << "Column is bounded.\n";)
            DEBUG_4ti2(*out << "Unbounded:\n";)
            DEBUG_4ti2(*out << projection->get_unbnd() << "\n";)
            DEBUG_4ti2(*out << "Bounded:\n";)
            DEBUG_4ti2(*out << projection->get_bnd() << "\n";)
            VectorArray tmp_cost(1,n,0);
            tmp_cost[0][column] = -1;
            tmp_cost.insert(cost);
            Completion algorithm;
            algorithm.compute(*projection, tmp_cost, gens, sols);
            if (sols[0][column] < 0)
            {
                *out << "Problem is infeasible.\n";
                *out << "Done. " << " Time: " << t << " / ";
                *out << Timer::global << " secs.\n";
                delete projection;
                delete next_projection;
                return -1;
            }
        }
        else
        {
            DEBUG_4ti2(*out << "Column is unbounded.\n";)
            DEBUG_4ti2(*out << "Ray:\n";)
            DEBUG_4ti2(*out << next_projection->get_ray() << "\n";)
            // TODO: Check for duplicates.
            // We should check for a ray first.
            gens.insert(next_projection->get_ray());
            make_feasible(sols, next_projection->get_ray());
        }

        Completion algorithm;
        algorithm.compute(*next_projection, costs, gens, sols);
        DEBUG_4ti2(*out << "Optimal Solution:\n" << sols[0] << "\n";)
        *out << "Objective = " << Vector::dot(sols[0], cost) << "\n";
        if (sols[0].is_non_negative(rs))
        {
            *out << "Solution is optimal.\n";
            DEBUG_4ti2(*out << sols[0] << "\n";)
            sol = sols[0];
            *out << "Objective = " << Vector::dot(sols[0], cost) << "\n";
            *out << "Done. " << " Time: " << t << " / ";
            *out << Timer::global << " secs.\n";
            delete projection;
            delete next_projection;
            return 0;
        }

        delete projection;
        projection = next_projection;
        next_projection = 0;
    }
    delete projection;
    std::cerr << "Software Error: Unexpected program execution.\n";
    exit(1);
    return 0;
}

// TODO: Deal with unbounded components even after cost introduced.
// We can project out these components.
int
Optimise::compute_feasible(
                Feasible& feasible,
                Vector& cost,
                Vector& sol)
{
    DEBUG_4ti2(*out << feasible << "\n";)
    DEBUG_4ti2(*out << "Cost Vector:\n" << cost << "\n";)
    DEBUG_4ti2(*out << "Solution:\n" << sol << "\n";)

    // Add the cost function to the matrix.
    const VectorArray& matrix = feasible.get_matrix();
    VectorArray ext_matrix(matrix.get_number(),matrix.get_size()+1,0);
    VectorArray::lift(matrix, 0, matrix.get_size(), ext_matrix);
    Vector cost_row(cost.get_size()+1);
    Vector::lift(cost, 0, cost.get_size(), cost_row);
    cost_row[cost.get_size()] = 1;
    ext_matrix.insert(cost_row);
    DEBUG_4ti2(*out << "Full Matrix:\n" << ext_matrix << "\n";)

    // Add the cost function to the lattice.
    const VectorArray& lattice = feasible.get_basis();
    VectorArray ext_lattice(lattice.get_number(),lattice.get_size()+1);
    VectorArray::lift(lattice, 0, lattice.get_size(), ext_lattice);
    Vector cost_col(lattice.get_number());
    VectorArray::dot(lattice, cost, cost_col);
    for (int i = 0; i < ext_lattice.get_number(); ++i)
    {   ext_lattice[i][lattice.get_size()] = -cost_col[i]; }
    DEBUG_4ti2(*out << "Full Lattice:\n" << ext_lattice << "\n";)

    // Extend the urs.
    const BitSet& urs = feasible.get_urs();
    BitSet ext_urs(urs.get_size()+1);
    BitSet::extend(urs, ext_urs);
    DEBUG_4ti2(*out << "Extended URS:\n" << ext_urs << "\n";)

    // Extend the solution.
    Vector ext_sol(sol.get_size()+1,0);
    Vector::lift(sol, 0, sol.get_size(), ext_sol);
    DEBUG_4ti2(*out << "Extended Solution:\n" << ext_sol << "\n";)

    // TODO: We do not transfer weights!!
    Feasible ext_feasible(&ext_lattice, &ext_matrix, &ext_urs, &ext_sol);

    IntegerType cost_offset = Vector::dot(cost, sol);

    int status = compute_feasible(ext_feasible, sol.get_size(), cost_offset, ext_sol);
    Vector::project(ext_sol, 0, sol.get_size(), sol);
    return status;
}

// TODO: Deal with unbounded components even after cost introduced.
// We can project out these components.
int
Optimise::compute_feasible(
                Feasible& feasible,
                int cost_col, IntegerType cost_offset,
                Vector& sol)
{
    assert(0 <= cost_col && cost_col < feasible.get_dimension());
    assert(feasible.get_dimension() == sol.get_size());
    DEBUG_4ti2(*out << feasible << "\n";)
    DEBUG_4ti2(*out << "Cost column is " << cost_col << "\n";)
    DEBUG_4ti2(*out << "Solution:\n" << sol << "\n";)
    *out << "Upper Bound = " << cost_offset - sol[cost_col] << "\n";

    Timer t;

    *out << "Optimizing.\n";

    int n = feasible.get_dimension();
    BitSet rs(feasible.get_urs());
    rs.set_complement();

    Vector cost(n,0);
    cost[cost_col] = -1;

    Vector rhs(feasible.get_matrix().get_number());
    VectorArray::dot(feasible.get_matrix(), sol, rhs);
    BitSet basic(n);
    RationalType objective;
    int status = lp_solve(feasible.get_matrix(), rhs, cost, feasible.get_urs(),
                    basic, objective);
    switch (status)
    {
        case 0:
            *out << "LP Objective value/Lower Bound = " << objective << "\n";
            DEBUG_4ti2(*out << "Basic variables:\n" << basic << "\n";)
            break;
        case -1:
            *out << "Problem is infeasible.\n";
            return -1;
        case 1:
            *out << "Problem is unbounded.\n";
            return 1;
        default:
            std::cerr << "Software Error: Unexpected LP solver output.\n";
            exit(1);
    }

    VectorArray costs(0, n);
    costs.insert(cost);

    // We project onto the non-basic variables.
    BitSet proj(feasible.get_urs());
    proj.set_union(basic);
    proj.set(cost_col);
    Feasible* projection = new Feasible(feasible, proj);
    DEBUG_4ti2(*out << *projection << "\n";)

    VectorArray gens(feasible.get_basis());

    assert(projection->get_bnd().empty());
    DEBUG_4ti2(*out << "Bounded:\n" << projection->get_bnd() << "\n";)
    DEBUG_4ti2(*out << "Unbounded:\n" << projection->get_unbnd() << "\n";)

    // We add a ray so that gens is a generating set.
    // TODO: We don't need to do the following step. We can just compute an
    // upper triangle form of the basis.
    gens.insert(projection->get_ray());

    *out << "Solving the group relaxation...\n";
    VectorArray sols(0, n);
    sols.insert(sol);
    Completion algorithm;
    algorithm.compute(*projection, costs, gens, sols);
    *out << "Optimal Solution:\n" << sols[0] << "\n";
    *out << "Objective = " << cost_offset - sols[0][cost_col] << "\n";
    if (sols[0].is_non_negative(rs))
    {
        *out << "Solution is optimal.\n";
        sol = sols[0];
        DEBUG_4ti2(*out << sols[0] << "\n";)
        *out << "Done. " << " Time: " << t << " / ";
        *out << Timer::global << " secs.\n";
        delete projection;
        return 0;
    }

#if 0
    VectorArray tmp_sols(0,n);
    tmp_sols.insert(sol);
    Minimize minimize;
    minimize.minimize(feasible, costs, gens, tmp_sols);
    *out << "Partial Optimal:\n" << tmp_sols[0] << "\n";
    *out << "Objective = " << cost_offset - tmp_sols[0][cost_col] << "\n";
#endif

    *out << "Solving the extended group relaxations...\n";
    proj.unset(cost_col);
    delete projection;
    projection = new Feasible(feasible, proj);
    Feasible* next_projection = 0;
    DEBUG_4ti2(*out << *projection << "\n";)

    BitSet remaining(basic);
    while (!remaining.empty())
    {
        DEBUG_4ti2(*out << "NEXT ITERATION\n";)
        DEBUG_4ti2(*out << "Left = " << remaining.count() << "\n";)
        DEBUG_4ti2(*out << "Remaining:\n" << remaining << "\n";)
        DEBUG_4ti2(*out << "Projection:\n" << proj << "\n";)

        int column = next_support(gens, remaining, sols[0]);
        DEBUG_4ti2(*out << "Next column is " << column << "\n";)
        proj.unset(column);
        remaining.unset(column);

        *out << "Computing Generating Set...\n";
        next_projection = new Feasible(feasible, proj);
        if (next_projection->get_bnd()[column])
        {
            DEBUG_4ti2(*out << "Column is bounded.\n";)
            DEBUG_4ti2(*out << "Unbounded:\n";)
            DEBUG_4ti2(*out << projection->get_unbnd() << "\n";)
            DEBUG_4ti2(*out << "Bounded:\n";)
            DEBUG_4ti2(*out << projection->get_bnd() << "\n";)
            VectorArray tmp_cost(1,n,0);
            tmp_cost[0][column] = -1;
            tmp_cost.insert(cost);
            Completion algorithm;
            algorithm.compute(*projection, tmp_cost, gens);
        }
        else
        {
            DEBUG_4ti2(*out << "Column is unbounded.\n";)
            DEBUG_4ti2(*out << "Ray:\n";)
            DEBUG_4ti2(*out << next_projection->get_ray() << "\n";)
            // TODO: Check for duplicates.
            // We should check for a ray first.
            // We add a ray so that gens is a generating set.
            gens.insert(next_projection->get_ray());
        }

        *out << "Computing Groebner basis...\n";
        DEBUG_4ti2(*out << *next_projection << "\n";)
        VectorArray sols(0, n);
        sols.insert(sol);
        Completion algorithm;
        algorithm.compute(*next_projection, costs, gens, sols);
        DEBUG_4ti2(*out << "Optimal Solution:\n" << sols[0] << "\n";)
        *out << "Objective = " << cost_offset - sols[0][cost_col] << "\n";
        if (sols[0].is_non_negative(rs))
        {
            *out << "Solution is optimal.\n";
            DEBUG_4ti2(*out << sols[0] << "\n";)
            sol = sols[0];
            *out << "Done. " << " Time: " << t << " / ";
            *out << Timer::global << " secs.\n";
            delete projection;
            delete next_projection;
            return 0;
        }

        delete projection;
        projection = next_projection;
        next_projection = 0;
    }
    delete projection;
    std::cerr << "Software Error: Unexpected program execution.\n";
    exit(1);
    return 0;
}

// TODO: Deal with unbounded components even after cost introduced.
// We can project out these components.
int
Optimise::compute_bounded(
                Feasible& feasible,
                int cost_col, IntegerType cost_offset,
                Vector& sol, IntegerType upper_bound)
{
    assert(0 <= cost_col && cost_col < feasible.get_dimension());
    assert(feasible.get_dimension() == sol.get_size());
    DEBUG_4ti2(*out << feasible << "\n";)
    DEBUG_4ti2(*out << "Cost column is " << cost_col << "\n";)
    DEBUG_4ti2(*out << "Solution:\n" << sol << "\n";)

    Timer t;

    *out << "Optimizing.\n";

    int n = feasible.get_dimension();
    BitSet rs(feasible.get_urs());
    rs.set_complement();

    Vector cost(n,0);
    cost[cost_col] = -1;

    Vector rhs(feasible.get_matrix().get_number());
    VectorArray::dot(feasible.get_matrix(), sol, rhs);
    BitSet basic(n);
    RationalType objective;
    int status = lp_solve(feasible.get_matrix(), rhs, cost, feasible.get_urs(),
                    basic, objective);
    switch (status)
    {
        case 0:
            *out << "LP Objective value = " << objective << "\n";
            DEBUG_4ti2(*out << "Basic variables:\n" << basic << "\n";)
            break;
        case -1:
            *out << "Problem is infeasible.\n";
            return -1;
        case 1:
            *out << "Problem is unbounded.\n";
            return 1;
        default:
            std::cerr << "Software Error: Unexpected LP solver output.\n";
            exit(1);
    }

    VectorArray costs(0, n);
    costs.insert(cost);
    VectorArray sols(0, n);
    sols.insert(sol);

    // We project onto the non-basic variables.
    BitSet proj(feasible.get_urs());
    proj.set_union(basic);
    proj.set(cost_col);
    Feasible* projection = new Feasible(feasible, proj);
    DEBUG_4ti2(*out << *projection << "\n";)

    VectorArray gens(feasible.get_basis());

    assert(projection->get_bnd().empty());
    DEBUG_4ti2(*out << "Bounded:\n" << projection->get_bnd() << "\n";)
    DEBUG_4ti2(*out << "Unbounded:\n" << projection->get_unbnd() << "\n";)

    // TODO: We don't need to do the following step. We can just compute an
    // upper triangle form of the basis.
    gens.insert(projection->get_ray());
    make_feasible(sols, projection->get_ray()); 
    DEBUG_4ti2(*out << "Solution:\n" << sols[0] << "\n";)

    *out << "Solving the group relaxation...\n";
    Completion algorithm;
    algorithm.compute(*projection, costs, gens, sols);
    *out << "Optimal Solution:\n" << sols[0] << "\n";
    *out << "Objective = " << cost_offset - sols[0][cost_col] << "\n";
    if (sols[0].is_non_negative(rs))
    {
        *out << "Solution is optimal.\n";
        sol = sols[0];
        DEBUG_4ti2(*out << sols[0] << "\n";)
        *out << "Done. " << " Time: " << t << " / ";
        *out << Timer::global << " secs.\n";
        delete projection;
        return 0;
    }

#if 0
    VectorArray tmp_sols(0,n);
    tmp_sols.insert(sol);
    Minimize minimize;
    minimize.minimize(feasible, costs, gens, tmp_sols);
    *out << "Partial Optimal:\n" << tmp_sols[0] << "\n";
    *out << "Objective = " << cost_offset - tmp_sols[0][cost_col] << "\n";
#endif

    *out << "Solving the extended group relaxations...\n";
    proj.unset(cost_col);
    delete projection;
    projection = new Feasible(feasible, proj);
    Feasible* next_projection = 0;
    DEBUG_4ti2(*out << *projection << "\n";)

    BitSet remaining(basic);
    while (!remaining.empty())
    {
        DEBUG_4ti2(*out << "NEXT ITERATION\n";)
        DEBUG_4ti2(*out << "Left = " << remaining.count() << "\n";)
        DEBUG_4ti2(*out << "Remaining:\n" << remaining << "\n";)
        DEBUG_4ti2(*out << "Projection:\n" << proj << "\n";)

        int column = next_support(gens, remaining, sols[0]);
        DEBUG_4ti2(*out << "Next column is " << column << "\n";)
        proj.unset(column);
        remaining.unset(column);

        *out << "Computing Generating Set...\n";
        next_projection = new Feasible(feasible, proj);
        if (next_projection->get_bnd()[column])
        {
            DEBUG_4ti2(*out << "Column is bounded.\n";)
            DEBUG_4ti2(*out << "Unbounded:\n";)
            DEBUG_4ti2(*out << projection->get_unbnd() << "\n";)
            DEBUG_4ti2(*out << "Bounded:\n";)
            DEBUG_4ti2(*out << projection->get_bnd() << "\n";)
            VectorArray tmp_cost(1,n,0);
            tmp_cost[0][column] = -1;
            tmp_cost.insert(cost);
            Completion algorithm;
            algorithm.compute(*projection, tmp_cost, gens, sols);
            if (sols[0][column] < 0)
            {
                *out << "Problem is infeasible.\n";
                *out << "Done. " << " Time: " << t << " / ";
                *out << Timer::global << " secs.\n";
                delete projection;
                delete next_projection;
                return -1;
            }
        }
        else
        {
            DEBUG_4ti2(*out << "Column is unbounded.\n";)
            DEBUG_4ti2(*out << "Ray:\n";)
            DEBUG_4ti2(*out << next_projection->get_ray() << "\n";)
            // TODO: Check for duplicates.
            // We should check for a ray first.
            gens.insert(next_projection->get_ray());
            make_feasible(sols, next_projection->get_ray());
        }

        *out << "Computing Groebner basis...\n";
        Completion algorithm;
        algorithm.compute(*next_projection, costs, gens, sols);
        DEBUG_4ti2(*out << "Optimal Solution:\n" << sols[0] << "\n";)
        *out << "Objective = " << cost_offset - sols[0][cost_col] << "\n";
        if (sols[0].is_non_negative(rs))
        {
            *out << "Solution is optimal.\n";
            DEBUG_4ti2(*out << sols[0] << "\n";)
            sol = sols[0];
            *out << "Done. " << " Time: " << t << " / ";
            *out << Timer::global << " secs.\n";
            delete projection;
            delete next_projection;
            return 0;
        }

        delete projection;
        projection = next_projection;
        next_projection = 0;
    }
    delete projection;
    std::cerr << "Software Error: Unexpected program execution.\n";
    exit(1);
    return 0;
}

int
Optimise::add_support(
                const VectorArray& gens,
                BitSet& fin)
{
    int num_lifts = 0;
    for (Index c = 0; c < gens.get_size(); ++c)
    {
        if (fin[c] == 1 && positive_count(gens, c) == 0)
        {
            fin.unset(c);
            ++num_lifts;
        }
    }
    if (num_lifts != 0) 
    {
        *out << "  Lifted already on " << num_lifts << " variable(s).";
        *out << std::endl;
    }
    return num_lifts;
}

int
Optimise::positive_count(
                const VectorArray& gens,
                int c)
{
    int count = 0;
    for (Index r = 0; r < gens.get_number(); ++r)
    {
        if (gens[r][c] > 0) ++count;
    }
    return count;
}

int
Optimise::next_support(
                const VectorArray& gens,
                const BitSet& fin,
                const Vector& sol)
{
#if 1
    int col = -1;
    IntegerType min = 0;
    for (Index c = 0; c < gens.get_size(); ++c)
    {
        if (fin[c] && sol[c] < min)
        {
            min = sol[c];
            col = c;
        }
    }
    assert(col != -1);
    return col;
#endif
#if 0
    int min_count = gens.get_number()+1;
    int col = -1;
    for (Index c = 0; c < gens.get_size(); ++c)
    {
        if (fin[c])
        {
            int count = positive_count(gens, c);
            if (count < min_count)
            {
                col = c;
                min_count = count;
            }
        }
    }
    assert(col != -1);
    assert(min_count != (int) gens.get_number()+1);
    return col;
#endif
}

void
Optimise::make_feasible(
                VectorArray& feasibles,
                const Vector& ray)
{
    DEBUG_4ti2(*out << "Making feasible.\n";)
    DEBUG_4ti2(*out << "Feasibles:\n" << feasibles << "\n";)
    DEBUG_4ti2(*out << "Ray:\n" << ray << "\n";)
    IntegerType factor = 0;
    for (int i = 0; i < feasibles.get_number(); ++i)
    {
        for (int j = 0; j < ray.get_size(); ++j)
        {
            if (feasibles[i][j] < 0 && ray[j] > 0)
            {
                IntegerType ratio = -feasibles[i][j]/ray[j] + 1;
                if (ratio > factor) { factor = ratio; }
            }
        }
        if (factor != 0) { feasibles[i].add(ray, factor); } 
    }
}
