/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Matthias Walter.

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

#ifndef _4ti2__LinearSystem_
#define _4ti2__LinearSystem_

#include <iostream>
#include <vector>
#include <cassert>

#include "Integer.h"
#include "VectorArray.hpp"
#include "Variables.hpp"
#include "Relation.hpp"
#include "Vector.hpp"

template <typename T> class LinearSystem : public VariableProperties <T>
{
protected:
    std::vector <Relation <T> *> m_relation_properties;
    size_t m_relations;
    VectorArray <T>* m_matrix;
    T* m_rhs;

public:
    LinearSystem (const VectorArray <T>& matrix, T* rhs, bool free, const T& lower, const T& upper) : VariableProperties <T> (matrix.width (), free, lower, upper)
    {
        m_matrix = new VectorArray <T> (matrix);
        m_rhs = copy_vector (rhs, matrix.height ());
        m_relations = m_matrix->height ();
        m_relation_properties.resize (m_relations);
        for (size_t i = 0; i < m_relations; i++)
        {
            m_relation_properties[i] = new Relation <T> (Relation <T> ::Equal);
        }

        assert (check_consistency ());
    }

    ~LinearSystem ()
    {
        delete m_matrix;
        delete_vector (m_rhs);
        for (size_t i = 0; i < m_relations; i++)
            delete m_relation_properties[i];
        m_relation_properties.clear ();
    }

    size_t relations () const
    {
        return m_relations;
    }

    size_t variables () const
    {
        return VariableProperties <T> :: m_variable_properties.size ();
    }

    Relation<T> & get_relation (const size_t index)
    {
        return *m_relation_properties[index];
    }

    VectorArray <T>& matrix () const
    {
        return *m_matrix;
    }

    T* rhs () const
    {
        return m_rhs;
    }

    bool is_homogeneous () const
    {
        return is_zero_vector (m_rhs, m_relations);
    }

    bool is_equality_system () const
    {
        for (size_t i = 0; i < m_relations; i++)
            if (! m_relation_properties[i]->is_equality ())
                return false;
        return true;
    }
    
    bool check_consistency () const
    {
        if (!m_matrix->check_consistency ())
            return false;
        
        if (!check_vector_consistency (m_rhs, m_relations))
            return false;

        if (m_matrix->height () != m_relations)
            return false;

        if (m_matrix->width () != VariableProperties <T> :: m_variable_properties.size ())
            return false;

        if (m_relations != m_relation_properties.size ())
            return false;

        for (size_t i = 0; i < m_relation_properties.size (); i++)
            if (!m_relation_properties[i]->check_consistency ())
                return false;

        return true;
    }
    
    template <typename X> friend std::ostream& operator<< (std::ostream& out, LinearSystem <X>& system);
};

template <typename T> std::ostream& operator<< (std::ostream& out, LinearSystem <T>& system)
{
    size_t vars = system.variables ();
    size_t rels = system.relations ();
    size_t* space = new size_t[vars + 2];

    for (size_t i = 0; i < vars; i++)
    {
        VariableProperty <T> & var = system.get_variable (i);
        space[i] = max (var.upper_space (), var.lower_space ());
        for (size_t j = 0; j < rels; j++)
        {
            space[i] = max <size_t> (space[i], integer_space(system.matrix () [j][i]));
        }
    }
    space[vars] = 1;
    space[vars+1] = 1;
    for (size_t i = 0; i < rels; i++)
    {
        Relation <T> & rel = system.get_relation (i);

        space[vars] = max <size_t> (space[vars], rel.space ());
        space[vars+1] = max <size_t> (space[vars+1], integer_space (system.rhs () [i]));
    }

    // print variables lines
    for (size_t i = 0; i < vars; i++)
    {
        VariableProperty <T> & var = system.get_variable (i);
        if (i > 0)
            out << " ";
        for (int j = space[i] - var.upper_space (); j > 0; j--)
            out << " ";
        var.upper (out);
    }
    out << "\n";
    for (size_t i = 0; i < vars; i++)
    {
        VariableProperty <T> & var = system.get_variable (i);
        if (i > 0)
            out << " ";
        for (int j = space[i] - var.lower_space (); j > 0; j--)
            out << " ";
        var.lower (out);
    }
    out << "\n";
    for (size_t i = 0; i < vars; i++)
    {
        VariableProperty <T> & var = system.get_variable (i);
        if (i > 0)
            out << " ";
        for (int j = space[i] - 1; j > 0; j--)
            out << " ";
        if (var.free ())
            out << "F";
        else if (var.lower () > 0 && var.upper () < 0)
            out << "G";
        else if (var.upper () < 0)
            out << "H";
        else if (var.lower () == 0 && var.upper () == 1)
            out << "B";
        else
            out << " ";
    }
    out << "\n";
    for (size_t i = 0; i < rels; i++)
    {
        out << "\n";
        for (size_t j = 0; j < vars; j++)
        {
            if (j > 0)
                out << " ";
            const T value = system.matrix () [i][j];
            for (int k = space[j] - integer_space (value); k > 0; k--)
                out << " ";
            out << value;
        }
        out << " ";
        Relation <T> & rel = system.get_relation (i);
        for (int k = space[vars] - rel.space (); k > 0; k--)
            out << " ";
        rel.print (out);
        out << " ";
        const T value = system.rhs () [i];
        for (int k = space[vars+1] - integer_space (value); k > 0; k--)
            out << " ";
        out << value;
    }

    out << "\n" << std::flush;

    delete[] space;

    return out;
}

template <typename T> LinearSystem <T>* homogenize_linear_system (LinearSystem <T>* other)
{
    int inequalities = 0;
    bool inhom = false;
    T* rhs = copy_vector <T> (other->rhs (), other->relations ());
    
    for (size_t i = 0; i < other->relations (); i++)
    {
        const Relation <T> rel = other->get_relation (i);
        rhs[i] += rel.get_adjustment ();
        if (! rel.is_equality ())
            inequalities++;
        if (rhs[i] != 0)
            inhom = true;
    }

    VectorArray <T> matrix (other->relations (), other->variables () + inequalities + (inhom ? 1 : 0));

    // fill old matrix

    for (size_t col = 0; col < other->matrix ().width (); col++)
    {
        for (size_t row = 0; row < other->matrix ().height (); row++)
        {
            matrix[row][col] = other->matrix() [row][col];
        }
    }

    // fill slack area

    size_t current = other->variables ();
    for (size_t i = 0; i < other->relations (); i++)
    {
        const Relation <T>& rel = other->get_relation (i);
        if (rel.is_equality ())
            continue;
        for (size_t j = 0; j < other->relations (); j++)
        {
            matrix[j][current] = (i == j) ? rel.get_slack_value () : 0;
        }
        current++;
    }

    // set rhs

    if (inhom)
    {
        for (size_t i = 0; i < other->relations (); i++)
        {
            matrix[i][current] = -rhs [i];
            rhs[i] = 0;
        }
    }

    // create linear system

    LinearSystem <T> * system = new LinearSystem <T> (matrix, rhs, true, 1, -1);

    for (size_t i = 0; i < other->variables (); i++)
    {
        system->get_variable (i).set (other->get_variable (i));
    }

    // slack columns
    current = other->variables ();
    for (size_t i = 0; i < other->relations (); i++)
    {
        const Relation <T>& rel = other->get_relation (i);
        if (rel.is_equality ())
            continue;
        VariableProperty <T>& var = system->get_variable (current);
        var.set (-1, false, rel.get_slack_bounds (true), rel.get_slack_bounds (false));
        current++;
    }
    // rhs-column
    if (inhom)
    {
        VariableProperty <T>& var = system->get_variable (current);
        var.set (-2, false, 0, 1);
    }

    delete_vector (rhs);

    return system;
}

#endif
