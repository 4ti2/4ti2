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

#ifndef _4ti2_zsolve__Variables_
#define _4ti2_zsolve__Variables_

#include <iostream>
#include <vector>
#include <cassert>

#include "zsolve/Integer.h"

namespace _4ti2_zsolve_
{

template <typename T> class VariableProperty
{
protected:
    int m_column_id;
    bool m_is_free;
    T m_upper_bound;
    T m_lower_bound;

public:
    VariableProperty (const int column, const bool free, const T& lower, const T& upper)
    {
        set (column, free, lower, upper);
    }

    VariableProperty (const VariableProperty& other)
    {
        set (other.m_column_id, other.m_is_free, other.m_lower_bound, other.m_upper_bound);
    }

    bool is_symmetric () const
    {
        return m_upper_bound + m_lower_bound == 0;
    }

    void set_bound (bool set_lower, const T& value)
    {
        if (m_is_free)
            return;

        if (set_lower)
            m_lower_bound = value;
        else
            m_upper_bound = value;
    }

    void set (const VariableProperty& other)
    {
        set (other.m_column_id, other.m_is_free, other.m_lower_bound, other.m_upper_bound);
    }

    void set (const int column, const bool free, const T& lower, const T& upper)
    {
        m_column_id = column;
        m_is_free = free;
        m_lower_bound = lower;
        m_upper_bound = upper;
        
    }
    
    void set (const bool free, const T& lower = 1, const T& upper = -1)
    {
        m_is_free = free;
        m_lower_bound = lower;
        m_upper_bound = upper;
        
    }
    
    bool free () const
    {
        return m_is_free;
    }

    int column () const
    {
        return m_column_id;
    }

    T lower () const
    {
        return m_lower_bound;
    }

    T upper () const
    {
        return m_upper_bound;
    }

    int compare (const VariableProperty& other) const
    {
        int max = (m_column_id > other.m_column_id ? m_column_id : other.m_column_id) + 1;
        int c1 = m_column_id < 0 ? max - m_column_id : m_column_id;
        int c2 = other.m_column_id < 0 ? max - other.m_column_id : other.m_column_id;

        return c1 - c2;
    }

    int count_infinity () const
    {
        int result = 2;
        if (m_lower_bound <= 0)
            result--;
        if (m_upper_bound >= 0)
            result--;
        return result;
    }

    T get_range ()
    {
        T result = 0;
        if (m_upper_bound > 0)
            result += m_upper_bound;
        if (m_lower_bound < 0)
            result -= m_lower_bound;
        return result;
    }

    bool check_bounds (const T& value)
    {
        if (m_lower_bound <= 0 && value < m_lower_bound)
            return false;
        if (m_upper_bound >= 0 && value > m_upper_bound)
            return false;
        return true;
    }

    bool check_consistency () const
    {
        return m_column_id >= -2;
    }

    int lower_space () const
    {
        return m_lower_bound < 0 ? integer_space (m_lower_bound) : 1;
    }
    
    int upper_space () const
    {
        return m_upper_bound > 0 ? integer_space (m_upper_bound) : 1;
    }

    std::ostream& upper (std::ostream& out) const
    {
        if (m_upper_bound >= 0)
            out << m_upper_bound;
        else
            out << "+";
        return out;
    }
    
    std::ostream& lower (std::ostream& out) const
    {
        if (m_lower_bound <= 0)
            out << m_lower_bound;
        else
            out << "-";
        return out;
    }

    std::ostream& dump (std::ostream& out) const
    {
        out << m_column_id;
        out << (m_is_free ? " 1 " : " 0 ");
        out << m_lower_bound;
        out << " ";
        out << m_upper_bound;
        return out;
    }
};

template <typename T> class VariableProperties
{
protected:
    std::vector <VariableProperty <T> *> m_variable_properties;

public:
    VariableProperties (size_t variables, bool free, const T& lower, const T& upper)
    {
        m_variable_properties.resize (variables);
        for (size_t i = 0; i < variables; i++)
        {
            m_variable_properties[i] = new VariableProperty <T> (i, free, lower, upper);
        }
    }

    VariableProperties (VariableProperties <T>* other)
    {
        m_variable_properties.resize (other->variables ());
        for (size_t i = 0; i < other->variables (); i++)
        {
            m_variable_properties[i] = new VariableProperty <T> (other->get_variable (i));
        }
    }

    ~VariableProperties ()
    {
        for (size_t i = 0; i < m_variable_properties.size (); i++)
            delete m_variable_properties[i];
        m_variable_properties.clear ();
    }

    VariableProperty <T>& get_variable (size_t index)
    {
        return *m_variable_properties[index];
    }

    size_t variables () const
    {
        return m_variable_properties.size ();
    }

    void swap (size_t a, size_t b)
    {
        VariableProperty <T> * temp = m_variable_properties[a];
        m_variable_properties[a] = m_variable_properties[b];
        m_variable_properties[b] = temp;
    }

    bool check_consistency () const
    {
        for (size_t i = 0; i < m_variable_properties.size (); i++)
            if (!m_variable_properties[i].check_consistency ())
                return false;
        return true;
    }
};

} // namespace _4ti2_zsolve_

#endif
