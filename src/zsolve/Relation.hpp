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

#ifndef _4ti2_zsolve__Relation_
#define _4ti2_zsolve__Relation_

#include <cassert>
#include <iostream>

#include "zsolve/Integer.h"

namespace _4ti2_zsolve_
{

template <typename T> class Relation
{
public:
    enum RelationType { Equal, Lesser, LesserEqual, Greater, GreaterEqual, Modulo};

protected:
    RelationType m_relation_type;
    T m_modulus;

public:

    Relation (RelationType type)
    {
        set (type);
    }

    Relation (const T& modulus)
    {
        set (modulus);
    }

    void set (RelationType type)
    {
        assert (type != Modulo);

        m_relation_type = type;
        m_modulus = 0;
    }

    void set (const T& modulus)
    {
        m_relation_type = Modulo;
        m_modulus = modulus;
    }

    bool is_equality () const
    {
        return m_relation_type == Equal;
    }

    T get_slack_value () const
    {
        switch (m_relation_type)
        {
        case Equal:
            return 0;
        case Modulo:
            return m_modulus;
        case Lesser:
        case LesserEqual:
            return 1;
        case Greater:
        case GreaterEqual:
            return -1;
        default:
            assert (false);
            return 0;
        }
    }

    int get_slack_bounds (bool is_lower) const
    {
        if (m_relation_type == Modulo)
            return is_lower ? 1 : -1;
        else if (m_relation_type == Equal)
            return 0;
        else
            return is_lower ? 0 : -1;
    }

    int get_adjustment () const
    {
        if (m_relation_type == Lesser)
            return -1;
        else if (m_relation_type == Greater)
            return 1;
        else
            return 0;
    }

    bool check_consistency () const
    {
        return true;
    }

    int space () const
    {
        if (m_relation_type == LesserEqual || m_relation_type == GreaterEqual)
            return 2;
        else
            return 1;
    }

    std::ostream& print (std::ostream& out) const
    {
        switch (m_relation_type)
        {
        case Equal:
            out << "=";
        break;
        case Modulo:
            out << "=";
        break;
        case Lesser:
            out << "<";
        break;
        case LesserEqual:
            out << "<=";
        break;
        case Greater:
            out << ">";
        break;
        case GreaterEqual:
            out << ">=";
        break;
        default:
            assert (false);
        }
        return out;
    }
};

} // namespace _4ti2_zsolve_
  
#endif
