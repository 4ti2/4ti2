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

#ifndef _4ti2__VectorArray_
#define _4ti2__VectorArray_

#include <iostream>
#include <vector>
#include <cassert>

#include "Integer.h"
#include "Vector.hpp"

template <typename T> class VectorArray
{
protected:
    std::vector <T*> m_data;
    size_t m_variables;
    size_t m_vectors;

public:
    VectorArray ()
    {
        m_vectors = m_variables = 0;
    }

    VectorArray (const size_t variables)
    {
        m_variables = variables;
        m_vectors = 0;
    }

    VectorArray (const size_t height, const size_t width)
    {
        m_vectors = height;
        m_variables = width;
        m_data.resize (height);
        for (size_t i = 0; i < height; i++)
            m_data[i] = create_vector <T> (width);
    }

    VectorArray (const VectorArray& other)
    {
        m_vectors = other.m_vectors;
        m_variables = other.m_variables;
        m_data.resize (m_vectors);
        for (size_t i = 0; i < m_vectors; i++)
        {
            m_data[i] = copy_vector (other[i], m_variables);
        }
    }

    ~VectorArray ()
    {
        clear ();
    }

    void clear ()
    {
        for (size_t i = 0; i < m_vectors; i++)
        {
            delete_vector (m_data[i]);
        }
        m_data.clear ();
        m_vectors = 0;
    }

    T* operator[] (const size_t index) const
    {
        assert (index >= 0 && index < m_vectors);

        return m_data [index];
    }
    
    size_t vectors () const
    {
        return m_vectors;

    }

    size_t height () const
    {
        return m_vectors;
    }

    size_t variables () const
    {
        return m_variables;
    }

    size_t width () const
    {
        return m_variables;
    }

    int append_vector (T* vector)
    {
        assert (vector != NULL);

        m_data.push_back (vector);
        m_vectors++;

        assert (m_vectors == m_data.size ());

        return m_vectors-1;
    }

    void remove_unsorted (size_t index)
    {
        delete[] m_data[index];
        m_data[index] = m_data[m_vectors-1];
        m_data.pop_back ();
        m_vectors--;
    }

    void append_negatives ()
    {
        for (int i = m_vectors - 1; i >= 0; i--)
        {
            T* vector = copy_vector (m_data[i], m_variables);
            negate_vector (vector, m_variables);
            append_vector (vector);
        }
    }

    void swap_rows (size_t a, size_t b)
    {
        assert (a < m_vectors);
        assert (b < m_vectors);

        T* tmp = m_data[a];
        m_data[a] = m_data[b];
        m_data[b] = tmp;
    }

    void swap_columns (size_t a, size_t b)
    {
        assert (a < m_variables);
        assert (b < m_variables);

        for (size_t i = 0; i < m_vectors; i++)
            swap_vector (m_data[i], a, b);
    }

    void set_identity (size_t s)
    {
        clear ();

        m_variables = s;
        m_vectors = s;
        m_data.resize (s);
        for (size_t i = 0; i < s; i++)
            m_data[i] = create_unit_vector <T> (s, i);
    }

    T gcd_row (size_t index, size_t start, size_t end)
    {
        T result = m_data[index][start++];
        while (start < end)
            result = gcd (result, m_data[index][start++]);
        return result;
    }
    
    T gcd_column (size_t index, size_t start, size_t end)
    {
        T result = m_data[start++][index];
        while (start < end)
            result = gcd (result, m_data[start++][index]);
        return result;
    }

    void negate_row (size_t index)
    {
        negate_vector (m_data[index]);
    }

    void negate_column (size_t index)
    {
        for (size_t i = 0; i < m_vectors; i++)
        {
            m_data[i][index] = -m_data[i][index];
        }
    }

    void combine_columns (size_t dest, const T& factor, size_t src)
    {
        for (size_t i = 0; i < m_vectors; i++)
        {
            m_data[i][dest] += factor * m_data[i][src];
        }
    }
    
    void combine_rows (size_t dest, const T& factor, size_t src)
    {
        for (size_t i = 0; i < m_variables; i++)
        {
            m_data[dest][i] += factor * m_data[src][i];
        }
    }

    bool check_consistency () const
    {
        if (m_variables == 0)
            return false;

        if (m_vectors != m_data.size ())
            return false;

        for (size_t vec = 0; vec < m_vectors; vec++)
            if (! check_vector_consistency (m_data[vec], m_variables))
                return false;

        return true;
    }

    void save (std::string name) const
    {
        std::ofstream file (name.c_str());
        file << *this;
    }

    template <typename X> friend std::ostream& operator<< (std::ostream& out, const VectorArray <X>& va);
    template <typename X> friend std::istream& operator>> (std::istream& in, VectorArray <X>& va);
};
    
template <typename T> std::ostream& operator<< (std::ostream& out, const VectorArray <T>& va)
{
    out << va.m_vectors << ' ' << va.m_variables << '\n';
    for (size_t i = 0; i < va.m_vectors; i++)
    {
        print_vector <T> (out, va.m_data[i], va.m_variables);
        out << '\n';
    }
    return out;
}

template <typename T> std::istream& operator>> (std::istream& in, VectorArray <T>& va)
{
    va.clear ();
    in >> va.m_vectors >> va.m_variables;
    va.m_data.resize (va.m_vectors);
    for (size_t i = 0; i < va.m_vectors; ++i)
    {
        va.m_data[i] = read_vector <T> (in, va.m_variables);   
    }
    return in;    
}


#endif
