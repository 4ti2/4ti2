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

#ifndef _4ti2__Vector_
#define _4ti2__Vector_

#include <cassert>
#include <iostream>

#include "zsolve/Integer.h"

template <typename T> T* create_vector (size_t size)
{
    assert (size > 0);

    return new T [size];
}

template <typename T> T* create_zero_vector (size_t size)
{
    assert (size > 0);

    T* result = new T [size];
    for (size_t i = 0; i < size; i++)
        result[i] = 0;
    return result;
}

template <typename T> T* create_unit_vector (size_t size, size_t index)
{
    assert (size > 0);

    T* result = new T [size];
    for (size_t i = 0; i < size; i++)
        result[i] = 0;
    result[index] = 1;
    return result;
}

template <typename T> T* copy_vector (T* other, size_t size)
{
    assert (size > 0);
    assert (other != NULL);

    T* result = new T [size];
    for (size_t i = 0; i < size; i++)
        result[i] = other[i];
    return result;
}

template <typename T> void delete_vector (T* vector)
{
    assert (vector != NULL);

    delete[] vector;
}

template <typename T> T* read_vector (std::istream& in, size_t size)
{
    assert (size > 0);

    T* result = create_vector <T> (size);
    for (size_t i = 0; i < size; i++)
    {
        in >> result [i];
    }
    return result;
}

template <typename T> std::ostream& print_vector (std::ostream& out, T* vector, size_t size)
{
    assert (vector != NULL);
    assert (size > 0);

    for (size_t i = 0; i < size; i++)
    {
        if (i > 0)
            out << " ";
        out << vector[i];
    }
    return out;
}

template <typename T> void add_vectors (T* a, T* b, T* r, size_t size)
{
    assert (a != NULL);
    assert (b != NULL);
    assert (r != NULL);
    assert (size > 0);

    for (size_t i = 0; i < size; i++)
    {
        r[i] = a[i] + b[i];
    }
}

template <typename T> void negate_vector (T* v, size_t size)
{
    assert (v != NULL);
    assert (size > 0);

    for (size_t i = 0; i < size; i++)
    {
        v[i] = -v[i];
    }
}

template <typename T> T norm_vector (T* v, size_t size)
{
    assert (v != NULL);

    T result = 0;
    for (size_t i = 0; i < size; i++)
    {
        result += abs (v[i]);
    }
    return result;
}

template <typename T> T gcd_vector (T* v, size_t size)
{
    assert (size > 0);

    T result =  v[0];

    for (size_t i = 1; i < size; i++)
        result = gcd (v[i], result);

    return result;
}

template <typename T> void swap_vector (T* v, size_t a, size_t b)
{
    assert (v != NULL);

    T tmp = v[a];
    v[a] = v[b];
    v[b] = tmp;
}

template <typename T> bool is_zero_vector (T* v, size_t size)
{
    assert (v != NULL);
    assert (size > 0);

    for (size_t i = 0; i < size; i++)
    {
        if (v[i] != 0)
            return false;
    }
    return true;
}

template <typename T> bool check_vector_consistency (T* v, size_t size)
{
    if (v == NULL || size == 0)
        return false;

    T result = 0;
    for (size_t i = 0; i < size; i++)
    {
        result += abs (v[i]);
    }

    return true;
}

template <typename T> int lex_compare_vector_with_negative (T* v, size_t size)
{
    size_t i = 0;
    while (i < size && v[i] == 0)
        i++;
    if (i == size)
        return 0;
    else if (v[i] < 0)
        return -1;
    else if (v[i] > 0)
        return 1;
    else
        return 0;
}


#endif
