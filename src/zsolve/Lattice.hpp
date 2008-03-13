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

#ifndef _4ti2__Lattice_
#define _4ti2__Lattice_

#include <iostream>
#include <cassert>

#include "zsolve/Exception.h"
#include "zsolve/LinearSystem.hpp"

template <typename T> class Lattice : public VectorArray <T>, public VariableProperties <T>
{
public:
    Lattice (VariableProperties<T>* properties) : VectorArray <T> (properties->variables ()), VariableProperties <T> (properties) { }

    Lattice (VectorArray <T> * vectors, bool free, const T& lower, const T& upper) : VectorArray <T> (*vectors), VariableProperties <T> (vectors->variables (), free, lower, upper) { }

    Lattice (Lattice <T>& other) : VectorArray <T> ((VectorArray <T> &)other), VariableProperties <T> (&other) {}

    bool check_consistency () const
    {
        if (!VectorArray <T> :: check_consistency ())
            return false;

        if (!VariableProperties <T> :: check_consistency ())
            return false;

        return true;
    }
 
    void sort_columns ()
    {
        for (size_t i = 0; i < VectorArray <T> :: m_variables; i++)
        {
            size_t k = i;
            for (size_t j = i+1; j < VectorArray <T> :: m_variables; j++)
            {
                if (VariableProperties <T> :: get_variable (j).compare (VariableProperties <T> :: get_variable (k)) < 0)
                    k = j;
            }
            swap_columns (i,k);
        } 
    }

    size_t variables ()
    {
        return VectorArray <T> :: m_variables;
    }

    void filter_bounds (size_t current)
    {
        for (size_t i = 0; i < VectorArray <T> :: m_vectors; i++)
        {
            if (! VariableProperties <T> :: m_variable_properties[current]->check_bounds (VectorArray <T> :: m_data[i][current]))
            {
                //std::cout << "Filtering [";
                //print_vector (std::cout, VectorArray <T> :: m_data[i], VectorArray <T> :: m_variables);
                //std::cout << "]" << std::endl;
                VectorArray <T> :: remove_unsorted (i);
                i--;
            }
        }
    }

    void swap_columns (size_t a, size_t b)
    {
        VectorArray <T> :: swap_columns (a, b);
        VariableProperties <T> :: swap (a, b);
    }

    void reduce_gaussian ()
    {
        for (size_t column = 0; column < VectorArray <T> :: m_variables; column++)
        {
            int current_index;
            T current_value;
            int best_index = column;
            T best_value = VectorArray <T> :: gcd_column (column, column, VectorArray <T> :: m_vectors);
            for (current_index = column+1; current_index < (int) VectorArray <T> :: m_variables; current_index++)
            {
                current_value = VectorArray <T> :: gcd_column (current_index, column, VectorArray <T> :: m_vectors);
                if (current_value < best_value)
                {
                    best_index = current_index;
                    best_value = current_value;
                }
            }
            swap_columns (column, best_index);

            while (true)
            {
                best_index = -1;
                best_value = 0;
                for (current_index = column; current_index < (int) VectorArray <T> :: m_vectors; current_index++)
                {
                    current_value = abs (VectorArray <T> :: m_data[current_index][column]);
                    if ( current_value > 0 && (best_index < 0 || current_value < best_index))
                    {
                        best_index = current_index;
                        best_value = current_value;
                    }
                }

                if (best_index < 0)
                    return;

                VectorArray <T> :: swap_rows (column, best_index);

                bool repeat = false;
                for (size_t i = 0; i < VectorArray <T> :: m_vectors; i++)
                {
                    if (i == column)
                        continue;
                    T factor = - VectorArray <T> :: m_data[i][column] / VectorArray <T> :: m_data[column][column];
                    if (factor != 0)
                    {
                        VectorArray <T> :: combine_rows (i, factor, column);
                        repeat = true;
                    }
                }
                if (!repeat)
                    break;
            }
        }

        for (size_t i = 0; i < VectorArray <T> :: m_vectors; i++)
        {
            if (is_zero_vector (VectorArray <T> :: m_data[i], VectorArray <T> :: m_variables))
            {
                VectorArray <T> :: remove_unsorted (i);
                i--;
            }
        }

    }

    size_t get_result_variables ()
    {
        size_t result = 0;
        for (size_t i = 0; i < VectorArray <T> :: m_variables; i++)
            if (VariableProperties <T> :: get_variable (i).column () >= 0)
                result++;
        return result;
    }

    int get_splitter ()
    {
        for (size_t i = 0; i < VectorArray <T> :: m_variables; i++)
            if (VariableProperties <T> :: get_variable (i).column () == -2)
                return i;
        return -1;
    }

    
    template <typename X> friend std::ostream& operator<< (std::ostream& out, Lattice <X>& lattice);
};

template <typename T> std::ostream& operator<< (std::ostream& out, Lattice <T>& lattice)
{
    size_t vars = lattice.variables ();
    size_t rels = lattice.vectors ();
    size_t* space = new size_t[vars];

    for (size_t i = 0; i < vars; i++)
    {
        VariableProperty <T> & var = lattice.get_variable (i);
        space[i] = max (var.upper_space (), var.lower_space ());
        for (size_t j = 0; j < rels; j++)
        {
            space[i] = max <size_t> (space[i], integer_space(lattice[j][i]));
        }
    }

    // print variables lines
    for (size_t i = 0; i < vars; i++)
    {
        VariableProperty <T> & var = lattice.get_variable (i);
        if (i > 0)
            out << " ";
        for (int j = space[i] - var.upper_space (); j > 0; j--)
            out << " ";
        var.upper (out);
    }
    out << "\n";
    for (size_t i = 0; i < vars; i++)
    {
        VariableProperty <T> & var = lattice.get_variable (i);
        if (i > 0)
            out << " ";
        for (int j = space[i] - var.lower_space (); j > 0; j--)
            out << " ";
        var.lower (out);
    }
    out << "\n";
    for (size_t i = 0; i < vars; i++)
    {
        VariableProperty <T> & var = lattice.get_variable (i);
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
            const T value = lattice[i][j];
            for (int k = space[j] - integer_space (value); k > 0; k--)
                out << " ";
            out << value;
        }
    }

    out << "\n" << std::flush;

    delete[] space;

    return out;
}

template <typename T> Lattice <T>* generate_lattice (LinearSystem <T>* system)
{
    assert (system->is_homogeneous ());
    assert (system->is_equality_system ());

    //std::cerr << "generate_lattice" << std::endl;

    size_t n,e;
    size_t identities = 0;
    VectorArray <T> H (system->matrix ());
    Lattice <T>* result = new Lattice <T> (system);
    
    n = system->variables ();
    e = system->relations ();
    for (int i = n-1; i >= 0; i--)
    {
        bool flag = false;
        for (size_t j = 0; j < identities; j++)
        {
            if (H[j+e-identities][i] != 0)
            {
                flag = true;
                break;
            }
        }
        if (flag)
            continue;
        int identity_row = -1;
        for (size_t j = 0; j < e-identities; j++)
        {
            if (identity_row == -1 && abs (H[j][i]) == 1)
                identity_row = j;
            else if (H[j][i] != 0)
            {
                identity_row = -2;
                break;
            }
        }
        if (identity_row >= 0)
        {
            H.swap_rows (identity_row, e-identities-1);
            H.swap_columns (i, n-identities-1);
            result->swap_columns (i, n-identities-1);
            identities++;
        }
    }

    VectorArray <T> C (H);

    n -= identities;
    e -= identities;

    //std::cerr << "removed " << identities << " identities" << std::endl;
    /*for (int i = 0; i < e; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cerr << H[i][j] << " ";
        }
        std::cerr << std::endl;
    }*/

    if (n == 0)
        return result;

    VectorArray <T> I;
    I.set_identity (n);

    size_t m = min (n, e);
    size_t i;
    for (i = 0; i < m; i++)
    {
        size_t best_index = i;
        T best_value = 0;
        size_t current_index;
        T current_value;

        for (current_index = i; current_index < e; current_index++)
        {
            current_value = H.gcd_row (current_index, i, n);
            if (current_value != 0 && (best_value == 0 || current_value < best_value))
            {
                best_index = current_index;
                best_value = current_value;
            }
        }
        if (best_value == 0)
            break;
        
        H.swap_rows (best_index, i);
        
        //std::cerr << "After swapping Best row to " << i << "\n\n" << H << std::endl;

        bool reduced;
        do
        {
            reduced = false;

            // find a reducer
            best_index = i;
            best_value = 0;
            for (current_index = i; current_index < n; current_index++)
            {
                current_value = abs (H[i][current_index]);
                if (current_value != 0 && (best_value <= 0 || current_value < best_value))
                {
                    best_index = current_index;
                    best_value = current_value;
                }
            }
            assert (best_value != 0);
        
            //std::cerr << "Reducer is column " << best_index << " with value " << H[i][best_index] << std::endl;

            // negative?
            if (H[i][best_index] < 0)
            {
                H.negate_column (best_index);
                I.negate_column (best_index);
            }

            // reduce others
            for (size_t j = 0; j < n; j++)
            {
                if (j != best_index)
                {
                    //if (H[i][j] % best_value != 0 && j > i)
                    //    throw CalcException ("Linear system cannot be transformed into lattice.");
                    // BUG bbb.mat needs this loop, because integer gaussian may have modulus on division!
                    T factor = - H[i][j] / best_value;
                    if (factor != 0)
                    {
                        H.combine_columns (j, factor, best_index);
                        I.combine_columns (j, factor, best_index);
                        reduced = true;
                    }
                }
            }

            // swap
            H.swap_columns (i, best_index);
            I.swap_columns (i, best_index);
        }
        while (reduced);
        
        //std::cerr << "After reduction and column swap:\n" << H << std::endl;
    }

    // fill vector array
    while (i < n)
    {
        T* vector = create_vector <T> (H.variables ());
        for (size_t j = 0; j < n; j++)
        {
            vector[j] = I[j][i];
        }
        for (size_t j = 0; j < identities; j++) 
        {
            T factor = 0;
            for (size_t k = 0; k < n; k++)
            {
                T temp = vector[k]*C[e+j][k];
                factor -= temp;
            }
            vector[n+j] = factor * H[e+j][n+j];
        }
        result->append_vector (vector);
        i++;
    }

    result->sort_columns ();

    return result;
}

#endif
