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

#ifndef _4ti2_qsolve__VectorArrayStream_
#define _4ti2_qsolve__VectorArrayStream_

#include "qsolve/VectorArray.h"
#include "qsolve/VectorStream.h"
#include <iostream>
#include <fstream>

namespace _4ti2_
{

template <class T>
std::ostream&
operator<<(std::ostream& out, const VectorArrayT<T>& vs);

template <class T>
std::istream&
operator>>(std::istream& in, VectorArrayT<T>& vs);

template <class T>
void
output(const char* filename, const VectorArrayT<T>& vs);

template <class T>
void
output(std::ostream& out, const VectorArrayT<T>& vs);

// Same as << but only outputs a projection of the vector array.
template <class T>
void
print(std::ostream& out, const VectorArrayT<T>& vs, int start, int end);

template <class T>
VectorArrayT<T>*
input_VectorArray(const char* filename);

template <class T>
VectorArrayT<T>*
input_VectorArray(int dim, const char* filename);

// Template function definitions.

template <class T>
std::ostream&
operator<<(std::ostream& out, const VectorArrayT<T>& vs)
{
    for (Index i = 0; i < vs.get_number(); ++i) {
        out << vs[i] << "\n";
    }
    return out;
}

template <class T>
std::istream&
operator>>(std::istream& in, VectorArrayT<T>& vs)
{
    for (Index i = 0; i < vs.get_number(); ++i) {
        in >> vs[i];
    }
    return in;
}

template <class T>
void
output(const char* filename, const VectorArrayT<T>& vs)
{
    std::ofstream file(filename);
    output(file, vs);
}

template <class T>
void
output(std::ostream& out, const VectorArrayT<T>& vs)
{
    out << vs.get_number() << " " << vs.get_size() << "\n" << vs;
}

template <class T>
void
print(std::ostream& out, const VectorArrayT<T>& vs, int start, int end)
{
    for (Index i = 0; i < vs.get_number(); ++i) {
        print(out, vs[i], start, end);
    }
}

template <class T>
VectorArrayT<T>*
input_VectorArray(const char* filename)
{
    std::ifstream file(filename);
    if (!file.good()) { return 0; }
    int m,n;
    m = n = 0;
    file >> m >> n;
    VectorArrayT<T>* vs_ptr = new VectorArrayT<T>(m,n);
    file >> *vs_ptr;
    if (file.fail() || file.bad()) {
        std::cerr << "INPUT ERROR: Badly formatted file " << filename << ".\n";
        std::cerr << "INPUT ERROR: Check the number of rows and columns.\n";
        std::cerr << "INPUT ERROR: Check there are only integers.";
        std::cerr << std::endl;
        exit(1);
    }
    return vs_ptr;
}

template <class T>
VectorArrayT<T>*
input_VectorArray(int dim, const char* filename)
{
    VectorArrayT<T>* vs = input_VectorArray<T>(filename);
    if (vs != 0 && vs->get_size() != dim) {
        std::cerr << "INPUT ERROR: Incorrect input size in " << filename << ".\n";
        std::cerr << "INPUT ERROR: Size is " << vs->get_size();
        std::cerr << ", but should be " << dim << ".\n";
        exit(1);
    }

    return vs;
}

}

#endif
