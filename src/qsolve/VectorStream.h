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

#ifndef _4ti2_qsolve__VectorTStream_
#define _4ti2_qsolve__VectorTStream_

#include "qsolve/Vector.h"
#include <iostream>
#include <fstream>

namespace _4ti2_
{

template <class T>
std::ostream&
operator<<(std::ostream& out, const VectorR<T>& v);

template <class T>
std::istream&
operator>>(std::istream& in, VectorR<T>& v);

template <class T>
void
output(const char* filename, const VectorR<T>& v);

template <class T>
void
output(std::ostream& out, const VectorR<T>& v);

// Same as << but only outputs a projection of the vector.
template <class T>
void
print(std::ostream& out, const VectorR<T>& v, int start, int end);

template <class T>
VectorT<T>*
input_VectorT(const char* filename);

template <class T>
VectorT<T>*
input_VectorT(int dim, const char* filename);

// Definitions of template functions

template <class T>
std::ostream&
operator<<(std::ostream& out, const VectorR<T>& v)
{
    for (Index i = 0; i < v.get_size(); ++i) {
        out.width(2);
        out << v[i] << " ";
    }
    return out;
}

// TODO: Add error checking.
template <class T>
std::istream&
operator>>(std::istream& in, VectorR<T>& v)
{
    T value;
    for (Index i = 0; i < v.get_size(); i++) { 
        in >> value;
        v.set(i,value);
    }
    return in;
}

template <class T>
void
output(const char* filename, const VectorR<T>& v)
{
    std::ofstream file(filename);
    output(file, v);
}

template <class T>
void
output(std::ostream& out, const VectorR<T>& v)
{
    out << v.get_size() << "\n" << v << "\n";
}

template <class T>
void
print(std::ostream& out, const VectorR<T>& v, int start, int end)
{
    assert(start >= 0 && start <= end && end <= v.get_size());
    for (Index i = start; i < end; ++i) {
        out.width(2);
        out << v[i] << " ";
    }
    out << "\n";
}

template <class T>
VectorT<T>*
input_Vector(const char* filename)
{
    std::ifstream file(filename);
    if (!file.good()) { return 0; }
    int n;
    file >> n;
    VectorT<T>* v = new VectorT<T>(n);
    file >> *v;
    if (file.fail() || file.bad()) {
        std::cerr << "INPUT ERROR: Badly formatted file " << filename << ".\n";
        std::cerr << "INPUT ERROR: Check the size.\n";
        std::cerr << "INPUT ERROR: Check there are only integers.";
        std::cerr << std::endl;
        exit(1);
    }
    return v;
}

template <class T>
VectorT<T>*
input_Vector(int dim, const char* filename)
{
    VectorT<T>* v = input_Vector<VectorT<T> >(filename);
    if (v != 0 && v->get_size() != dim) {
        std::cerr << "INPUT ERROR: Incorrect input size in " << filename << ".\n";
        std::cerr << "INPUT ERROR: Size is " << v->get_size();
        std::cerr << ", but should be " << dim << ".\n";
        exit(1);
    }
    return v;
}

} // namespace 4ti2

#endif
