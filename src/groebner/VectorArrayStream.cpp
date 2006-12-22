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

#include "VectorArrayStream.h"
#include "VectorStream.h"
#include <fstream>

using namespace _4ti2_;

std::ostream&
_4ti2_::operator<<(std::ostream& out, const VectorArray& vs)
{
    for (Index i = 0; i < vs.get_number(); ++i)
    {
        out << vs[i] << "\n";
    }
    return out;
}

std::istream&
_4ti2_::operator>>(std::istream& in, VectorArray& vs)
{
    for (Index i = 0; i < vs.get_number(); ++i)
    {
        in >> vs[i];
    }
    return in;
}

void
_4ti2_::output(const char* filename, const VectorArray& vs)
{
    std::ofstream file(filename);
    output(file, vs);
}

void
_4ti2_::output(std::ostream& out, const VectorArray& vs)
{
    out << vs.get_number() << " " << vs.get_size() << "\n" << vs;
}

void
_4ti2_::print(std::ostream& out, const VectorArray& vs, int start, int end)
{
    for (Index i = 0; i < vs.get_number(); ++i)
    {
        print(out, vs[i], start, end);
    }
}

VectorArray*
_4ti2_::input_VectorArray(const char* filename)
{
    std::ifstream file(filename);
    if (!file.good()) { return 0; }
    int m,n;
    file >> m >> n;
    VectorArray* vs_ptr = new VectorArray(m,n);
    file >> *vs_ptr;
    if (file.fail() || file.bad())
    {
        std::cerr << "INPUT ERROR: Badly formatted file " << filename << ".\n";
        std::cerr << "INPUT ERROR: Check the number of rows and columns.\n";
        std::cerr << "INPUT ERROR: Check there are only integers.";
        std::cerr << std::endl;
        exit(1);
    }
    return vs_ptr;
}

VectorArray*
_4ti2_::input_VectorArray(int dim, const char* filename)
{
    VectorArray* vs = input_VectorArray(filename);
    if (vs != 0 && vs->get_size() != dim)
    {
        std::cerr << "INPUT ERROR: Incorrect input size in " << filename << ".\n";
        std::cerr << "INPUT ERROR: Size is " << vs->get_size();
        std::cerr << ", but should be " << dim << ".\n";
        exit(1);
    }

    return vs;
}
