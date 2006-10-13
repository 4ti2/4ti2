/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.

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

#include "VectorStream.h"

using namespace _4ti2_;

std::ostream&
_4ti2_::operator<<(std::ostream& out, const Vector& v)
{
    for (Index i = 0; i < v.get_size(); ++i)
    {
        out.width(2);
        out << v[i] << " ";
    }
    return out;
}

// TODO: Add error checking.
std::istream&
_4ti2_::operator>>(std::istream& in, Vector& v)
{
    for (Index i = 0; i < v.get_size(); i++) { in >> v[i]; }
    return in;
}

void
_4ti2_::output(const char* filename, const Vector& v)
{
    std::ofstream file(filename);
    output(file, v);
}

void
_4ti2_::output(std::ostream& out, const Vector& v)
{
    out << v.get_size() << "\n" << v << "\n";
}

void
_4ti2_::print(std::ostream& out, const Vector& v, int start, int end)
{
    assert(start >= 0 && start <= end && end <= v.get_size());
    for (Index i = start; i < end; ++i)
    {
        out.width(2);
        out << v[i] << " ";
    }
    out << "\n";
}

Vector*
_4ti2_::input_Vector(const char* filename)
{
    std::ifstream file(filename);
    if (!file.good()) { return 0; }
    int n;
    file >> n;
    Vector* v = new Vector(n);
    file >> *v;
    if (file.fail() || file.bad())
    {
        std::cerr << "INPUT ERROR: Badly formatted file " << filename << ".\n";
        std::cerr << "INPUT ERROR: Check the size.\n";
        std::cerr << "INPUT ERROR: Check there are only integers.";
        std::cerr << std::endl;
        exit(1);
    }
    return v;
}

Vector*
_4ti2_::input_Vector(int dim, const char* filename)
{
    Vector* v = input_Vector(filename);
    if (v != 0 && v->get_size() != dim)
    {
        std::cerr << "INPUT ERROR: Incorrect input size in " << filename << ".\n";
        std::cerr << "INPUT ERROR: Size is " << v->get_size();
        std::cerr << ", but should be " << dim << ".\n";
        exit(1);
    }
    return v;
}
