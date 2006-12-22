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

#include "LongDenseIndexSetStream.h"

using namespace _4ti2_;

std::ostream&
_4ti2_::operator<<(std::ostream& out, const LongDenseIndexSet& b)
{
    for (Index i = 0; i < b.get_size(); ++i)
    {
        out.width(2);
        out << b[i] << " ";
    }
    return out;
}

std::istream&
_4ti2_::operator>>(std::istream& in, LongDenseIndexSet& b)
{
    for (Index i = 0; i < b.get_size(); ++i)
    {
        bool bit;
        in >> bit;
        if (bit == true) b.set(i);
        else b.unset(i);
    }
    return in;
}

void
_4ti2_::output(const char* filename, const LongDenseIndexSet& bs)
{
    std::ofstream file(filename);
    output(file, bs);
}

void
_4ti2_::output(std::ostream& out, const LongDenseIndexSet& bs)
{
    out << bs.get_size() << "\n" << bs << "\n";
}

LongDenseIndexSet*
_4ti2_::input_LongDenseIndexSet(const char* filename)
{
    std::ifstream file(filename);
    if (!file.good()) { return 0; }
    int n;
    file >> n;
    LongDenseIndexSet* bs = new LongDenseIndexSet(n);
    file >> *bs;
    if (file.fail() || file.bad())
    {
        std::cerr << "ERROR: Badly formatted file " << filename << ".\n";
        std::cerr << "ERROR: Check the size.\n";
        std::cerr << "ERROR: Check there are 0 or 1 entries.";
        std::cerr << std::endl;
        exit(1);
    }
    return bs;
}
