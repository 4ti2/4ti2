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

#include "BitSetStream.h"

using namespace _4ti2_;

BitSet*
_4ti2_::input_BitSet(const char* filename)
{
    std::ifstream file(filename);
    if (!file.good()) { return 0; }
    int n;
    file >> n;
    BitSet* bs = new BitSet(n);
    file >> *bs;
    if (file.fail() || file.bad())
    {
        std::cerr << "INPUT ERROR: Badly formatted file " << filename << ".\n";
        std::cerr << "INPUT ERROR: Check the size.\n";
        std::cerr << "INPUT ERROR: Check there are 0 or 1 entries.";
        std::cerr << std::endl;
        exit(1);
    }
    return bs;
}

BitSet*
_4ti2_::input_BitSet(int dim, const char* filename)
{
    BitSet* bs = input_BitSet(filename);
    if (bs != 0 && bs->get_size() != dim)
    {
        std::cerr << "INPUT ERROR: Incorrect input size in " << filename << ".\n";
        std::cerr << "INPUT ERROR: Size is " << bs->get_size();
        std::cerr << ", but should be " << dim << ".\n";
        exit(1);
    }
    return bs;
}
