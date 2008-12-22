/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2008 4ti2 team.
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

#include "zsolve/RelAPI.h"

namespace _4ti2_zsolve_ {

RelAPI::RelAPI(int num_rows, int num_cols)
    : VectorArrayAPI<int>(num_rows, num_cols)
{
    if (num_rows != 1) {
        throw IOException("Relations matrix must have height of 1.");
    }
}

void
RelAPI::read(std::istream& in)
{
    assert(VectorArrayAPI<int>::data.height() == 1);
    if (!in.good()) { throw IOException("Unreadable istream for relations."); }
    std::string s;
    for (size_t i = 0; i < VectorArrayAPI<int>::data.width(); ++i) {
        in >> s;
        if (in.fail()) { throw IOException("Unreadable istream for relations."); }
        if (s == "<" || s == "-1") { VectorArrayAPI<int>::data[0][i] = -1; }
        else if (s == ">" || s == "1") { VectorArrayAPI<int>::data[0][i] = 1; }
        else if (s == "=" || s == "0") { VectorArrayAPI<int>::data[0][i] = 0; }
        else { throw IOException("Unrecognised input for relations: " + s); }
    }
}

} // namspace _4ti2_

