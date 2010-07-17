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

#include "zsolve/SignAPI.h"

namespace _4ti2_zsolve_ {

SignAPI::SignAPI(int num_rows, int num_cols)
    : VectorArrayAPI<int>(num_rows, num_cols)
{
    if (num_rows != 1) {
        throw IOException("Sign matrix must have height of 1.");
    }
}

void
SignAPI::read(std::istream& in)
{
    assert(VectorArrayAPI<int>::data.height() == 1);
    if (!in.good()) { throw IOException("Unreadable input stream for sign."); }
    std::string s;
    for (size_t i = 0; i < VectorArrayAPI<int>::data.width(); ++i) {
        in >> s;
        if (in.fail()) { throw IOException("Unreadable input stream for sign."); }
        if (s == "0") { VectorArrayAPI<int>::data[0][i] = 0; }
        else if (s == "1" || s == "+1") { VectorArrayAPI<int>::data[0][i] = 1; }
        else if (s == "-1") { VectorArrayAPI<int>::data[0][i] = -1; }
        else if (s == "2") { VectorArrayAPI<int>::data[0][i] = 2; }
        // Deprecated symbol processing.
        else if (s == "free" || s == "f") { VectorArrayAPI<int>::data[0][i] = 0; deprecated(); }
        else if (s == "hilbert" || s == "h" || s == "+h" || s == "+") { VectorArrayAPI<int>::data[0][i] = 1; deprecated(); }
        else if (s == "-hilbert" || s == "-h" || s == "-") { VectorArrayAPI<int>::data[0][i] = -1; deprecated(); }
        else if (s == "graver" || s == "g") { VectorArrayAPI<int>::data[0][i] = 2; deprecated(); }
        else { throw IOException("Unrecognised input for sign: " + s); }
    }
}

void
SignAPI::deprecated()
{
    std::cerr << "WARNING: Use of deprecated symbols for sign.\n";
    std::cerr << "WARNING: Please use `0' for free variable, `1' for >=, `-1' for <= and `2' for graver variable instead.";
}

} // namspace _4ti2_

