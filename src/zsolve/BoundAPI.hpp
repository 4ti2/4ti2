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

#ifndef _4ti2_zsolve__BoundAPI_
#define _4ti2_zsolve__BoundAPI_

#include "zsolve/VectorArrayAPI.hpp"

namespace _4ti2_zsolve_ {

template <class T>
class BoundAPI : public VectorArrayAPI<T> {
public:
    BoundAPI(int num_rows, int num_cols, bool _is_lower);

    virtual void read(std::istream& in);

protected:
    bool is_lower;
};


template <class T>
BoundAPI<T>::BoundAPI(int num_rows, int num_cols, bool _is_lower)
    : VectorArrayAPI<T>(num_rows, num_cols), is_lower(_is_lower)
{
    if (num_rows != 1) {
        throw IOException("Bounds matrix must have height of 1.");
    }
}

template <class T>
void
BoundAPI<T>::read(std::istream& in)
{
    assert(VectorArrayAPI<T>::data.height() == 1);
    if (!in.good()) { throw IOException("Unreadable istream for bounds."); }
    T v;
    std::string s;
    for (size_t i = 0; i < VectorArrayAPI<T>::data.width(); ++i) {
        in >> v;
        if (!in.fail()) { VectorArrayAPI<T>::data[0][i] = v; }
        else {
            in.clear(); // Clear the error state flags.
            in >> s; 
            if (in.fail()) { throw IOException("Unreadable istream for bounds."); }
            if (s == "*") { VectorArrayAPI<T>::data[0][i] = (is_lower ? 1: -1); }
            else { throw IOException("Unrecognised input for bounds: " + s); }
        }
    }
}

} // namspace _4ti2_

#endif

