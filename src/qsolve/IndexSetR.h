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

#ifndef _4ti2_qsolve__IndexSetR_
#define _4ti2_qsolve__IndexSetR_

#include "qsolve/Size.h"
#include "qsolve/Index.h"

#include <cassert>

namespace _4ti2_
{

class IndexSetR
{
public:
    IndexSetR(Index _l, Index _u) : l(_l), u(_u) {}
    IndexSetR(const IndexSetR& is) : l(is.l), u(is.u) {}
    IndexSetR& operator=(const IndexSetR& is) { l = is.l; u = is.u; return *this; }

    typedef Index Iter;
/*
    class Iter {
    public:
        operator Index() const { return i; }
        Index operator*() const { return i; }
        void operator++() { ++i; }
        void operator=(Iter it) { i = it.i; }
        bool operator!=(Iter it) { return i != it.i; }
    private:
        Iter(Index _i) : i(_i) {}
        Index i;
        friend class IndexSetR;
    };
*/

    Size count() const;

    Iter begin() const { return l; }
    Iter end() const { return u; }

protected:
    IndexSetR();

    // Represents the set of integers [l,u).
    Iter l;
    Iter u;
};

inline
Size
IndexSetR::count() const
{
    return u-l;
}

} // namespace _4ti2_

#endif
