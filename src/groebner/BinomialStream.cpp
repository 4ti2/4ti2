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

#include "BinomialStream.h"

using namespace _4ti2_;

std::ostream&
_4ti2_::operator<<(std::ostream& out, const Binomial& b)
{
    // Output the bounded and restricted in sign components.
    for (Index i = 0; i < b.bnd_end; ++i)
    {
        out.width(2);
        out << " " << b[i];
    }
    out << "| ";
    // Output the unbounded and restricted in sign components.
    for (Index i = b.bnd_end; i < b.rs_end; ++i)
    {
        out.width(2);
        out << " " << b[i];
    }
    out << "| ";
    // Output the unrestricted in sign components.
    for (Index i = b.rs_end; i < b.urs_end; ++i)
    {
        out.width(2);
        out << " " << b[i];
    }
    out << "| ";
    // Output the cost components.
    for (Index i = b.cost_start; i < b.cost_end; ++i)
    {
        out.width(2);
        out << " " << b[i];
    }
    out << "| ";
    // Output the rest of the binomial.
    for (Index i = b.cost_end; i < b.size; ++i)
    {
        out.width(2);
        out << " " << b[i];
    }
    return out;
}
