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

#include "BinomialSetStream.h"
#include "BinomialStream.h"

using namespace _4ti2_;

std::ostream&
_4ti2_::operator<<(std::ostream& out, const BinomialSet& bs)
{
    for (Index i = 0; i < (Index) bs.get_number(); ++i)
    {
        out << "(" << i << ") " << bs[i] << "\n";
    }
    return out;
}
