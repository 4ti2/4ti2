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

#ifndef _4ti2_groebner__LongDenseIndexSetStream_
#define _4ti2_groebner__LongDenseIndexSetStream_

#include "groebner/LongDenseIndexSetStream.h"
#include "groebner/ShortDenseIndexSetStream.h"

#include <iostream>
#include <fstream>
#include "groebner/LongDenseIndexSet.h"

namespace _4ti2_
{

std::ostream&
operator<<(std::ostream& out, const LongDenseIndexSet& b);

std::istream&
operator>>(std::istream& in, LongDenseIndexSet& b);

void
output(const char* filename, const LongDenseIndexSet& bs);

void
output(std::ostream& out, const LongDenseIndexSet& bs);

LongDenseIndexSet*
input_LongDenseIndexSet(const char* filename);

} // namespace _4ti2_

#endif
