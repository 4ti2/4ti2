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

#include "EliminateAlgorithm.h"
#include "HermiteAlgorithm.h"
#include "Globals.h"

#include "VectorArrayStream.h"
#include "VectorStream.h"
#include "BitSetStream.h"

//#define DEBUG_4ti2(X) X
#include "Debug.h"

using namespace _4ti2_;

void
_4ti2_::eliminate(VectorArray& vs, const BitSet& bs)
{
    BitSet tmp(bs);
    tmp.set_complement();
    Index rows = upper_triangle(vs, tmp);
    vs.remove(0, rows);
}

void
_4ti2_::eliminate(VectorArray& vs, Index num_vars)
{
    Index rows = upper_triangle(vs, vs.get_number(), num_vars);
    DEBUG_4ti2(*out << "After Hermite:\n" << vs;)
    vs.remove(0, rows);
}

