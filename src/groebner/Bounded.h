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

#ifndef _4ti2__Bounded_
#define _4ti2__Bounded_

#include "VectorArray.h"
#include "BitSet.h"

namespace _4ti2_ {

void bounded(   const VectorArray& matrix,
                const VectorArray& lattice,
                const BitSet& urs,
                BitSet& bounded,
                Vector& grading,
                BitSet& unbounded,
                Vector& ray);

bool bounded(   const VectorArray& matrix,
                const VectorArray& lattice,
                const BitSet& urs,
                const VectorArray& cost,
                const BitSet& bnd,
                const Vector& grading,
                const BitSet& unbnd,
                const Vector& ray,
                BitSet& cost_unbnd);

void bounded_projection(
                const VectorArray& matrix,
                const VectorArray& lattice,
                const BitSet& urs,
                const Vector& rhs,
                BitSet& proj);

void lp_weight_l1(
                const VectorArray& lattice,
                const BitSet& urs,
                const Vector& rhs,
                Vector& weight);

void lp_weight_l2(
                const VectorArray& lattice,
                const BitSet& urs,
                const Vector& rhs,
                Vector& weight);

bool lp_feasible(
                const VectorArray& lattice,
                const Vector& rhs);

bool ip_feasible(
                const VectorArray& lattice,
                const Vector& rhs);

// Solves min{cx : Ax = rhs, x >= 0}
int lp_solve(
                const VectorArray& matrix,
                const Vector& rhs,
                const Vector& cost,
                const BitSet& urs,
                BitSet& basic,
                RationalType& objective);

} // namespace _4ti2_

#endif
