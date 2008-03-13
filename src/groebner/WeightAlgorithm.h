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

#ifndef _4ti2_groebner__WeightAlgorithm_
#define _4ti2_groebner__WeightAlgorithm_

#include "groebner/VectorArray.h"
#include "groebner/BitSet.h"
#include "groebner/Weight.h"

namespace _4ti2_
{

// TODO: Allow for urs varaibles.
class WeightAlgorithm
{
public:
    static bool check_weights(
                    const VectorArray& matrix,
                    const VectorArray& lattice,
                    const BitSet& urs,
                    VectorArray& weights);

    static bool get_weights(
                    const VectorArray& matrix,
                    const VectorArray& lattice,
                    const BitSet& urs,
                    VectorArray& weights);

    static void strip_weights(
                VectorArray* weights,
                Weight* max_weight,
                const BitSet& urs);

private:
    static bool is_candidate(
                    const Vector& v,
                    const BitSet& urs,
                    const BitSet& mask);
    static bool violates_urs(const Vector& v, const BitSet& mask);
    static int positive_count(const Vector& v, const BitSet& mask);
    static void update_mask(BitSet& mask, const Vector& v);
    static bool get_weights(
                    const VectorArray& matrix,
                    const BitSet& urs,
                    BitSet& mask,
                    VectorArray& weights);
};

} // namespace _4ti2_

#endif
