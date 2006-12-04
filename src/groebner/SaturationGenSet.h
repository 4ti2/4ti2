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

#ifndef _4ti2__SaturationGenSet_
#define _4ti2__SaturationGenSet_

#include "VectorArray.h"
#include "BitSet.h"
#include "Feasible.h"

namespace _4ti2_
{

class SaturationGenSet
{
public:
    SaturationGenSet();
    virtual ~SaturationGenSet();

    virtual void compute(
                    Feasible& feasible,
                    VectorArray& gens,
                    BitSet& sat,
                    bool minimal = true);

protected:
    virtual void compute_bounded(
                    Feasible& feasible,
                    VectorArray& gens,
                    BitSet& sat,
                    bool minimal = true);

    bool is_saturated(
                    const BitSet& sat,
                    const BitSet& urs);

    int compute_saturations(
                    const VectorArray& gens,
                    const BitSet& sat,
                    const BitSet& urs,
                    VectorArray& useful_gens);

    int saturate(
                    const VectorArray& gens,
                    BitSet& sat,
                    const BitSet& urs);
    int saturate(
                    const VectorArray& gens,
                    BitSet& sat,
                    const BitSet& urs,
                    VectorArray& useful);

    int next_saturation(
                    const VectorArray& gens,
                    BitSet& sat,
                    const BitSet& urs);

    int add_support(
                    const Vector& p,
                    BitSet& sat,
                    const BitSet& urs);

    void support_count(
                    const Vector& p,
                    BitSet& sat,
                    const BitSet& urs,
                    int& pos,
                    int& neg);

    void saturate_zero_columns(
                    const VectorArray& gens,
                    BitSet& sat,
                    const BitSet& urs);

    bool is_column_zero(
                    const VectorArray& gens,
                    int col);
};
    
}

#endif
