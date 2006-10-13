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

#ifndef _4ti2__BinomialFactory_
#define _4ti2__BinomialFactory_

#include "Binomial.h"
#include "Vector.h"
#include "VectorArray.h"
#include "DataType.h"
#include "Statistics.h"
#include "BitSet.h"
#include "Feasible.h"
#include "Weight.h"
#include "Grading.h"
#include "Permutation.h"
#include "VectorStream.h"
#include "Index.h"
#include "Size.h"
#include "Bounded.h"
#include "Filter.h"
#include "BinomialArray.h"
#include "BinomialSet.h"
#include <iostream>

namespace _4ti2_
{

class BinomialFactory
{
public:
    BinomialFactory(Feasible& feasible,
                    const VectorArray& cost);
    BinomialFactory(Feasible& feasible,
                    const VectorArray& cost,
                    const BitSet& sat);
    ~BinomialFactory();

    void convert(const Binomial& b, Vector& v) const;
    void convert(const Vector& v, Binomial& b) const;
    void convert(const VectorArray& vs, BinomialCollection& bc, bool orientate=true) const;
    void convert(const BinomialArray& bs, VectorArray& vs) const;
    void convert(const BinomialSet& bs, VectorArray& vs) const;

    void add_weight(
                    const Vector& weight,
                    IntegerType max_weight);

protected:
    void initialise(
                    int dim,
                    const VectorArray& lattice,
                    const VectorArray& cost,
                    const BitSet& urs,
                    const BitSet& bnd,
                    const BitSet& unbnd,
                    const Vector& grading,
                    const VectorArray* weights = 0,
                    const Vector* max_weights = 0,
                    const Vector* rhs = 0);

    void set_weights(
                    const VectorArray* weights,
                    const Vector* max);
    void set_truncated(
                    const VectorArray& lattice,
                    const Vector* _rhs);

protected:
    void check_cost(Feasible feasible, VectorArray& cost);
    void initialise_permutation(
                    const BitSet& bnd_mask,
                    const BitSet& urs_mask);
    Permutation* permutation;
    VectorArray* costs;
    BitSet* orig_bnd;
};

} // namespace _4ti2_

#endif
