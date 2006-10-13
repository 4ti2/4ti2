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

#include "Binomial.h"
#include "WeightAlgorithm.h"
#include <iostream>
#include "BitSetStream.h"
#include "VectorArrayStream.h"
#include "VectorStream.h"
#include "BinomialStream.h"
#include "Globals.h"

//#define DEBUG_4ti2(X) X
#include "Debug.h"

using namespace _4ti2_;

Index Binomial::rs_end = 0;
Index Binomial::bnd_end = 0;
Index Binomial::urs_end = 0;
Index Binomial::cost_start = 0;
Index Binomial::cost_end = 0;
Size Binomial::size = 0;
VectorArray* Binomial::weights = 0;
Weight* Binomial::max_weights = 0;
Vector* Binomial::rhs = 0;
VectorArray* Binomial::lattice = 0;
Grading* Binomial::grading = 0;

bool
Binomial::truncated(const Binomial& b)
{
    if (rhs != 0)
    {
        assert(lattice != 0);
        Vector v(rhs->get_size());
        for (int i = 0; i < bnd_end; ++i)
        {
            if (b.data[i] > 0) { v[i] = (*rhs)[i] - b.data[i]; }
            else { v[i] = (*rhs)[i]; }
        }
        DEBUG_4ti2(*out << "Binomial:\n" << b << "\n";)
        DEBUG_4ti2(*out << "RHS:\n" << *rhs << "\n";)
        DEBUG_4ti2(*out << "Vector:\n" << v << "\n";)
        if (Globals::truncation == Globals::IP)
        {
            return !(ip_feasible(*lattice, v));
        }
        else
        {
            return !(lp_feasible(*lattice, v));
        }
    }
    return false;
}
