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

#ifndef _4ti2__RaySupportAlgorithm_
#define _4ti2__RaySupportAlgorithm_

#include "groebner/RayImplementation.h"
#include "groebner/VectorArray.h"

namespace _4ti2_
{

template <class IndexSet>
class RaySupportAlgorithm : public RayImplementation<IndexSet>
{
public:
    RaySupportAlgorithm();
    virtual ~RaySupportAlgorithm();

    virtual IndexSet compute(
                    VectorArray& matrix,
                    VectorArray& vs,
                    const IndexSet& rs);
    virtual IndexSet compute(
                    VectorArray& matrix,
                    VectorArray& vs,
                    std::vector<IndexSet>& supports,
                    const IndexSet& rs);

protected:
    IndexSet compute0(VectorArray& matrix,
                    VectorArray& vs,
                    std::vector<IndexSet>& supports,
                    const IndexSet& rs);
    IndexSet compute1(VectorArray& matrix,
                    VectorArray& vs,
                    std::vector<IndexSet>& supports,
                    const IndexSet& rs);
    IndexSet compute2(VectorArray& matrix,
                    VectorArray& vs,
                    std::vector<IndexSet>& supports,
                    const IndexSet& rs);
    IndexSet compute3(VectorArray& matrix,
                    VectorArray& vs,
                    std::vector<IndexSet>& supports,
                    const IndexSet& rs);
    IndexSet compute4(VectorArray& matrix,
                    VectorArray& vs,
                    std::vector<IndexSet>& supports,
                    const IndexSet& rs);
};

} // namespace _4ti2_

// Definitions of template class functions.
#include "groebner/RaySupportAlgorithm.tpp"

#endif
