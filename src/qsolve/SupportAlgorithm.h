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

#ifndef _4ti2_qsolve__SupportAlgorithm_
#define _4ti2_qsolve__SupportAlgorithm_

#include "qsolve/QSolveAlgorithm.h"
#include "qsolve/Vector.h"
#include "qsolve/VectorArray.h"
#include "qsolve/QSolveConsOrder.h"
#include "qsolve/RayState.h"
#include "qsolve/SubAlgorithm.h"

namespace _4ti2_
{

template <class IndexSet>
class SupportSubAlgorithmBase : public SubAlgorithm
{
public:
    SupportSubAlgorithmBase(RayStateAPI<IndexSet>& state);
    virtual ~SupportSubAlgorithmBase();

    void compute_rays(Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end);
    void compute_cirs(Index r1_start, Index r1_end, Index r2_start, Index r2_index, Index r2_end);
    void transfer();

protected:
    RayStateAPI<IndexSet>& state;
    RaySubStateAPI<IndexSet>& helper;
};

template <class IndexSet>
class SupportRayAlgorithm : public SupportSubAlgorithmBase<IndexSet>
{
public:
    SupportRayAlgorithm(RayStateAPI<IndexSet>& state, IndexRanges& indices);
    SupportRayAlgorithm* clone();

    virtual void compute();

protected:
    IndexRanges& indices;
};

template <class IndexSet>
class SupportCirAlgorithm : public SupportSubAlgorithmBase<IndexSet>
{
public:
    SupportCirAlgorithm(RayStateAPI<IndexSet>& state, IndexRanges& indices);
    SupportCirAlgorithm* clone();

    virtual void compute();

protected:
    IndexRanges& indices;
};


} // namespace _4ti2_

// Definitions of template class functions.
#include "qsolve/SupportAlgorithm.hpp"

#endif
