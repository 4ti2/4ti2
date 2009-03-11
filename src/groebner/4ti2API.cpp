/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2008 4ti2 team.
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

#include "4ti2/4ti2.h"
#include "4ti2/4ti2xx.h"
#include "groebner/QSolveAPI.h"
#include "groebner/RaysAPI.h"
#include "groebner/CircuitsAPI.h"

// TODO: Handle different precision.
// TODO: Handle errors.

using namespace _4ti2_;

extern "C" 
{

_4ti2_state*
_4ti2_qsolve_create_state(_4ti2_precision prec)
{
    return new _4ti2_::QSolveAPI();
}

_4ti2_state*
_4ti2_rays_create_state(_4ti2_precision prec)
{
    return new _4ti2_::QSolveAPI();
}

_4ti2_state*
_4ti2_circuits_create_state(_4ti2_precision prec)
{
    return new _4ti2_::QSolveAPI();
}

}

