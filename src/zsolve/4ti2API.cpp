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
#include "zsolve/ZSolveAPI.hpp"
#include "zsolve/HilbertAPI.hpp"
#include "zsolve/GraverAPI.hpp"

// TODO: Handle errors.

using namespace _4ti2_zsolve_;

extern "C" 
{

_4ti2_state*
_4ti2_zsolve_create_state(_4ti2_precision prec)
{
    switch (prec) {
    case _4ti2_PREC_INT_32:
        return new ZSolveAPI<int32_t>();
    case _4ti2_PREC_INT_64:
        return new ZSolveAPI<int64_t>();
#ifdef _4ti2_GMP_
    case _4ti2_PREC_INT_ARB:
        return new ZSolveAPI<mpz_class>();
#endif
    default: 
        std::cerr << "ERROR: Undefined precision.\n";
        exit(1);
    }
}

_4ti2_state*
_4ti2_hilbert_create_state(_4ti2_precision prec)
{
    switch (prec) {
    case _4ti2_PREC_INT_32:
        return new HilbertAPI<int32_t>();
    case _4ti2_PREC_INT_64:
        return new HilbertAPI<int64_t>();
#ifdef _4ti2_GMP_
    case _4ti2_PREC_INT_ARB:
        return new HilbertAPI<mpz_class>();
#endif
    default: 
        std::cerr << "ERROR: Undefined precision.\n";
        exit(1);
    }
}

_4ti2_state*
_4ti2_graver_create_state(_4ti2_precision prec)
{
    switch (prec) {
    case _4ti2_PREC_INT_32:
        return new GraverAPI<int32_t>();
    case _4ti2_PREC_INT_64:
        return new GraverAPI<int64_t>();
#ifdef _4ti2_GMP_
    case _4ti2_PREC_INT_ARB:
        return new GraverAPI<mpz_class>();
#endif
    default: 
        std::cerr << "ERROR: Undefined precision.\n";
        exit(1);
    }
}

}

