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

#ifndef _4ti2_groebner__DataType_
#define _4ti2_groebner__DataType_
#include <inttypes.h>

#include "4ti2/4ti2_config.h"

#ifdef _4ti2_GMP_

#include <gmp.h>
#include <gmpxx.h>
typedef mpz_class IntegerType;
typedef mpq_class RationalType;

#elif defined(_4ti2_INT64_)

typedef int64_t IntegerType;
typedef double RationalType;

#elif defined(_4ti2_INT32_)

typedef int32_t IntegerType;
typedef float RationalType;

#elif defined(_4ti2_INT16_)

typedef int16_t IntegerType;
typedef float RationalType;

#endif

#endif
