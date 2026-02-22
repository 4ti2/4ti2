/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Matthias Walter

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

#ifndef __4ti2_zsolve__Integer_
#define __4ti2_zsolve__Integer_

#ifdef __cplusplus
#include <cstdint>
#else
#include <inttypes.h>
#endif

#include "4ti2/4ti2_config.h"

#ifdef _4ti2_HAVE_GMP
#include "4ti2/gmp_integer.h"
#endif

#include <sstream>
#include <stdexcept>

namespace _4ti2_zsolve_
{

// sign 

inline int sgn (int64_t a)
{
    if (a > 0)
        return 1;
    else if (a < 0)
        return -1;
    else
        return 0;
}

inline int sgn (int32_t a)
{
    if (a > 0)
        return 1;
    else if (a < 0)
        return -1;
    else
        return 0;
}

// absolute value

template <typename T>
inline T abs (T a)
{
    if (a >= 0)
        return a;
    else
        return -a;
}

template <typename T>
inline T gcd (T a, T b)
{
    T tmp;
    while (b != 0)
    {
        tmp = a % b;
        a = b;
        b = tmp;
    }
    return a < 0 ? -a : a;
}

#ifdef _4ti2_HAVE_GMP
template <>
inline _4ti2_GMP_INTEGER::Integer gcd (_4ti2_GMP_INTEGER::Integer a, _4ti2_GMP_INTEGER::Integer b)
{
    return _4ti2_GMP_INTEGER::gcd(a, b);
}
#endif

inline int calcPrecision (int32_t n)
{
    if (n < 0)
	n = -n;
    int result = 0;
    while (n != 0)
    {
	result++;
	n /= 2;
    }
    return result;
}

inline int calcPrecision (int64_t n)
{
    if (n < 0)
        n = -n;
    int result = 0;
    while (n != 0)
	{
	    result++;
	    n /= 2;
	}
    return result;
}
 
#ifdef _4ti2_HAVE_GMP
inline int calcPrecision (const _4ti2_GMP_INTEGER::Integer& n)
{
    return _4ti2_GMP_INTEGER::calc_precision(n);
}
#endif

inline int maxPrecision (int32_t n)
{
    return 32;
}

inline int maxPrecision (int64_t n)
{
    return 64;
}

#ifdef _4ti2_HAVE_GMP
inline int maxPrecision (const _4ti2_GMP_INTEGER::Integer& n)
{
    (void)n;
    return -1;
}
#endif



// maximum

template <typename T> T max (T a, T b)
{
    return a > b ? a : b;
}

// minimum

inline int min (int a, int b)
{
    return a < b ? a : b;
}

// integer parsing

template <typename T> void parse_integer (const std::string& string, T& result)
{
    std::istringstream iss (string);
    iss >> result;
}

// integer space

template <typename T> int integer_space (const T& number)
{
    std::ostringstream oss;
    oss << number;
    return oss.str ().size ();
}

} // namespace _4ti2_zsolve_

#endif
