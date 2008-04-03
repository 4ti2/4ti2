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

#ifdef _4ti2_GMP_
#include <gmpxx.h>
#endif

#include <sstream>

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

inline int64_t abs (int64_t a)
{
    if (a >= 0)
        return a;
    else
        return -a;
}

inline int32_t gcd (int32_t a, int32_t b)
{
    int32_t tmp;
    while (b != 0)
    {
        tmp = a % b;
        a = b;
        b = tmp;
    }
    return a < 0 ? -a : a;
}

inline int64_t gcd (int64_t a, int64_t b)
{
    int64_t tmp;
    while (b != 0)
    {
        tmp = a % b;
        a = b;
        b = tmp;
    }
    return a < 0 ? -a : a;
}

#ifdef _4ti2_GMP_
inline mpz_class gcd (const mpz_class& a, const mpz_class& b)
{
    mpz_class result;

    mpz_gcd (result.get_mpz_t (), a.get_mpz_t (), b.get_mpz_t ());

    return result;
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

#ifdef _4ti2_GMP_

template <typename T> void parse_integer (const std::string& string, mpz_class& result)
{
    result = string.c_str();
}

#endif

// integer space

template <typename T> int integer_space (const T& number)
{
    std::ostringstream oss;
    oss << number;
    return oss.str ().size ();
}

} // namespace _4ti2_zsolve_

#endif
