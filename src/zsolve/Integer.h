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

#ifdef HAVE_GMP
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

#ifdef HAVE_GMP
inline mpz_class gcd (const mpz_class& a, const mpz_class& b)
{
    mpz_class result;

    mpz_gcd (result.get_mpz_t (), a.get_mpz_t (), b.get_mpz_t ());

    return result;
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

#ifdef HAVE_GMP
inline int calcPrecision (const mpz_class& n)
{
    return mpz_sizeinbase (n.get_mpz_t (), 2);
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

#ifdef HAVE_GMP
inline int maxPrecision (const mpz_class& n)
{
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

#ifdef HAVE_GMP

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
