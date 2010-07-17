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

#ifndef __4ti2_qsolve__TypeConversion_
#define __4ti2_qsolve__TypeConversion_

namespace _4ti2_
{

// TODO: Probably should not have the default case.
template <class T1, class T2>
inline
void
type_conversion(const T1& v1, T2& v2)
{
    v2 = v1;
}

template <>
inline
void
type_conversion(const int64_t& v1, int32_t& v2)
{
    if (v1 < INT32_MIN || v1 > INT32_MAX) {
        std::cerr << "ERROR: number " << v1 << " out of range.\n";
        std::cerr << "ERROR: range is (" << INT32_MIN << "," << INT32_MAX << ").\n";
        exit(1);    
    }
    v2 = v1;
}

#ifdef _4ti2_GMP_

template <>
inline
void
type_conversion(const int32_t& v1, mpz_class& v2)
{
    v2 = v1;
}

template <>
inline
void
type_conversion(const int64_t& v1, mpz_class& v2)
{
    v2 = v1;
}

template <>
inline
void
type_conversion(const mpz_class& v1, int& v2)
{
    if (!v1.fits_sint_p()) {
        std::cerr << "ERROR: number " << v1 << " out of range.\n";
        std::cerr << "ERROR: range is (" << INT_MIN << "," << INT_MAX << ").\n";
        exit(1);    
    }
    v2 = v1.get_si();
}

template <>
inline
void
type_conversion(const mpz_class& v1, long int& v2)
{
    if (!v1.fits_slong_p()) {
        std::cerr << "ERROR: number " << v1 << " out of range.\n";
        std::cerr << "ERROR: range is (" << LONG_MIN << "," << LONG_MAX << ").\n";
        exit(1);    
    }
    v2 = v1.get_si();
}

#endif

}

#endif
