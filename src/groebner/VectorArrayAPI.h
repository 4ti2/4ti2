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

#ifndef _4ti2_groebner__VectorArrayAPI_
#define _4ti2_groebner__VectorArrayAPI_

#include "groebner/4ti2API.h"
#include "groebner/VectorArray.h"

namespace _4ti2_ {

class VectorArrayAPI : public _4ti2_matrix {
public:
    VectorArrayAPI(int num_rows, int num_cols);
    virtual ~VectorArrayAPI();

    virtual int get_num_rows() const;
    virtual int get_num_cols() const;

    virtual void write(const char* filename) const;
    virtual void write(std::ostream& out) const; 
    virtual void read(std::istream& in);

    virtual void set_entry_int32_t(int r, int c, const int32_t& value); 
    virtual void get_entry_int32_t(int r, int c, int32_t& value) const;
    virtual void set_entry_int64_t(int r, int c, const int64_t& value);
    virtual void get_entry_int64_t(int r, int c, int64_t& value) const;

#ifdef _4ti2_GMP_
    virtual void set_entry_mpz_class(int r, int c, const mpz_class& value);
    virtual void get_entry_mpz_class(int r, int c, mpz_class& value) const;
#endif

protected:
    template <class T1, class T2>
    static void convert(const T1&, T2&);

public:
    VectorArray data;
};

template <class T1, class T2>
inline
void
VectorArrayAPI::convert(const T1& v1, T2& v2)
{
    v2 = v1;
}

template <>
inline
void
VectorArrayAPI::convert(const int64_t& v1, int32_t& v2)
{
    if (v1 < INT32_MIN || v2 > INT32_MAX) {
        std::cerr << "ERROR: number " << v1 << " out of range.\n";
        std::cerr << "ERROR: range is (" << INT32_MIN << "," << INT32_MAX << ").\n";
        exit(1);    
    }
    v2 = v1;
}

template <>
inline
void
VectorArrayAPI::convert(const int32_t& v1, int64_t& v2)
{
    v2 = v1;
}

#ifdef _4ti2_GMP_

template <>
inline
void
VectorArrayAPI::convert(const int32_t& v1, mpz_class& v2)
{
    v2 = v1;
}

template <>
inline
void
VectorArrayAPI::convert(const int64_t& v1, mpz_class& v2)
{
    v2 = v1;
}

template <>
inline
void
VectorArrayAPI::convert(const mpz_class& v1, int& v2)
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
VectorArrayAPI::convert(const mpz_class& v1, long int& v2)
{
    if (!v1.fits_slong_p()) {
        std::cerr << "ERROR: number " << v1 << " out of range.\n";
        std::cerr << "ERROR: range is (" << LONG_MIN << "," << LONG_MAX << ").\n";
        exit(1);    
    }
    v2 = v1.get_si();
}

#endif


} // namspace _4ti2_

#endif

