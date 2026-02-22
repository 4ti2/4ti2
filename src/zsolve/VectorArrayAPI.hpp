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

#ifndef _4ti2_zsolve__VectorArrayAPI_
#define _4ti2_zsolve__VectorArrayAPI_

#include "4ti2/4ti2xx.h"
#include "zsolve/VectorArray.hpp"
#include "zsolve/Exception.h"
#ifdef _4ti2_HAVE_GMP
#include "4ti2/gmp_integer.h"
#endif
#include <fstream>
#include <cstdint>
#include <cstdlib>
#include <climits>

namespace _4ti2_zsolve_ {

template <class T>
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

#ifdef _4ti2_HAVE_GMP
    virtual void set_entry_mpz_ptr(int r, int c, mpz_srcptr value);
    virtual void get_entry_mpz_ptr(int r, int c, mpz_ptr value) const;
#endif

protected:
    //template <class T1, class T2>
    //static void convert(const T1&, T2&);

public:
    VectorArray<T> data;
};

template <class T1, class T2>
inline
void
convert(const T1& v1, T2& v2)
{
    v2 = v1;
}

template <>
inline
void
convert(const int64_t& v1, int32_t& v2)
{
    if (v1 < INT32_MIN || v1 > INT32_MAX) {
        uint64_t magnitude;
        if (v1 < 0) {
            magnitude = static_cast<uint64_t>(-(v1 + 1)) + 1;
        }
        else {
            magnitude = static_cast<uint64_t>(v1);
        }

        int required_precision = 0;
        do {
            ++required_precision;
            magnitude >>= 1;
        } while (magnitude != 0);

        throw PrecisionException(required_precision);
    }
    v2 = v1;
}

#ifdef _4ti2_HAVE_GMP

template <>
inline
void
convert(const mpz_srcptr& v1, int64_t& v2)
{
    // TODO: Better precision exception information.
    if (!mpz_fits_int64_p(v1)) {
        throw PrecisionException(0);
    }
    v2 = mpz_get_int64(v1);
}

template <>
inline
void
convert(const mpz_srcptr& v1, int32_t& v2)
{
    if (!mpz_fits_sint_p(v1)) {
        throw PrecisionException(0);
    }
    v2 = static_cast<int32_t>(mpz_get_si(v1));
}

template <>
inline
void
convert(const int64_t& v1, mpz_ptr& v2)
{
    mpz_set_int64(v2, v1);
}

template <>
inline
void
convert(const int32_t& v1, mpz_ptr& v2)
{
    mpz_set_si(v2, static_cast<long>(v1));
}

template <>
inline
void
convert(const mpz_srcptr& v1, _4ti2_GMP_INTEGER::Integer& v2)
{
    v2.set_mpz(v1);
}

template <>
inline
void
convert(const _4ti2_GMP_INTEGER::Integer& v1, mpz_ptr& v2)
{
    mpz_set(v2, v1.get_mpz_t());
}

template <>
inline
void
convert(const _4ti2_GMP_INTEGER::Integer& v1, int64_t& v2)
{
    if (!mpz_fits_int64_p(v1.get_mpz_t())) {
        throw PrecisionException(0);
    }
    v2 = mpz_get_int64(v1.get_mpz_t());
}

template <>
inline
void
convert(const _4ti2_GMP_INTEGER::Integer& v1, int32_t& v2)
{
    if (!mpz_fits_sint_p(v1.get_mpz_t())) {
        throw PrecisionException(0);
    }
    v2 = static_cast<int32_t>(mpz_get_si(v1.get_mpz_t()));
}

#endif


template <class T>
VectorArrayAPI<T>::VectorArrayAPI(int num_rows, int num_cols)
    : data(num_rows, num_cols, 0)
{
}

template <class T>
VectorArrayAPI<T>::~VectorArrayAPI()
{
}

template <class T>
int
VectorArrayAPI<T>::get_num_rows() const
{
    return data.height();
}

template <class T>
int
VectorArrayAPI<T>::get_num_cols() const
{
    return data.width();
}

template <class T>
void
VectorArrayAPI<T>::write(const char* filename) const
{
    std::ofstream out(filename);
    if (!out.good()) { throw IOException(std::string("Could not open file ") + filename); }
    data.write(out);
}

template <class T>
void
VectorArrayAPI<T>::write(std::ostream& out) const
{
    data.write(out);
}

template <class T>
void
VectorArrayAPI<T>::read(std::istream& in)
{
    data.read(in, false);
}

template <class T>
void
VectorArrayAPI<T>::set_entry_int32_t(int r, int c, const int32_t& value)
{
    convert(value, data[r][c]);
}

template <class T>
void
VectorArrayAPI<T>::get_entry_int32_t(int r, int c, int32_t& value) const
{
    convert(data[r][c], value);
}

template <class T>
void
VectorArrayAPI<T>::set_entry_int64_t(int r, int c, const int64_t& value)
{
    convert(value, data[r][c]);
}

template <class T>
void
VectorArrayAPI<T>::get_entry_int64_t(int r, int c, int64_t& value) const
{
    convert(data[r][c], value);
}

#ifdef _4ti2_HAVE_GMP
template <class T>
void
VectorArrayAPI<T>::set_entry_mpz_ptr(int r, int c, mpz_srcptr value)
{
    convert(value, data[r][c]);
}

template <class T>
void
VectorArrayAPI<T>::get_entry_mpz_ptr(int r, int c, mpz_ptr value) const
{
    convert(data[r][c], value);
}
#endif


} // namspace _4ti2_

#endif

