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
#include <fstream>

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
    virtual void set_entry_mpz_class(int r, int c, const mpz_class& value);
    virtual void get_entry_mpz_class(int r, int c, mpz_class& value) const;
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
    // TODO: Better precision exception information.
    if (v1 < INT32_MIN || v2 > INT32_MAX) {
        throw PrecisionException(32);
    }
    v2 = v1;
}

#ifdef _4ti2_HAVE_GMP

template <>
inline
void
convert(const mpz_class& v1, int& v2)
{
    // TODO: Better precision exception information.
    if (!v1.fits_sint_p()) {
        throw PrecisionException(0);
    }
    v2 = v1.get_si();
}

template <>
inline
void
convert(const mpz_class& v1, long int& v2)
{
    // TODO: Better precision exception information.
    if (!v1.fits_slong_p()) {
        throw PrecisionException(0);
    }
    v2 = v1.get_si();
}

#ifndef _4ti2_HAVE_MPZ_INT64_CONVERSION
template <>
inline
void
convert(const mpz_class& v1, int64_t& v2)
{
  std::cerr << "UNIMPLEMENTED: Need to convert from mpz to int64_t" << std::endl;
  exit(1);
}

template <>
inline
void
convert(const int64_t& v1, mpz_class &v2)
{
  std::cerr << "UNIMPLEMENTED: Need to convert from int64_t to mpz" << std::endl;
  exit(1);
}
#endif

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
VectorArrayAPI<T>::set_entry_mpz_class(int r, int c, const mpz_class& value)
{
    convert(value, data[r][c]);
}

template <class T>
void
VectorArrayAPI<T>::get_entry_mpz_class(int r, int c, mpz_class& value) const
{
    convert(data[r][c], value);
}
#endif


} // namspace _4ti2_

#endif

