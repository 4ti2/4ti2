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

#ifndef _4ti2_qsolve__ConeAPI_
#define _4ti2_qsolve__ConeAPI_

#include "4ti2/4ti2xx.h"
#include "qsolve/Cone.h"
#include "qsolve/VectorArray.h"
#include "qsolve/VectorArrayStream.h"
#include "qsolve/TypeConversion.h"

namespace _4ti2_ {

// Wrapper for a VectorArray object.
template <class T>
class ConeAPI : public _4ti2_matrix {
public:
    ConeAPI();
    ConeAPI(Size m, Size n);
    void init(Size m, Size n);
    virtual ~ConeAPI();

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

public:
    ConeT<T> cone;
};

template <class T>
ConeAPI<T>::ConeAPI()
    : cone()
{
}

template <class T>
ConeAPI<T>::ConeAPI(Size m, Size n)
    : cone(m, n)
{
}

template <class T>
void
ConeAPI<T>::init(Size m, Size n)
{
    cone.init(m, n);
}

template <class T>
ConeAPI<T>::~ConeAPI()
{
}

template <class T>
int
ConeAPI<T>::get_num_rows() const
{
    return cone.num_cons();
}

template <class T>
int
ConeAPI<T>::get_num_cols() const
{
    return cone.num_vars();
}

template <class T>
void
ConeAPI<T>::write(const char* filename) const
{
    std::ofstream out(filename);
    write(out);
}

template <class T>
void
ConeAPI<T>::write(std::ostream& out) const
{
    output(out, cone.get_matrix());
}

template <class T>
void
ConeAPI<T>::read(std::istream& in)
{
    in >> cone.get_matrix();
}

template <class T>
void
ConeAPI<T>::set_entry_int32_t(int r, int c, const int32_t& value)
{
    type_conversion(value, cone.get_matrix()[r][c]);
}

template <class T>
void
ConeAPI<T>::get_entry_int32_t(int r, int c, int32_t& value) const
{
    type_conversion(cone.get_matrix()[r][c], value);
}

template <class T>
void
ConeAPI<T>::set_entry_int64_t(int r, int c, const int64_t& value)
{
    type_conversion(value, cone.get_matrix()[r][c]);
}

template <class T>
void
ConeAPI<T>::get_entry_int64_t(int r, int c, int64_t& value) const
{
    type_conversion(cone.get_matrix()[r][c], value);
}

#ifdef _4ti2_GMP_
template <class T>
void
ConeAPI<T>::set_entry_mpz_class(int r, int c, const mpz_class& value)
{
    type_conversion(value, cone.get_matrix()[r][c]);
}

template <class T>
void
ConeAPI<T>::get_entry_mpz_class(int r, int c, mpz_class& value) const
{
    type_conversion(cone.get_matrix()[r][c], value);
}
#endif

} // namspace _4ti2_

#endif

