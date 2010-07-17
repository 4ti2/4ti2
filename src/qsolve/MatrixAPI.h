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

#ifndef _4ti2_qsolve__MatrixAPI_
#define _4ti2_qsolve__MatrixAPI_

#include "4ti2/4ti2xx.h"
#include "qsolve/VectorArray.h"
#include <iostream>

namespace _4ti2_ {

// Wrapper for a VectorArray object.
template <class MatrixT>
class MatrixAPI : public _4ti2_matrix {
public:
    MatrixAPI();
    MatrixAPI(Size m, Size n);
    virtual void resize(Size m, Size n);
    virtual ~MatrixAPI();

    virtual int get_num_rows() const;
    virtual int get_num_cols() const;

    virtual void write(std::ostream& out) const; 
    virtual void read(std::istream& in);
    virtual void assign(const _4ti2_matrix& m);
    virtual void swap(_4ti2_matrix& m);

#define _MatrixAPI_declare_(TYPE) \
    virtual void set_entry(int r, int c, const TYPE& v);  \
    virtual void get_entry(int r, int c, TYPE& v) const;

    _MatrixAPI_declare_(int32_t)
    _MatrixAPI_declare_(int64_t)
#ifdef _4ti2_GMP_
    _MatrixAPI_declare_(mpz_class)
#endif

    virtual void set_entry(int r, int c, std::istream& in);
    virtual void get_entry(int r, int c, std::ostream& out) const;

public:
    MatrixT data;
};

} // namspace _4ti2_

#include "qsolve/MatrixAPI.hpp"

#endif

