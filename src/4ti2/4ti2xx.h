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

#ifndef _4ti2API_
#define _4ti2API_

#include <iostream>
#include <fstream>

#include "4ti2/4ti2_config.h"

#ifdef _4ti2_HAVE_GMP
#include <gmp.h>
#include <gmpxx.h>
#endif

class _4ti2_matrix {
public:
    _4ti2_matrix() {}
    virtual ~_4ti2_matrix() {}
    virtual void resize(int m, int n) = 0;

    virtual int get_num_rows() const = 0;
    virtual int get_num_cols() const = 0;

    virtual void write(std::ostream& out) const = 0; 
    virtual void write(const char* filename) const;
    virtual void read(std::istream& in) = 0; 
    virtual void assign(const _4ti2_matrix& m) = 0;
    virtual void swap(_4ti2_matrix& m) = 0;

// TODO
#define _4ti2_matrix_declare(TYPE) \
    virtual void set_entry(int r, int c, const TYPE& value) = 0;  \
    virtual void set_row(int r, const TYPE* row) {}; \
    virtual void set_col(int c, const TYPE* col) {}; \
    virtual void set_matrix_row_major(const TYPE* mat) {}; \
    virtual void set_matrix_col_major(const TYPE* mat) {}; \
    virtual void get_entry(int r, int c, TYPE& value) const = 0; \
    virtual void get_row(int r, TYPE* row) {}; \
    virtual void get_col(int c, TYPE* col) {}; \
    virtual void get_matrix_row_major(TYPE* mat) {}; \
    virtual void get_matrix_col_major(TYPE* mat) {};

    _4ti2_matrix_declare(int32_t)
    _4ti2_matrix_declare(int64_t)
#ifdef _4ti2_HAVE_GMP
    _4ti2_matrix_declare(mpz_class)
#endif
};

class _4ti2_state {
public:
    _4ti2_state() {}
    virtual ~_4ti2_state() {}

    virtual void compute() = 0;

    virtual void set_options(int argc, char** argv) = 0; 

    virtual void read(const char* project) = 0;
    virtual void write(const char* project) = 0;

    virtual _4ti2_matrix* create_matrix(int num_rows, int num_cols, const char* name) = 0;
    virtual _4ti2_matrix* create_matrix(const char* filename, const char* name) = 0;
    virtual _4ti2_matrix* create_matrix(std::istream& in, const char* name) = 0;

    virtual _4ti2_matrix* get_matrix(const char* name) = 0;
};

inline
void
_4ti2_matrix::write(const char* filename) const
{
    std::ofstream file(filename);
    write(file);
}

#endif
