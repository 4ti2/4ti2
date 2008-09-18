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
#ifdef _4ti2_GMP_
#include <gmp.h>
#include <gmpxx.h>
#endif

struct _4ti2_matrix {
public:
    _4ti2_matrix() {}
    virtual ~_4ti2_matrix() {}

    virtual int get_num_rows() const = 0;
    virtual int get_num_cols() const = 0;

    virtual void write(const char* filename) const = 0;
    virtual void write(std::ostream& out) const = 0; 
    virtual void read(std::istream& in) = 0; 

    virtual void set_entry_int32_t(int r, int c, const int32_t& value) = 0; 
    virtual void get_entry_int32_t(int r, int c, int32_t& value) const = 0;
    virtual void set_entry_int64_t(int r, int c, const int64_t& value) = 0;
    virtual void get_entry_int64_t(int r, int c, int64_t& value) const = 0;

#ifdef _4ti2_GMP_
    virtual void set_entry_mpz_class(int r, int c, const mpz_class& value) = 0;
    virtual void get_entry_mpz_class(int r, int c, mpz_class& value) const = 0;
#endif
};

struct _4ti2_state {
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

#endif
