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

#include <fstream>
#include "groebner/VectorArrayAPI.h"
#include "groebner/VectorArrayStream.h"

using namespace _4ti2_;

VectorArrayAPI::VectorArrayAPI(int num_rows, int num_cols)
    : data(num_rows, num_cols, 0)
{
}

VectorArrayAPI::~VectorArrayAPI()
{
}

int
VectorArrayAPI::get_num_rows() const
{
    return data.get_number();
}

int
VectorArrayAPI::get_num_cols() const
{
    return data.get_size();
}

void
VectorArrayAPI::write(const char* filename) const
{
    std::ofstream out(filename);
    write(out);
}

void
VectorArrayAPI::write(std::ostream& out) const
{
    output(out, data);
}

void
VectorArrayAPI::read(std::istream& in)
{
    in >> data;
}

void
VectorArrayAPI::set_entry_int32_t(int r, int c, const int32_t& value)
{
    convert(value, data[r][c]);
}

void
VectorArrayAPI::get_entry_int32_t(int r, int c, int32_t& value) const
{
    convert(data[r][c], value);
}

void
VectorArrayAPI::set_entry_int64_t(int r, int c, const int64_t& value)
{
    convert(value, data[r][c]);
}

void
VectorArrayAPI::get_entry_int64_t(int r, int c, int64_t& value) const
{
    convert(data[r][c], value);
}

#ifdef _4ti2_GMP_
void
VectorArrayAPI::set_entry_mpz_class(int r, int c, const mpz_class& value)
{
    convert(value, data[r][c]);
}

void
VectorArrayAPI::get_entry_mpz_class(int r, int c, mpz_class& value) const
{
    convert(data[r][c], value);
}
#endif
