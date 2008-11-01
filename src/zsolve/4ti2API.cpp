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

#include "4ti2/4ti2.h"
#include "4ti2/4ti2xx.h"
#include "zsolve/ZSolveAPI.hpp"
#include "zsolve/HilbertAPI.hpp"
#include "zsolve/GraverAPI.hpp"

// TODO: Handle errors.

using namespace _4ti2_zsolve_;

extern "C" 
{

_4ti2_state*
_4ti2_zsolve_create_state(_4ti2_precision prec)
{
    switch (prec) {
    case _4ti2_PREC_INT_32:
        return new ZSolveAPI<int32_t>();
    case _4ti2_PREC_INT_64:
        return new ZSolveAPI<int64_t>();
#ifdef _4ti2_GMP_
    case _4ti2_PREC_INT_ARB:
        return new ZSolveAPI<mpz_class>();
#endif
    default: 
        std::cerr << "ERROR: Undefined precision.\n";
        exit(1);
    }
}

_4ti2_state*
_4ti2_hilbert_create_state(_4ti2_precision prec)
{
    switch (prec) {
    case _4ti2_PREC_INT_32:
        return new HilbertAPI<int32_t>();
    case _4ti2_PREC_INT_64:
        return new HilbertAPI<int64_t>();
#ifdef _4ti2_GMP_
    case _4ti2_PREC_INT_ARB:
        return new HilbertAPI<mpz_class>();
#endif
    default: 
        std::cerr << "ERROR: Undefined precision.\n";
        exit(1);
    }
}

_4ti2_state*
_4ti2_graver_create_state(_4ti2_precision prec)
{
    switch (prec) {
    case _4ti2_PREC_INT_32:
        return new GraverAPI<int32_t>();
    case _4ti2_PREC_INT_64:
        return new GraverAPI<int64_t>();
#ifdef _4ti2_GMP_
    case _4ti2_PREC_INT_ARB:
        return new GraverAPI<mpz_class>();
#endif
    default: 
        std::cerr << "ERROR: Undefined precision.\n";
        exit(1);
    }
}

_4ti2_status
_4ti2_state_set_options(_4ti2_state* state, int argc, char** argv)
{
    state->set_options(argc, argv);
    return _4ti2_OK;
}

void
_4ti2_state_delete(_4ti2_state* state)
{
    delete state;
}

_4ti2_status
_4ti2_state_compute(_4ti2_state* state)
{
    state->compute();
    return _4ti2_OK;
}

_4ti2_status
_4ti2_state_create_matrix(
                _4ti2_state* state,
                int num_rows,
                int num_cols,
                const char* name,
                _4ti2_matrix** matrix)
{
    *matrix = state->create_matrix(num_rows, num_cols, name);
    return _4ti2_OK;
}

_4ti2_status
_4ti2_state_read_matrix(
                _4ti2_state* state, 
                const char* filename, 
                const char* name, 
                _4ti2_matrix** matrix)
{
    *matrix = state->create_matrix(filename, name);
    return _4ti2_OK;
}

_4ti2_status
_4ti2_state_get_matrix(
                _4ti2_state* state, 
                const char* name, 
                _4ti2_matrix** matrix)
{
    *matrix = state->get_matrix(name);
    return _4ti2_OK;
}

int
_4ti2_matrix_get_num_rows(const _4ti2_matrix*  matrix)
{
    return matrix->get_num_rows();
}

int
_4ti2_matrix_get_num_cols(const _4ti2_matrix*  matrix)
{
    return matrix->get_num_cols();
}

void
_4ti2_matrix_write_to_stdout(const _4ti2_matrix*  matrix)
{
    matrix->write(std::cout);
}

void
_4ti2_matrix_write_to_stderr(const _4ti2_matrix*  matrix)
{
    matrix->write(std::cerr);
}

void
_4ti2_matrix_write_to_file(const _4ti2_matrix*  matrix, const char* filename)
{
    matrix->write(filename);
}

_4ti2_status
_4ti2_matrix_set_entry_int32_t(_4ti2_matrix*  matrix, int r, int c, int32_t value)
{
    matrix->set_entry_int32_t(r, c, value);
    return _4ti2_OK;
}

_4ti2_status
_4ti2_matrix_get_entry_int32_t(const _4ti2_matrix*  matrix, int r, int c, int32_t* value)
{
    matrix->get_entry_int32_t(r, c, *value);
    return _4ti2_OK;
}

_4ti2_status
_4ti2_matrix_set_entry_int64_t(_4ti2_matrix*  matrix, int r, int c, int64_t value)
{
    matrix->set_entry_int64_t(r, c, value);
    return _4ti2_OK;
}

_4ti2_status
_4ti2_matrix_get_entry_int64_t(const _4ti2_matrix*  matrix, int r, int c, int64_t* value)
{
    matrix->get_entry_int64_t(r, c, *value);
    return _4ti2_OK;
}

#ifdef _4ti2_GMP_
_4ti2_status
_4ti2_matrix_set_entry_mpz_ptr(_4ti2_matrix*  matrix, int r, int c, mpz_ptr value)
{
    mpz_class z(value);
    matrix->set_entry_mpz_class(r, c, z);
    return _4ti2_OK;
}

_4ti2_status
_4ti2_matrix_get_entry_mpz_ptr(const _4ti2_matrix*  matrix, int r, int c, mpz_ptr value)
{
    mpz_class z;
    matrix->get_entry_mpz_class(r, c, z);
    mpz_set(value, z.get_mpz_t());
    return _4ti2_OK;
}
#endif

}

