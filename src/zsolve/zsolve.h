/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Matthias Walter.

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

#ifndef _4ti2_zsolve__zsolve_
#define _4ti2_zsolve__zsolve_

#include <inttypes.h>

#if 0
typedef int int32_t;
typedef long long int64_t;
#endif

#include "4ti2/4ti2_config.h"

#ifdef _4ti2_HAVE_GMP
#include <gmp.h>
#endif

#ifdef __cplusplus
extern "C" 
{
#endif

#define ZSOLVE_RELATION_EQUAL 0
#define ZSOLVE_RELATION_LESSER 1
#define ZSOLVE_RELATION_LESSER_EQUAL 2
#define ZSOLVE_RELATION_GREATER 3
#define ZSOLVE_RELATION_GREATER_EQUAL 4

    typedef void* ZSolveMatrix;
    typedef void* ZSolveState;

    ZSolveState zsolve_state_create (int height, int width, int precision); // creates a height x width system with precision 32 = 32bit, 64 = 64bit or GMP else.
    void zsolve_state_delete (ZSolveState state); // deletes a system with all its associates matrices! Any matrix returned by zsolve_state_matrix is invalid after this call.
    int zsolve_state_compute (ZSolveState state); // computes the solutions for the system. returns 0 for success and 1, if precision is too low.
    ZSolveMatrix zsolve_state_matrix (ZSolveState state, char* name); // returns a matrix. name can be one of ("mat", "rhs", "rel", "sign", "lb", "ub", "zinhom", "zhom", "zfree").

    int zsolve_matrix_width (ZSolveMatrix matrix); // returns the width of a matrix.
    int zsolve_matrix_height (ZSolveMatrix matrix); // returns the height of a matrix.
    int zsolve_matrix_read_only (ZSolveMatrix matrix); // returns, whether the matrix is read-only.
    void zsolve_matrix_delete (ZSolveMatrix matrix); // deletes a matrix. should be called before zsolve_state_delete!
    int zsolve_matrix_set_32 (ZSolveMatrix matrix, int r, int c, int32_t value); // sets the entry at r,c to the 32bit value.
    int zsolve_matrix_set_64 (ZSolveMatrix matrix, int r, int c, int64_t value); // sets the entry at r,c to the 64bit value.
    int32_t zsolve_matrix_get_32 (ZSolveMatrix matrix, int r, int c); // returns the 32bit entry at r,c.
    int64_t zsolve_matrix_get_64 (ZSolveMatrix matrix, int r, int c); // returns the 64bit entry at r,c.
    void zsolve_matrix_print_32 (ZSolveMatrix matrix); // prints the matrix with 32bit values.
    void zsolve_matrix_print_64 (ZSolveMatrix matrix); // prints the matrix with 64bit values.

#ifdef _4ti2_HAVE_GMP
    int zsolve_matrix_set_gmp (ZSolveMatrix matrix, int r, int c, mpz_srcptr value); // sets the entry at r,c to the mpz_t value.
    void zsolve_matrix_get_gmp (ZSolveMatrix matrix, mpz_ptr result, int r, int c); // returns the mpz_t entry at r,c.
    void zsolve_matrix_print_gmp (ZSolveMatrix matrix); // prints the matrix with gmp values.
#endif

#ifdef __cplusplus
}
#endif

#endif
