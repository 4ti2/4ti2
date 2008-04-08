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

typedef int int32_t;
typedef long long int64_t;

#ifdef __cplusplus
extern "C" 
{
#endif

#define ZSOLVE_RELATION_EQUAL 0
#define ZSOLVE_RELATION_LESSER 1
#define ZSOLVE_RELATION_LESSER_EQUAL 2
#define ZSOLVE_RELATION_GREATER 3
#define ZSOLVE_RELATION_GREATER_EQUAL 4

    typedef void* ZSolveSystem32;

    ZSolveSystem32 zsolve_system_create_32 (int height, int width);
    void zsolve_system_set_matrix_32 (ZSolveSystem32 system, int column, int row, int32_t value);
    void zsolve_system_set_rhs_32 (ZSolveSystem32 system, int row, int32_t value);
    void zsolve_system_set_rel_32 (ZSolveSystem32 system, int row, int relation);
    void zsolve_system_set_free_32 (ZSolveSystem32 system, int column);
    void zsolve_system_set_hilbert_32 (ZSolveSystem32 system, int column);
    void zsolve_system_set_graver_32 (ZSolveSystem32 system, int column);
    void zsolve_system_set_bounded_32 (ZSolveSystem32 system, int column, int32_t lower, int32_t upper);
    void zsolve_system_delete_32 (ZSolveSystem32 system);
    void zsolve_system_print_32 (ZSolveSystem32 system);

    typedef void* ZSolveMatrix32;

    ZSolveMatrix32 zsolve_matrix_create_32 (int height, int width);
    int zsolve_matrix_width_32 (ZSolveMatrix32 matrix);
    int zsolve_matrix_height_32 (ZSolveMatrix32 matrix);
    void zsolve_matrix_set_32 (ZSolveMatrix32 matrix, int column, int row, int32_t value);
    int32_t zsolve_matrix_get_32 (ZSolveMatrix32 matrix, int column, int row);
    void zsolve_matrix_delete_32 (ZSolveMatrix32 matrix);
    void zsolve_matrix_print_32 (ZSolveMatrix32 matrix);

    typedef void* ZSolveState32;

    ZSolveState32 zsolve_state_create_32 (ZSolveSystem32 system);
    void zsolve_state_compute_32 (ZSolveState32 state, ZSolveMatrix32* inhom, ZSolveMatrix32* hom, ZSolveMatrix32* free);
    ZSolveMatrix32 zsolve_state_extract_inhom_32 (ZSolveState32 state);
    ZSolveMatrix32 zsolve_state_extract_hom_32 (ZSolveState32 state);
    ZSolveMatrix32 zsolve_state_extract_free_32 (ZSolveState32 state);
    void zsolve_state_delete_32 (ZSolveState32 state);

#ifdef __cplusplus
}
#endif

#endif
