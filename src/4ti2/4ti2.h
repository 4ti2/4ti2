/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright(C) 2006 4ti2 team.
Main author(s): Matthias Walter.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 
*/

#ifndef _4ti2_H
#define _4ti2_H

#include <inttypes.h>

#ifdef _4ti2_GMP_
#include <gmp.h>
#endif

#ifdef __cplusplus
extern "C" 
{
#endif

// Enum representing the possible arithmetic precision settings available.
typedef enum { _4ti2_PREC_INT_32 = 32, _4ti2_PREC_INT_64 = 64, _4ti2_PREC_INT_ARB = 0 } _4ti2_precision;

// Enum representing values describing the constraints on a variable or row of the constraint matrix.
typedef enum { _4ti2_FR = 0, _4ti2_LB = 1, _4ti2_UB = -1, _4ti2_DB = 2, _4ti2_FX = 3 } _4ti2_constraint;

// Enum representing the exit status of an API call to 4ti2.
typedef enum { _4ti2_OK = 0, _4ti2_ERROR = 1 } _4ti2_status;

// 4ti2 data structures.
typedef struct _4ti2_state _4ti2_state;
typedef struct _4ti2_matrix _4ti2_matrix;

// Create a QSolve 4ti2 state object.
_4ti2_status _4ti2_qsolve_create_state(_4ti2_precision prec, _4ti2_state** state);

// Create a QSolve 4ti2 rays object.
_4ti2_status _4ti2_rays_create_state(_4ti2_precision prec, _4ti2_state** state);

// Create a QSolve 4ti2 circuits object.
_4ti2_status _4ti2_circuits_create_state(_4ti2_precision prec, _4ti2_state** state);

// Read in options for the 4ti2 state object.
// These options are exactly the same as the command line options without the project filename at the end.
// Note that argv[0] is ignored!
_4ti2_status _4ti2_state_set_options(_4ti2_state* state, int argc, char** argv);

// Deletes a 4ti2 state object.
void _4ti2_state_delete(_4ti2_state* state);

// Runs the main algorithm of the 4ti2 state object.
_4ti2_status _4ti2_state_compute(_4ti2_state* state);

// Create a 4ti2 matrix. Previous matrix is deleted if it exists.  Pointer is 0 if "name" is not valid.
_4ti2_status _4ti2_state_create_matrix(_4ti2_state* state, int num_rows, int num_cols, const char* name, _4ti2_matrix** matrix);

// Read a 4ti2 matrix from a file.  Previous matrix is deleted if it exists. Returns 0 if "name" is not valid.
_4ti2_status _4ti2_state_read_matrix(_4ti2_state* state, const char* filename, const char* name, _4ti2_matrix** matrix);

// Get a 4ti2 matrix.  Returns 0 if "name" is not valid or if matrix has not been created.
_4ti2_status _4ti2_state_get_matrix(_4ti2_state* state, const char* name, _4ti2_matrix** matrix);

// Returns the number of rows of the matrix.
int _4ti2_matrix_get_num_rows(const _4ti2_matrix*  matrix);

// Returns the number of columns of the matrix.
int _4ti2_matrix_get_num_cols(const _4ti2_matrix*  matrix);

// Write the 4ti2 matrix to stdout.
void _4ti2_matrix_write_to_stdout(const _4ti2_matrix*  matrix);

// Write the 4ti2 matrix to sterr.
void _4ti2_matrix_write_to_stderr(const _4ti2_matrix*  matrix);

// Write the 4ti2 matrix to the file called "filename".
void _4ti2_matrix_write_to_file(const _4ti2_matrix*  matrix, const char* filename);

// Operations on the matrix.
_4ti2_status _4ti2_matrix_set_entry_int32_t(_4ti2_matrix*  matrix, int r, int c, int32_t value);

_4ti2_status _4ti2_matrix_get_entry_int32_t(const _4ti2_matrix*  matrix, int r, int c, int32_t* value);

_4ti2_status _4ti2_matrix_set_entry_int64_t(_4ti2_matrix*  matrix, int r, int c, int64_t value);

_4ti2_status _4ti2_matrix_get_entry_int64_t(const _4ti2_matrix*  matrix, int r, int c, int64_t* value);

#ifdef _4ti2_GMP_
_4ti2_status _4ti2_matrix_set_entry_mpz_ptr(_4ti2_matrix*  matrix, int r, int c, mpz_ptr value);

_4ti2_status _4ti2_matrix_get_entry_mpz_ptr(const _4ti2_matrix*  matrix, int r, int c, mpz_ptr value);
#endif

#ifdef __cplusplus
} // extern "C"
#endif

#endif
