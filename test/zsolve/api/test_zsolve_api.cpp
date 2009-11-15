/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
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

#include <iostream>
#include <fstream>

#include "4ti2/4ti2.h"
#include "4ti2/4ti2xx.h"

int
main()
{
    // Input data.
    const int m = 4;
    const int n = 3;
    int64_t mat[m][n] = {
                { 2,  3, -6 },
                { 2, -1, -4 },
                { 1,  2, -11 },
                { 1, -1,  1 }
            };
    int64_t rel[m] = { _4ti2_LB, _4ti2_UB, _4ti2_UB, _4ti2_LB };
    int64_t sign[n] = { _4ti2_LB, _4ti2_LB, _4ti2_LB };

///     // Output data.
///     const int k = 4;
///     int64_t zhom[k][n] = {
///                 { 3, 4, 1 },
///                 { 3, 8, 5 },
///                 { 9, 2, 4 },
///                 {19,18, 5 }
///             };

    _4ti2_state* zsolve_api = _4ti2_zsolve_create_state(_4ti2_PREC_INT_64);
    const int argc = 2;
    char*argv[2] = { "zsolve", "-q" };
    _4ti2_state_set_options(zsolve_api, argc, argv);

    _4ti2_matrix* cons_matrix;
    _4ti2_state_create_matrix(zsolve_api, m, n, "mat", &cons_matrix);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            _4ti2_matrix_set_entry_int64_t(cons_matrix, i, j,  mat[i][j]);
        }
    }
    //_4ti2_matrix_write_to_stdout(cons_matrix);

    _4ti2_matrix* rel_matrix;
    _4ti2_state_create_matrix(zsolve_api, 1, m, "rel", &rel_matrix);
    for (int i = 0; i < m; ++i) {
        _4ti2_matrix_set_entry_int64_t(rel_matrix, 0, i, rel[i]);
    }
    //_4ti2_matrix_write_to_stdout(rel_matrix);

    _4ti2_matrix* sign_matrix;
    _4ti2_state_create_matrix(zsolve_api, 1, n, "sign", &sign_matrix);
    for (int i = 0; i < n; ++i) {
        _4ti2_matrix_set_entry_int64_t(sign_matrix, 0, i, sign[i]);
    }
    //_4ti2_matrix_write_to_stdout(sign_matrix);

    _4ti2_state_compute(zsolve_api);

    _4ti2_matrix* zhom_matrix;
    _4ti2_state_get_matrix(zsolve_api, "zhom", &zhom_matrix);
    //_4ti2_matrix_write_to_stdout(zhom_matrix);

    // Check the output
///     if (_4ti2_matrix_get_num_rows(zhom_matrix) != k) { return 1; }
///     if (_4ti2_matrix_get_num_cols(zhom_matrix) != n) { return 1; }
///     for (int i = 0; i < k; ++i) {
///         for (int j = 0; j < n; ++j) {
///             int64_t value;
///             _4ti2_matrix_get_entry_int64_t(zhom_matrix, i, j, &value);
///             if (value != zhom[i][j]) { return 1; }
///         }
///     }

    _4ti2_state_delete(zsolve_api);
    return 0;
}   

