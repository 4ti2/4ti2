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

#include "4ti2/4ti2.h"

int
main()
{
    // Input data, like test/circuits/ppi3
    const int m = 1;
    const int n = 3;
    _4ti2_int64_t mat[m][n] = {
                { 1,  2, 3 }
            };

    // Output data.
    const int k = 3;
    _4ti2_int64_t cir[k][n] = {
                { 0, 3, -2 },
                { 2, -1, 0 },
                { 3, 0, -1 }
            };

    _4ti2_state* circuits_api = _4ti2_circuits_create_state(_4ti2_PREC_INT_64);
    const int argc = 2;
    char*argv[2] = { (char*) "circuits", (char*) "-q" };
    _4ti2_state_set_options(circuits_api, argc, argv);

    _4ti2_matrix* cons_matrix;
    _4ti2_state_create_matrix(circuits_api, m, n, "mat", &cons_matrix);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            _4ti2_matrix_set_entry_int64_t(cons_matrix, i, j,  mat[i][j]);
        }
    }
    //_4ti2_matrix_write_to_stdout(cons_matrix);

    _4ti2_state_compute(circuits_api);

    _4ti2_matrix* cir_matrix;
    _4ti2_state_get_matrix(circuits_api, "cir", &cir_matrix);
    //_4ti2_matrix_write_to_stdout(cir_matrix);

    // Check the output
    if (_4ti2_matrix_get_num_rows(cir_matrix) != k) { return 1; }
    if (_4ti2_matrix_get_num_cols(cir_matrix) != n) { return 1; }
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < n; ++j) {
            _4ti2_int64_t value;
            _4ti2_matrix_get_entry_int64_t(cir_matrix, i, j, &value);
            if (value != cir[i][j]) { return 1; }
        }
    }

    _4ti2_state_delete(circuits_api);
    return 0;
}   

