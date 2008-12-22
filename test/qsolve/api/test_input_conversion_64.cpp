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

#include <iostream>
#include <fstream>
#include <sstream>

#include "4ti2/4ti2.h"
#include "4ti2/4ti2xx.h"

using namespace std;

#define APITYPE int64_t

static void
usage()
{
  cerr << "usage: ..." << endl;
  exit(1);
}

#define CHECK_STATUS(expr)			\
if ((expr) != _4ti2_OK ) {			\
  cerr << "Error on " << #expr << endl;		\
  exit (1);					\
}

static void check_matrix(_4ti2_state* qsolve_api, const char *name, APITYPE x)
{
  _4ti2_matrix* matrix;
  APITYPE y;
  CHECK_STATUS(_4ti2_state_create_matrix(qsolve_api, 1, 1, name, &matrix));
  CHECK_STATUS(_4ti2_matrix_set_entry_int64_t(matrix, 0, 0, x));
  CHECK_STATUS(_4ti2_matrix_get_entry_int64_t(matrix, 0, 0, &y));
  if (x != y) {
    cerr << "Data conversion failed: " << name << ": " << x << " != " << y << endl;
    exit(1);
  }
}

// usage: test_input_conversion_APIPRECISION INTERNAL-PRECISION MATRIX-NAME DATA
int
main(int argc, char **argv)
{
    // Input data.
    if (argc != 4) usage();

    _4ti2_precision prec;
    switch (atoi(argv[1])) {
    case 0:
      prec = _4ti2_PREC_INT_ARB;
      break;
    case 32:
      prec = _4ti2_PREC_INT_32;
      break;
    case 64:
      prec = _4ti2_PREC_INT_64;
      break;
    default:
      usage();
    }
    
    _4ti2_state* qsolve_api = _4ti2_qsolve_create_state(prec);
    const int qsolve_argc = 2;
    char *qsolve_argv[2] = { "qsolve", "-q" };
    _4ti2_state_set_options(qsolve_api, qsolve_argc, qsolve_argv);

    APITYPE x;
    stringstream datastream(argv[3]);
    datastream >> x;
    if (datastream.bad()) {
      cerr << "Error reading data" << endl;
      exit(1);
    }
    check_matrix(qsolve_api, argv[2], x);
    
    _4ti2_state_delete(qsolve_api);
    return 0;
}   

