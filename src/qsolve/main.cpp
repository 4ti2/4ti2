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
#include <string>
#include <cstring>

#include "qsolve/DataType.h"
#include "qsolve/QSolveAPI.h"
#include "qsolve/RaysAPI.h"
#include "qsolve/CircuitsAPI.h"

//#define DEBUG_4ti2(X) X
#include "qsolve/Debug.h"
#include "qsolve/Globals.h"

using namespace _4ti2_;

int
main(int argc, char **argv)
{
    if (argc == 1) {
        std::cerr << "Usage: 4ti2 <function> ...\n";
        exit(1);
    }

    _4ti2_state* state = 0;
    std::string exec_name(argv[1]);
    if (exec_name == "qsolve") {
        state = new QSolveAPI();
    } else if (exec_name == "rays") {
        state = new RaysAPI();
    } else if (exec_name == "circuits") {
        state = new CircuitsAPI();
    } else {
        std::cerr << "ERROR: Unrecognized 4ti2 function name: " << argv[1] << "\n";
        exit(1);
    }

    state->set_options(argc-1, argv+1);
    // TODO
    state->read(argv[argc-1]);
    state->compute();
    state->write(argv[argc-1]);

    delete state;
    return 0;
}


