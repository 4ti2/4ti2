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

#include "groebner/qsolve_main.h"
#include "groebner/Vector.h"
#include "groebner/VectorStream.h"
#include "groebner/VectorArray.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/BitSet.h"
#include "groebner/BitSetStream.h"
#include "groebner/QSolveAlgorithm.h"
#include "groebner/LatticeBasis.h"
#include "groebner/CircuitOptions.h"

#include "groebner/QSolveAPI.h"
#include "groebner/RaysAPI.h"
#include "groebner/CircuitsAPI.h"

//#define DEBUG_4ti2(X) X
#include "groebner/Debug.h"
#include "groebner/Globals.h"

using namespace _4ti2_;

int
_4ti2_::qsolve_main(int argc, char **argv)
{
    QSolveAPI* qsolve_api = 0;
    if (!strcmp("qsolve", argv[0])) { qsolve_api = new QSolveAPI; } 
    else if (!strcmp("rays", argv[0])) { qsolve_api = new RaysAPI; }
    else if (!strcmp("circuits", argv[0])) { qsolve_api = new CircuitsAPI; }
    else {
        std::cerr << "ERROR: Unrecognized executable name " << argv[0] << ".\n";
        exit(1);
    }
    qsolve_api->set_options(argc-1, argv);
    qsolve_api->read(argv[argc-1]);
    qsolve_api->compute();
    qsolve_api->write(argv[argc-1]);

    delete qsolve_api;
    return 0;
}   
