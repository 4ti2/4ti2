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
#include "qsolve/Options.h"

//#define DEBUG_4ti2(X) X
#include "qsolve/Debug.h"
#include "qsolve/Globals.h"

using namespace _4ti2_;

template <template <class> class API>
_4ti2_state*
create_state_precision(_4ti2_precision prec)
{
    switch (prec) {
    case _4ti2_PREC_INT_32:
        return new API<int32_t>();
    case _4ti2_PREC_INT_64:
        return new API<int64_t>();
    case _4ti2_PREC_INT_ARB:
#ifdef _4ti2_GMP_
        return new API<mpz_class>();
#else
        return 0;
#endif
    }
    return 0;
}

int
main(int argc, char **argv)
{
    if (argc == 1) {
        std::cerr << "Usage: 4ti2 <exec> ...\n";
        exit(1);
    }
    Options options(argc-1, argv+1);

    _4ti2_state* state = 0;
    std::string exec_name(argv[1]);
    if (exec_name == "qsolve") {
        state = create_state_precision<QSolveAPI>(options.precision);
    } else if (exec_name == "rays") {
        state = create_state_precision<RaysAPI>(options.precision);
    } else if (exec_name == "circuits") {
        state = create_state_precision<CircuitsAPI>(options.precision);
    } else {
        std::cerr << "ERROR: Unrecoginized executable name: " << argv[1] << "\n";
        exit(1);
    }

    state->set_options(argc-1, argv+1);
    state->read(options.filename.c_str());
    state->compute();
    state->write(options.filename.c_str());

    delete state;
    return 0;
}


