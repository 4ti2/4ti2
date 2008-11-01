/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Matthias Walter, Peter Malkin.

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

#include "../banner.h"

#include "4ti2/4ti2.h"
#include "4ti2/4ti2xx.h"

#include "zsolve/Vector.hpp"
#include "zsolve/VectorArray.hpp"
#include "zsolve/Variables.hpp"
#include "zsolve/Relation.hpp"
#include "zsolve/LinearSystem.hpp"
#include "zsolve/Lattice.hpp"
#include "zsolve/Options.h"
#include "zsolve/Exception.h"
#include "zsolve/Algorithm.hpp"
#include "zsolve/DefaultController.hpp"

#include <iostream>
#include <fstream>

using namespace _4ti2_zsolve_;

int main (int argc, char **argv)
{
    Options options (argc, argv);
    if (options.verbosity () != 0) { 
        Options::print_banner ();
        options.print_precision ();
    }

    try {
        _4ti2_state* state = 0;
        if (options.hilbert()) { state = _4ti2_hilbert_create_state(options.precision()); }
        else if (options.graver()) { state = _4ti2_graver_create_state(options.precision()); } 
        else { state = _4ti2_zsolve_create_state(options.precision()); }
        state->set_options(argc, argv);
        state->read(options.project().c_str());
        state->compute();
        state->write(options.project().c_str());
    }
    catch (PrecisionException e)
    {
	    std::cerr << "Results were near maximum precision (" << e.precision () << "bit).\n";
        std::cerr << "Please restart with higher precision!" << std::endl;
	    return 1;
    }
    catch (IOException e)
    {
        std::cerr << e;
        return 1;
    }
}

