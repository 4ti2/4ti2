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

#include "groebner/minimize_main.h"
#include "groebner/Vector.h"
#include "groebner/VectorStream.h"
#include "groebner/VectorArray.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/BitSet.h"
#include "groebner/BitSetStream.h"
#include "groebner/Feasible.h"
#include "groebner/FeasibleStream.h"
#include "groebner/GeneratingSet.h"
#include "groebner/GroebnerBasis.h"
#include "groebner/Minimize.h"
#include "groebner/MinimizeOptions.h"
#include "groebner/Optimise.h"

#include <string>
#include <iostream>
#include <fstream>

using namespace _4ti2_;

int
_4ti2_::minimize_main(int argc, char **argv)
{
    MinimizeOptions::instance()->process_options(argc, argv);

    // Read in the sets of fibers.
    Feasible* feasible = input_Feasible(MinimizeOptions::instance()->filename.c_str());

    // Read in the file with the cost vectors.
    std::string cost_filename(MinimizeOptions::instance()->filename + ".cost");
    VectorArray* cost = input_VectorArray(feasible->get_dimension(), cost_filename.c_str());
    if (cost == 0 || cost->get_number() != 1)
    {
        std::cerr << "INPUT ERROR: There should be a single cost function.\n";
        exit(1);
    }

    if (feasible->get_rhs() == 0)
    {
        std::cerr << "INPUT ERROR: Could not find fiber file.\n";
        exit(1);
    }

    Vector sol(*feasible->get_rhs());
    Optimise opt;
    opt.compute(*feasible, (*cost)[0], sol);

    // Output the optimal solutions.
    std::string opt_filename(MinimizeOptions::instance()->filename + ".min");
    std::ofstream opt_file(opt_filename.c_str());
    output(opt_file, sol);

    delete feasible;
    delete cost;

    return 0;
}
