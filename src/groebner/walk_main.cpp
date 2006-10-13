/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.

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

#include "walk_main.h"
#include "Vector.h"
#include "VectorStream.h"
#include "VectorArray.h"
#include "VectorArrayStream.h"
#include "WalkOptions.h"
#include "WalkAlgorithm.h"
#include "Feasible.h"
#include "FeasibleStream.h"

#include <iostream>
#include <fstream>
#include <string>

// TODO: Memory allocation error handling.

using namespace _4ti2_;

int
_4ti2_::walk_main(int argc, char **argv)
{
    WalkOptions::instance()->process_options(argc, argv);

    // Read in the sets of fibers.
    Feasible* feasible = input_Feasible(WalkOptions::instance()->filename.c_str());

    // Read in the file with the old groebner basis.
    std::string groold_filename(WalkOptions::instance()->filename + ".gro.start");
    VectorArray* gro = input_VectorArray(feasible->get_dimension(), groold_filename.c_str());
    if (gro == 0)
    {
        std::cerr << "Input Error: could not find " << groold_filename << "\n";
        exit(1);
    }

    // Read in the file with the old cost vectors.
    std::string coststart_filename(WalkOptions::instance()->filename + ".cost.start");
    VectorArray* coststart = input_VectorArray(feasible->get_dimension(), coststart_filename.c_str());
    if (coststart == 0)
    {
        std::cerr << "Input Error: could not find " << groold_filename << "\n";
        exit(1);
    }

    // Read in the file with the new cost vectors.
    std::string cost_filename(WalkOptions::instance()->filename + ".cost");
    VectorArray* cost = input_VectorArray(feasible->get_dimension(), cost_filename.c_str());
    if (cost == 0)
    {
        std::cerr << "Input Error: could not find " << cost_filename << "\n";
        exit(1);
    }

    WalkAlgorithm algorithm;
    algorithm.compute(*feasible, *coststart, *gro, *cost);

    // Output the Groebner basis.
    std::string groebner_filename(WalkOptions::instance()->filename + ".gro");
    output(groebner_filename.c_str(), *gro);

    delete feasible;
    delete gro;
    delete cost;
    delete coststart;

    return 0;
}
