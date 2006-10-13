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

#include <iostream>
#include <fstream>

#include "groebner_main.h"
#include "Vector.h"
#include "VectorStream.h"
#include "VectorArray.h"
#include "VectorArrayStream.h"
#include "BitSet.h"
#include "BitSetStream.h"
#include "Feasible.h"
#include "FeasibleStream.h"
#include "GeneratingSet.h"
#include "GroebnerBasis.h"
#include "Globals.h"
#include "Options.h"

#include <string>

using namespace _4ti2_;

int
_4ti2_::groebner_main(int argc, char **argv)
{
    Options::instance()->process_options(argc, argv);

    // Read in the sets of fibers.
    Feasible* feasible = input_Feasible(Options::instance()->filename.c_str());

    // Read in the file with the generating set.
    std::string gens_filename(Options::instance()->filename + ".mar");
    VectorArray* gens = input_VectorArray(feasible->get_dimension(), gens_filename.c_str());

    // Construct the generating set.
    // Do not bother to compute a minimal generating set.
    Globals::minimal = false;
    GeneratingSet gs(*feasible, gens);

    // Read in the file with the cost vectors.
    std::string cost_filename(Options::instance()->filename + ".cost");
    VectorArray* cost = input_VectorArray(feasible->get_dimension(), cost_filename.c_str());

    // Construct the Groebner basis.
    GroebnerBasis gb(gs, cost);
    delete cost;

    // Output the Groebner basis.
    std::string groebner_filename(Options::instance()->filename + ".gro");
    output(groebner_filename.c_str(), gb.get_groebner_basis());

    delete feasible;

    return 0;
}
