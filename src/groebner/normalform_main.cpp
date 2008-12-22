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

#include "groebner/normalform_main.h"
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
#include "groebner/BasicOptions.h"
#include "groebner/Globals.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace _4ti2_;

int
_4ti2_::normalform_main(int argc, char **argv)
{
    BasicOptions::instance()->process_options(argc, argv);

    print_banner();

    // Read in the sets of fibers.
    Feasible* feasible = input_Feasible(BasicOptions::instance()->filename.c_str());

    // Read in the file with the groebner basis.
    std::string groebner_filename(BasicOptions::instance()->filename + ".gro");
    VectorArray* groebner = input_VectorArray(groebner_filename.c_str());
    // There should be a groebner basis.
    if (groebner == 0)
    {
        std::cerr << "ERROR: Could not find " << groebner_filename;
        std::cerr << "." << std::endl;
        exit(1);
    }

    // Read in the file with the cost vectors.
    std::string cost_filename(BasicOptions::instance()->filename + ".cost");
    VectorArray* cost =
            input_VectorArray(feasible->get_dimension(), cost_filename.c_str());
    if (cost == 0) { cost = new VectorArray(0,feasible->get_dimension()); }

    // Read in the file with the feasible points.
    std::string feas_filename(BasicOptions::instance()->filename + ".feas");
    VectorArray* feas
            = input_VectorArray(feasible->get_dimension(), feas_filename.c_str());
    // There should either be a feasible file.
    if (feas == 0)
    {
        std::cerr << "ERROR: Could not find " << feas_filename;
        std::cerr << "." << std::endl;
        exit(1);
    }

    // Compute the normal form.
    Minimize algorithm;
    algorithm.minimize(*feasible, *cost, *groebner, *feas);

    // Output the optimal solutions.
    std::string nf_filename(BasicOptions::instance()->filename + ".nf");
    std::ofstream nf_file(nf_filename.c_str());
    output(nf_file, *feas);

    delete feasible;
    delete cost;
    delete groebner;
    delete feas;

    return 0;
}
