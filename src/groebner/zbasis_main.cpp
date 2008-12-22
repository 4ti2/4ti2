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

#include "groebner/zbasis_main.h"
#include "groebner/VectorArray.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/LatticeBasis.h"
#include "groebner/BasicOptions.h"
#include "groebner/Globals.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace _4ti2_;

int
_4ti2_::zbasis_main(int argc, char **argv)
{
    BasicOptions::instance()->process_options(argc, argv);

    print_banner();

    // Read in the file with the matrix.
    std::string project_filename(BasicOptions::instance()->filename);
    VectorArray* project = input_VectorArray(project_filename.c_str());
    std::string matrix_filename(project_filename + ".mat");
    VectorArray* matrix = input_VectorArray(matrix_filename.c_str());
    if (matrix != 0 && project != 0)
    {
        std::cerr << "Input Error: Both " << project_filename << " and ";
        std::cerr << matrix_filename << " exist.\n";
        std::cerr << "Input Error: Only one of them allowed (preferably ";
        std::cerr << matrix_filename << ").\n";
        exit(1);
    }
    if (project != 0)
    {
        std::cout << "WARNING: Please specify the matrix in the file '";
        std::cout << matrix_filename << "' instead of '";
        std::cout << project_filename << "'.\n";
        matrix = project;
    }

    VectorArray zbasis(0, matrix->get_size());
    lattice_basis(*matrix, zbasis);

    std::string zbasis_filename(project_filename + ".lat");
    output(zbasis_filename.c_str(), zbasis);

    return 0;
}
