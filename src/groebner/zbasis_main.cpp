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

#include "zbasis_main.h"
#include "VectorArray.h"
#include "VectorArrayStream.h"
#include "LatticeBasis.h"
#include "BasicOptions.h"

#include <string>
#include <iostream>
#include <fstream>

using namespace _4ti2_;

int
_4ti2_::zbasis_main(int argc, char **argv)
{
    BasicOptions::instance()->process_options(argc, argv);

    // Read in the file with the matrix.
    std::string matrix_filename(BasicOptions::instance()->filename.c_str());
    std::ifstream matrix_file(matrix_filename.c_str());
    int m=0,n=0;
    if (!matrix_file.good())
    {
        std::cerr << "File not found: " << matrix_filename << std::endl;
        exit(1);
    }
    matrix_file >> m >> n;
    VectorArray matrix(m, n);
    matrix_file >> matrix;

    VectorArray zbasis(0, n);
    lattice_basis(matrix, zbasis);

    std::string zbasis_filename(matrix_filename + ".lat");
    output(zbasis_filename.c_str(), zbasis);

    return 0;
}
