/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2008 4ti2 team.
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

#include <cstring>
#include <iostream>
#include "4ti2/4ti2.h"
#include "groebner/QSolveAlgorithm.h"
#include "groebner/VectorArrayAPI.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/LatticeBasis.h"
#include "groebner/RaysAPI.h"

using namespace _4ti2_;

RaysAPI::RaysAPI()
{
}

RaysAPI::~RaysAPI()
{
}

void
RaysAPI::compute()
{
    // Consistency and default value checking.
    if (!matrix && !lat) {
        std::cerr << "ERROR: No matrix specified.\n";
        exit(1);
    }
    if (!matrix) {
        matrix = new VectorArrayAPI(0, lat->get_num_cols());
        lattice_basis(lat->data, matrix->data);
    }
    if (!lat) {
        lat = new VectorArrayAPI(0, matrix->get_num_cols());
        lattice_basis(matrix->data, lat->data);
    }
    if (!sign) {
        sign = new VectorArrayAPI(1, matrix->get_num_cols());
        for (Index i = 0; i < sign->get_num_cols(); ++i) { sign->data[0][i] = 1; }
    }
    if (!rel) {
        rel = new VectorArrayAPI(1, matrix->get_num_cols());
        for (Index i = 0; i < rel->get_num_cols(); ++i) { rel->data[0][i] = 0; }
    }
    assert(sign->get_number() == 1);
    assert(matrix->get_num_cols() == sign->get_num_cols());

    std::cout << "Matrix:\n";
    matrix->write(std::cout);
    std::cout << "Sign:\n";
    sign->write(std::cout);
    std::cout << "Rel:\n";
    rel->write(std::cout);
    std::cout << "Lat:\n";
    lat->write(std::cout);

    // Delete previous computation.
    delete ray; delete qfree;
    ray = new VectorArrayAPI(0, matrix->get_num_cols());
    ray->data.insert(lat->data);
    qfree = new VectorArrayAPI(0, matrix->get_num_cols());

    QSolveAlgorithm alg(algorithm, order);
    alg.compute(matrix->data, ray->data, qfree->data, rel->data[0], sign->data[0]); 

    ray->data.sort();
    qfree->data.sort();
}

void
RaysAPI::write_usage()
{
    std::cerr << "Usage: rays [options] <PROJECT>\n\n";
    std::cerr << "Computes the extreme rays of a cone.\n";
    write_input_files();
    write_output_files();
    write_options();
}

void
RaysAPI::write_output_files()
{
    std::cerr << "\
Output Files:\n\
  PROJECT.ray         The extreme rays of the cone.\n\
  PROJECT.qfree       A basis for the linear subspace of the cone.\n\
                      If this file does not exist then the linear subspace \n\
                      is trivial.\n\n";
}

void
RaysAPI::write(const char* basename_c_str)
{
    std::string basename(basename_c_str);

    std::string ray_filename(basename + ".ray");
    ray->write(ray_filename.c_str());

    std::string qfree_filename(basename + ".qfree");
    qfree->write(qfree_filename.c_str());
}
