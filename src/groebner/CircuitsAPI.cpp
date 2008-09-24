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
#include "groebner/CircuitsAPI.h"
#include "groebner/Globals.h"
#include "groebner/Debug.h"

using namespace _4ti2_;

CircuitsAPI::CircuitsAPI()
{
}

CircuitsAPI::~CircuitsAPI()
{
}

void
CircuitsAPI::compute()
{
    print_banner();

    // Consistency and default value checking.
    // TODO: More consistency checking.
    if (!mat) {
        std::cerr << "ERROR: No constraint matrix specified.\n";
        exit(1);
    }
    if (!sign) {
        sign = new VectorArrayAPI(1, mat->get_num_cols());
        for (Index i = 0; i < sign->get_num_cols(); ++i) { sign->data[0][i] = 2; }
    }
    if (!rel) {
        rel = new VectorArrayAPI(1, mat->get_num_cols());
        for (Index i = 0; i < rel->get_num_cols(); ++i) { rel->data[0][i] = 0; }
    }
    assert(sign->get_number() == 1);
    assert(mat->get_num_cols() == sign->get_num_cols());

    DEBUG_4ti2(std::cout << "Matrix:\n";)
    DEBUG_4ti2(mat->write(std::cout);)
    DEBUG_4ti2(std::cout << "Sign:\n";)
    DEBUG_4ti2(sign->write(std::cout);)
    DEBUG_4ti2(std::cout << "Rel:\n";)
    DEBUG_4ti2(rel->write(std::cout);)

    // Delete previous computation.
    delete ray; delete cir; delete qhom; delete qfree;
    ray = new VectorArrayAPI(0, mat->get_num_cols());
    cir = new VectorArrayAPI(0, mat->get_num_cols());
    qhom = new VectorArrayAPI(0, mat->get_num_cols());
    qfree = new VectorArrayAPI(0, mat->get_num_cols());

    QSolveAlgorithm alg(algorithm, order);
    alg.compute(mat->data, ray->data, cir->data, qfree->data, rel->data[0], sign->data[0]); 

    ray->data.sort();
    cir->data.sort();
    qfree->data.sort();
    VectorArray::transfer(ray->data, 0, ray->data.get_number(), cir->data, 0);
}

void
CircuitsAPI::write_usage()
{
    std::cerr << "Usage: circuits [options] <PROJECT>\n\n";
    std::cerr << "Computes the circuits of a cone.\n";
    write_input_files();
    write_output_files();
    write_options();
}

void
CircuitsAPI::write_input_files()
{
    std::cerr << "\
Input Files:\n\
  PROJECT.mat         A matrix (compulsory).\n\
  PROJECT.sign        The sign constraints of the variables ('1' means\n\
                      non-negative, '0' means a free variable, and '2' means\n\
                      both non-negative and non-positive).\n\
                      It is optional, and the default is both.\n\
  PROJECT.rel         The relations on the matrix rows ('<','>','=').\n\
                      It is optional and the default is all '='.\n\
                      The mat must be given with this file.\n";
}

void
CircuitsAPI::write_output_files()
{
    std::cerr << "\
Output Files:\n\
  PROJECT.cir         The circuits of the cone.\n\
  PROJECT.qfree       A basis for the linear subspace of the cone.\n\
                      If this file does not exist then the linear subspace \n\
                      is trivial.\n\n";
}

void
CircuitsAPI::write(const char* basename_c_str)
{
    std::string basename(basename_c_str);

    std::string cir_filename(basename + ".cir");
    cir->write(cir_filename.c_str());

    std::string qfree_filename(basename + ".qfree");
    qfree->write(qfree_filename.c_str());
}
