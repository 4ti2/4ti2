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
#include "qsolve/QSolveAlgorithm.h"
#include "qsolve/VectorArrayStream.h"
#include "qsolve/RaysAPI.h"
#include "qsolve/Globals.h"
#include "qsolve/Debug.h"

using namespace _4ti2_;

RaysAPI::RaysAPI()
    : QSolveAPI()
{
    QSolveAPI::sign_default = _4ti2_LB;
    QSolveAPI::rel_default = _4ti2_EQ;
}

RaysAPI::~RaysAPI()
{
}

_4ti2_matrix*
RaysAPI::get_matrix(const char* name)
{
    if (!strcmp(name, "mat")) { return mat; }
    if (!strcmp(name, "sign")) { return sign; }
    if (!strcmp(name, "rel")) { return rel; }
    if (!strcmp(name, "ray")) { return ray; }
    if (!strcmp(name, "qfree")) { return qfree; }
    std::cerr << "ERROR: Unrecognised mat type " << name << ".\n";
    return 0;
}

void
RaysAPI::post_compute()
{
    //QSolveAPI::ray.data.sort();
    //QSolveAPI::qfree.data.sort();
}

void
RaysAPI::write_usage()
{
    std::cerr << "Usage: rays [options] <PROJECT>\n\n";
    std::cerr << "Computes the extreme rays of a cone.\n";
    write_input_files();
    write_output_files();
    QSolveAPI::write_options();
}

void
RaysAPI::write_input_files()
{
    std::cerr << "\
Input Files:\n\
  PROJECT.mat         A matrix (compulsory).\n\
  PROJECT.sign        The sign constraints of the variables ('1' means\n\
                      non-negative, '0' means a free variable, and '2' means\n\
                      both non-negative and non-positive).\n\
                      It is optional, and the default is all non-negative.\n\
  PROJECT.rel         The relations on the matrix rows ('<','>','=').\n\
                      It is optional and the default is all '='.\n\
                      The mat must be given with this file.\n";
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
    QSolveAPI::ray->write(ray_filename.c_str());

    std::string qfree_filename(basename + ".qfree");
    QSolveAPI::qfree->write(qfree_filename.c_str());
}
