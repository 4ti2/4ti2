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
#include "qsolve/VectorArrayAPI.h"
#include "qsolve/VectorArrayStream.h"
#include "qsolve/CircuitsAPI.h"
#include "qsolve/Globals.h"
#include "qsolve/Debug.h"

using namespace _4ti2_;

template <class T>
CircuitsAPI<T>::CircuitsAPI()
    : QSolveAPI<T>()
{
    QSolveAPI<T>::sign_default = _4ti2_DB;
    QSolveAPI<T>::rel_default = _4ti2_EQ;
}

template <class T>
CircuitsAPI<T>::~CircuitsAPI()
{
}

template <class T>
void
CircuitsAPI<T>::post_compute()
{
    QSolveAPI<T>::ray.data.sort();
    VectorT<T> zero(QSolveAPI<T>::cir.data.get_size(), 0);
    for (int i = 0; i < QSolveAPI<T>::cir.data.get_number(); ++i) {
        if (QSolveAPI<T>::cir.data[i] < zero) { QSolveAPI<T>::cir.data[i].muleq(-1); }
    }
    QSolveAPI<T>::cir.data.sort();
    QSolveAPI<T>::qfree.data.sort();
    // TODO: Should we transfer the rays into the circuits?
    QSolveAPI<T>::cir.data.transfer(QSolveAPI<T>::ray.data, 0, QSolveAPI<T>::ray.data.get_number(), 0);
}

template <class T>
void
CircuitsAPI<T>::write_usage()
{
    std::cerr << "Usage: circuits [options] <PROJECT>\n\n";
    std::cerr << "Computes the circuits of a cone.\n";
    write_input_files();
    write_output_files();
    QSolveAPI<T>::write_options();
}

template <class T>
void
CircuitsAPI<T>::write_input_files()
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

template <class T>
void
CircuitsAPI<T>::write_output_files()
{
    std::cerr << "\
Output Files:\n\
  PROJECT.cir         The circuits of the cone.\n\
  PROJECT.qfree       A basis for the linear subspace of the cone.\n\
                      If this file does not exist then the linear subspace \n\
                      is trivial.\n\n";
}

template <class T>
void
CircuitsAPI<T>::write(const char* basename_c_str)
{
    std::string basename(basename_c_str);

    std::string cir_filename(basename + ".cir");
    QSolveAPI<T>::cir.write(cir_filename.c_str());

    std::string qfree_filename(basename + ".qfree");
    QSolveAPI<T>::qfree.write(qfree_filename.c_str());
}
