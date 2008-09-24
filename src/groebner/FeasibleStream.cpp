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

#include "groebner/FeasibleStream.h"
#include "groebner/VectorStream.h"
#include "groebner/VectorArrayStream.h"
#include "groebner/BitSetStream.h"
#include "groebner/Globals.h"

using namespace _4ti2_;

std::ostream&
_4ti2_::operator<<(std::ostream& out, Feasible& f)
{
    out << "Feasible:\n";
    out << "Matrix:\n" << f.get_matrix();
    out << "Basis:\n" << f.get_basis();
    out << "URS:\n" << f.get_urs() << "\n";
    out << "Bounded:\n" << f.get_bnd() << "\n";
    out << "Unbounded:\n" << f.get_unbnd() << "\n";
    out << "GRADING:\n" << f.get_grading() << "\n";
    out << "RAY:\n" << f.get_ray() << "\n";
    if (f.get_weights() != 0)
    {
        out << "WEIGHTS:\n" << *f.get_weights();
    }
    if (f.get_max_weights() != 0)
    {
        out << "MAX WEIGHTS: " << *f.get_max_weights() << "\n";
    }

    return out;
}

Feasible*
_4ti2_::input_Feasible(const char* filename)
{
    // Read in the file with the matrix.
    std::string project_filename(filename);
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
        *err << "WARNING: Please specify the matrix in the file '";
        *err << matrix_filename << "' instead of '";
        *err << project_filename << "'.\n";
        matrix = project;
    }

    // Read in the file with the basis.
    std::string basis_filename(project_filename + ".lat");
    VectorArray* basis = input_VectorArray(basis_filename.c_str());

    // There should either be a matrix file or a basis file.
    if (matrix == 0 && basis == 0)
    {
        std::cerr << "Input Error: Could not find either " << matrix_filename;
        std::cerr << " or " << basis_filename << ".\n";
        exit(1);
    }
    if (matrix != 0 && basis != 0 && matrix->get_size() != basis->get_size())
    {
        std::cerr << "Input Error: Size mismatch in files " << basis_filename;
        if (project != 0) { std::cerr << " and " << project_filename << ".\n"; }
        else { std::cerr << " and " << matrix_filename << ".\n"; }
        exit(1);
    }

    int dim;
    if (matrix != 0) { dim = matrix->get_size(); }
    else { dim = basis->get_size(); }

    // Read in the file with the sign.
    // Defaults to no urs variables if no file found.
    std::string sign_filename(project_filename + ".sign");
    VectorArray* sign = input_VectorArray(dim, sign_filename.c_str());
    BitSet urs(dim);
    if (sign != 0)
    {
        if (sign->get_number() != 1)
        {
            std::cerr << "Input Error: Expected one vector in ";
            std::cerr << sign_filename << "\n";
            exit(1);
        }
        for (int i = 0; i < dim; ++i) 
        {
            IntegerType value = (*sign)[0][i];
            if (value == 0) { urs.set(i); }
            else if (value == 1) { } // Nonnegative variable.
            else if (value == 2 || value == -1)
            {
                std::cerr << "Input Error in file " << sign_filename << "\n";
                std::cerr << "The value " << value << " is not yet supported in sign vector.\n";
                exit(1);
            }
            else
            {
                std::cerr << "Input Error in file " << sign_filename << "\n";
                std::cerr << "Unsupport number " << value << " in sign vector.\n";
                exit(1);
            }
        }
    }

    // Read in the file with the weight vectors.
    std::string weights_filename(project_filename + ".weights");
    VectorArray* weights = input_VectorArray(dim, weights_filename.c_str());

    Vector* max_weights = 0;
    if (weights != 0)
    {
        // Read in the file with the maximum weight vectors.
        std::string max_weights_filename(project_filename + ".weights.max");
        VectorArray* tmp_max_weights =
                input_VectorArray(weights->get_number(), max_weights_filename.c_str());
        if (tmp_max_weights == 0)
        {
            std::cerr << "Input Error: Could not find file " << max_weights_filename << ".\n";
            std::cerr << "Input Error: It is required with " << weights_filename << ".\n";
            exit(1);
        }
        if (tmp_max_weights->get_number() != 1)
        {
            std::cerr << "Input Error: Expected a single row matrix in file ";
            std::cerr << tmp_max_weights << "\n";
            exit(1);
        }
        max_weights = new Vector((*tmp_max_weights)[0]);
        delete tmp_max_weights;
    }

    // Read in the file with the zsol.
    std::string zsol_filename(project_filename + ".zsol");
    VectorArray* tmp_rhs = input_VectorArray(dim, zsol_filename.c_str());
    Vector* rhs = 0;
    if (tmp_rhs != 0)
    {
        if (tmp_rhs->get_number() != 1)
        {
            std::cerr << "Input Error: Expected a single row matrix in file ";
            std::cerr << zsol_filename << "\n";
            exit(1);
        }
        rhs = new Vector((*tmp_rhs)[0]);
        delete tmp_rhs;
    }

    // Check for the rel file.  It is not supported yet.
    std::string rel_filename(project_filename + ".rel");
    std::ifstream rel_file(rel_filename.c_str());
    if (rel_file.good())
    {
        std::cerr << "Input Error: The file " << rel_filename << " is not yet supported.\n";
        std::cerr << "Input Error: You need to move or remove it.\n";
        exit(1);
    }

    // Check for the rhs file.  It is not supported yet.
    std::string rhs_filename(project_filename + ".rhs");
    std::ifstream rhs_file(rhs_filename.c_str());
    if (rhs_file.good())
    {
        std::cerr << "Input Error: The file " << rhs_filename << " is not yet supported.\n";
        std::cerr << "Input Error: You need to move or remove it.\n";
        exit(1);
    }

    // Construct the feasible sets.
    Feasible* feasible = new Feasible(basis, matrix, &urs, rhs, weights, max_weights);

    delete basis;
    delete matrix;
    delete sign;
    delete max_weights;
    delete weights;
    delete rhs;

    return feasible;
}
