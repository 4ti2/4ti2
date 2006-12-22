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

#include "qsolve_main.h"
#include "Vector.h"
#include "VectorStream.h"
#include "VectorArray.h"
#include "VectorArrayStream.h"
#include "BitSet.h"
#include "BitSetStream.h"
#include "CircuitAlgorithm.h"
#include "RayAlgorithm.h"
#include "LatticeBasis.h"
#include "CircuitOptions.h"

//#define DEBUG_4ti2(X) X
#include "Debug.h"
#include "Globals.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace _4ti2_;

bool
input_Sign(const char* sign_filename, BitSet& rs, BitSet& cirs, BitSet& urs);
bool
input_Rel(const char* rel_filename, BitSet& eq, BitSet& gt, BitSet& lt);

int
_4ti2_::qsolve_main(int argc, char **argv)
{
    CircuitOptions::instance()->process_options(argc, argv);

    std::string rhs_filename(CircuitOptions::instance()->filename+".rhs");
    std::ifstream rhs_file(rhs_filename.c_str());
    if (rhs_file.good())
    {
        std::cerr << "ERROR: The file " << rhs_filename;
        std::cerr << " is not allowed with " << Globals::exec << ".\n";
        std::cerr << "ERROR: Remove it and call ";
        std::cerr << Globals::exec << " again.\n";
        exit(1);
    }

    // Read in the file with the matrix.
    std::string project_filename(CircuitOptions::instance()->filename);
    VectorArray* project = input_VectorArray(project_filename.c_str());
    std::string matrix_filename(project_filename + ".mat");
    VectorArray* matrix = input_VectorArray(matrix_filename.c_str());
    if (matrix != 0 && project != 0)
    {
        std::cerr << "INPUT ERROR: Both " << project_filename << " and ";
        std::cerr << matrix_filename << " exist.\n";
        std::cerr << "INPUT ERROR: Only one of them allowed (preferably ";
        std::cerr << matrix_filename << ").\n";
        exit(1);
    }
    if (project != 0)
    {
        std::cout << "WARNING: Please specify the matrix in '" << matrix_filename;
        std::cout << "' instead of '" << project_filename << "'.\n";
        matrix = project;
    }

    // Read in the file with the lattice basis.
    std::string lattice_filename(CircuitOptions::instance()->filename + ".lat");
    VectorArray* rays = input_VectorArray(lattice_filename.c_str());

    if (matrix == 0 && rays == 0)
    {
        std::cerr << "INPUT ERROR: Could not find either " << project_filename;
        std::cerr << " or " << lattice_filename << ".\n";
        exit(1);
    }
    if (matrix != 0 && rays != 0 && matrix->get_size() != rays->get_size())
    {
        std::cerr << "INPUT ERROR: Size mismatch in files " << lattice_filename;
        if (project != 0) { std::cerr << " and " << project_filename << ".\n"; }
        else { std::cerr << " and " << matrix_filename << ".\n"; }
        exit(1);
    }

    int dim = 0;
    if (matrix != 0) { dim = matrix->get_size(); }
    else if (rays != 0) { dim = rays->get_size(); }

    // Read in the sign of the components.
    std::string sign_filename(project_filename + ".sign");
    BitSet* urs = new BitSet(dim);
    BitSet* cirs = new BitSet(dim);
    BitSet* rs = new BitSet(dim);
    if (Globals::exec == "rays") { rs->one(); }
    else if (Globals::exec == "circuits") { cirs->one(); }
    else { urs->one(); }
    input_Sign(sign_filename.c_str(), *rs, *cirs, *urs);
    if (Globals::exec == "rays" && !cirs->empty())
    {
        std::cerr << "INPUT ERROR: Circuit components for `rays' executable.\n";
        std::cerr << "INPUT ERROR: Use the `circuits' executable instead.\n";
        exit(1);
    }

    if (matrix != 0)
    {
        // Read in the type of constraints.
        int rank = matrix->get_number();
        std::string rel_filename(project_filename + ".rel");
        BitSet eq(rank, true); // By default everything is an equality constraint.
        BitSet lt(rank);
        BitSet gt(rank);
        input_Rel(rel_filename.c_str(), eq, gt, lt);
        int extra = gt.count() + lt.count();
        if (extra > 0)
        {
            VectorArray* ext_matrix =
                    new VectorArray(matrix->get_number(), matrix->get_size()+extra, 0);
            VectorArray::lift(*matrix, 0, matrix->get_size(), *ext_matrix);
            int col = matrix->get_size();
            for (int i = 0; i < matrix->get_number(); ++i)
            {
                if (gt[i]) { (*ext_matrix)[i][col] = -1; ++col; }
                else if (lt[i]) { (*ext_matrix)[i][col] = 1; ++col; }
            }
            delete matrix; matrix = ext_matrix;
            BitSet* ext_urs = new BitSet(dim+extra);
            BitSet::extend(*urs, *ext_urs);
            BitSet* ext_cirs = new BitSet(dim+extra);
            BitSet::extend(*cirs, *ext_cirs);
            BitSet* ext_rs = new BitSet(dim+extra);
            BitSet::extend(*rs, *ext_rs);
            for (int i = dim; i < dim+extra; ++i) { ext_rs->set(i); }
            delete urs; urs = ext_urs;
            delete cirs; cirs = ext_cirs;
            delete rs; rs = ext_rs;
    
            delete rays; rays = 0;
        }

        if (rays == 0)
        {
            rays = new VectorArray(0, matrix->get_size());
            lattice_basis(*matrix, *rays);
        }
    }
    else //(matrix == 0)
    {
        matrix = new VectorArray(0, rays->get_size());
        lattice_basis(*rays, *matrix);
    }

    DEBUG_4ti2(*out << "MATRIX:\n" << *matrix << "\n";)
    DEBUG_4ti2(*out << "URS:\n" << *urs << "\n";)
    DEBUG_4ti2(*out << "CIRS:\n" << *cirs << "\n";)
    DEBUG_4ti2(*out << "RS:\n" << *rs << "\n";)

    VectorArray* subspace = new VectorArray(0, matrix->get_size());
    VectorArray* circuits = new VectorArray(0, matrix->get_size());
    if (!cirs->empty())
    {
        CircuitAlgorithm algorithm;
        algorithm.compute(*matrix, *rays, *circuits, *subspace, *rs, *cirs);
    }
    else
    {
        RayAlgorithm algorithm;
        algorithm.compute(*matrix, *rays, *subspace, *rs);
    }

    std::string output_filename(project_filename);
    if (Globals::exec == "rays")
    {
        output_filename += ".ray";
        rays->sort();
        std::ofstream file(output_filename.c_str());
        file << rays->get_number() << " " << dim << "\n";
        print(file, *rays, 0, dim);
    }
    else if (Globals::exec == "circuits")
    {
        output_filename += ".cir";
        rays->sort();
        circuits->sort();
        std::ofstream file(output_filename.c_str());
        file << rays->get_number()+circuits->get_number() << " " << dim << "\n";
        print(file, *rays, 0, dim);
        print(file, *circuits, 0, dim);
    }
    else //if (Globals::exec == "qsolve")
    {
        output_filename += ".qhom";
        rays->sort();
        circuits->sort();
        std::ofstream file(output_filename.c_str());
        file << rays->get_number()+2*circuits->get_number() << " " << dim << "\n";
        print(file, *rays, 0, dim);
        print(file, *circuits, 0, dim);
        circuits->mul(-1);
        print(file, *circuits, 0, dim);
    }

    if (!urs->empty())
    {
        subspace->sort();
        std::string subspace_filename(project_filename + ".qfree");
        std::ofstream subspace_file(subspace_filename.c_str());
        subspace_file << subspace->get_number() << " " << dim << "\n";
        print(subspace_file, *subspace, 0, dim);
    }

    delete matrix;
    delete rays;
    delete subspace;
    delete urs;
    delete cirs;
    delete rs;

    return 0;
}

bool
input_Sign(const char* sign_filename, BitSet& rs, BitSet& cirs, BitSet& urs)
{
    assert(rs.get_size() == cirs.get_size() && rs.get_size() == urs.get_size());
    int dim = rs.get_size();
    VectorArray* sign = input_VectorArray(dim, sign_filename);
    if (sign != 0)
    {
        rs.zero(); cirs.zero(); urs.zero();
        if (sign->get_number() != 1)
        {
            std::cerr << "INPUT ERROR: Expected one vector in ";
            std::cerr << sign_filename << "\n";
            exit(1);
        }
        for (int i = 0; i < dim; ++i) 
        {
            IntegerType value = (*sign)[0][i];
            if (value == 0) { urs.set(i); }
            else if (value == 1) { rs.set(i); } // Nonnegative component.
            else if (value == 2) { cirs.set(i); } // Circuit component.
            else if (value == -1)
            {
                std::cerr << "INPUT ERROR: In file " << sign_filename << "\n";
                std::cerr << "INPUT ERROR: the value " << value;
                std::cerr << " is not yet supported.\n";
                exit(1);
            }
            else
            {
                std::cerr << "INPUT ERROR in file " << sign_filename << "\n";
                std::cerr << "Unsupport number " << value << " in sign vector.\n";
                exit(1);
            }
        }
        delete sign;
    }
    DEBUG_4ti2(*out << "CIRCUIT COMPONENTS:\n" << cirs << "\n";)
    DEBUG_4ti2(*out << "RAY COMPONENTS:\n" << rs << "\n";)
    DEBUG_4ti2(*out << "URS COMPONENTS:\n" << urs << "\n";)
    if (sign == 0) { return false; }
    return true;
}

bool
input_Rel(const char* rel_filename, BitSet& eq, BitSet& gt, BitSet& lt)
{
    assert(eq.get_size() == lt.get_size() && eq.get_size() == gt.get_size());
    int rank = eq.get_size();
    std::ifstream file(rel_filename);
    if (!file.good()) { return 0; }
    eq.zero(); lt.zero(); gt.zero();
    int m, n;
    file >> m >> n;
    if (m != 1)
    {
        std::cerr << "INPUT ERROR: Expected one vector in ";
        std::cerr << rel_filename << "\n";
        exit(1);
    }
    if (n != rank)
    {
        std::cerr << "INPUT ERROR: Expected length of " << rank;
        std::cerr << " instead of " << n << " in " << rel_filename << "\n";
        exit(1);
    }
    char c;
    for (int i = 0; i < n; ++i)
    {
        file >> c;
        switch (c)
        {
        case '=':
            eq.set(i);
            break;
        case '<':
            lt.set(i);
            break;
        case '>':
            gt.set(i);
            break;
        default:
            std::cerr << "INPUT ERROR in file " << rel_filename << ".\n";
            std::cerr << "Unrecognised character " << c << "\n";
            std::cerr << "Character should be one of =, <, and >.\n";
        }
    }
    if (file.fail() || file.bad())
    {
        std::cerr << "INPUT ERROR: Badly formatted file " << rel_filename << ".\n";
        std::cerr << "INPUT ERROR: Check the size.\n";
        std::cerr << "INPUT ERROR: Check there are only the characters =, <, and >.\n";
        exit(1);
    }
}
