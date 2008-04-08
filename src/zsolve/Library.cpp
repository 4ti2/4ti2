/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Matthias Walter.

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

#include "zsolve.h"

#include "Vector.hpp"
#include "VectorArray.hpp"
#include "Variables.hpp"
#include "Relation.hpp"
#include "LinearSystem.hpp"
#include "Lattice.hpp"
#include "Options.h"
#include "Exception.h"
#include "Algorithm.hpp"
#include "DefaultController.hpp"

using namespace _4ti2_zsolve_;

extern "C"
{
    // 32 bit system

    ZSolveSystem32 zsolve_system_create_32 (int height, int width)
    {
	VectorArray <int32_t> array (height, width);
	int32_t* rhs = create_vector <int32_t> (height);

	LinearSystem <int32_t>* system = new LinearSystem <int32_t> (array, rhs, true, 1, -1);
	delete_vector <int32_t> (rhs);

	return (ZSolveSystem32) system;
    }

    void zsolve_system_set_matrix_32 (ZSolveSystem32 system, int column, int row, int32_t value)
    {
	LinearSystem <int32_t>& ls = *(LinearSystem <int32_t>*) system;
	VectorArray <int32_t>& matrix = ls.matrix ();
	matrix[row][column] = value;
    }

    void zsolve_system_set_rhs_32 (ZSolveSystem32 system, int row, int32_t value)
    {
	LinearSystem <int32_t>& ls = *(LinearSystem <int32_t>*) system;
	ls.rhs ()[row] = value;
    }

    void zsolve_system_set_rel_32 (ZSolveSystem32 system, int row, int relation)
    {
	LinearSystem <int32_t>& ls = *(LinearSystem <int32_t>*) system;
        Relation <int32_t>& rel = ls.get_relation (row);
        rel.set ((Relation <int32_t>::RelationType) relation);
    }

    void zsolve_system_set_free_32 (ZSolveSystem32 system, int column)
    {
	LinearSystem <int32_t>& ls = *(LinearSystem <int32_t>*) system;
	ls.get_variable (column).set (true);
    }

    void zsolve_system_set_hilbert_32 (ZSolveSystem32 system, int column)
    {
        LinearSystem <int32_t>& ls = *(LinearSystem <int32_t>*) system;
        ls.get_variable (column).set (false, 0, -1);
    }

    void zsolve_system_set_graver_32 (ZSolveSystem32 system, int column)
    {
        LinearSystem <int32_t>& ls = *(LinearSystem <int32_t>*) system;
        ls.get_variable (column).set (false, 1, -1);
    }

    void zsolve_system_set_bounded_32 (ZSolveSystem32 system, int column, int32_t lower, int32_t upper)
    {
        LinearSystem <int32_t>& ls = *(LinearSystem <int32_t>*) system;
        ls.get_variable (column).set (false, lower, upper);
    }

    void zsolve_system_delete_32 (ZSolveSystem32 system)
    {
	delete ((LinearSystem <int32_t>*) system);
    }    

    void zsolve_system_print_32 (ZSolveSystem32 system)
    {
	LinearSystem <int32_t>& ls = *(LinearSystem <int32_t>*) system;
	std::cout << ls << std::endl;
    }

    // 32 bit matrix

    ZSolveMatrix32 zsolve_matrix_create_32 (int height, int width)
    {
	VectorArray <int32_t>* array = new VectorArray <int32_t> (height, width);
	return (ZSolveMatrix32) array;
    }

    void zsolve_matrix_set_32 (ZSolveMatrix32 matrix, int column, int row, int32_t value)
    {
	VectorArray <int32_t>& array = *(VectorArray <int32_t>*) matrix;
	array[row][column] = value;
    }

    int32_t zsolve_matrix_get_32 (ZSolveMatrix32 matrix, int column, int row)
    {
	VectorArray <int32_t>& array = *(VectorArray <int32_t>*) matrix;
	return array[row][column];
    }

    int zsolve_matrix_width_32 (ZSolveMatrix32 matrix)
    {
	VectorArray <int32_t>& array = *(VectorArray <int32_t>*) matrix;
        return array.width ();
    }

    int zsolve_matrix_height_32 (ZSolveMatrix32 matrix)
    {
	VectorArray <int32_t>& array = *(VectorArray <int32_t>*) matrix;
        return array.height ();
    }

    void zsolve_matrix_delete_32 (ZSolveMatrix32 matrix)
    {
	delete ((VectorArray <int32_t>*) matrix);
    }

    void zsolve_matrix_print_32 (ZSolveMatrix32 matrix)
    {
	VectorArray <int32_t>& array = *(VectorArray <int32_t>*) matrix;
	std::cout << array << std::endl;
    }

    // 32 bit state
    
    ZSolveState32 zsolve_state_create_32 (ZSolveSystem32 system)
    {
	LinearSystem <int32_t>* ls = (LinearSystem <int32_t>*) system;
	Algorithm <int32_t>* algorithm = new Algorithm <int32_t> (ls, NULL);

	return (ZSolveState32) algorithm;
    }

    void zsolve_state_compute_32 (ZSolveState32 state, ZSolveMatrix32* inhom, ZSolveMatrix32* hom, ZSolveMatrix32* free)
    {
	Algorithm <int32_t>* algorithm = (Algorithm <int32_t>*) state;
	algorithm->compute ();

	int variables = algorithm->get_result_variables ();

	VectorArray <int32_t>* inhoms = new VectorArray <int32_t> (variables);
	VectorArray <int32_t>* homs = new VectorArray <int32_t> (variables);
	VectorArray <int32_t>* frees = new VectorArray <int32_t> (variables);

	algorithm->extract_zsolve_results (*inhoms, *homs, *frees);

	*inhom = (ZSolveMatrix32) (inhoms);
	*hom = (ZSolveMatrix32) (homs);
	*free = (ZSolveMatrix32) (frees);
    }

    void zsolve_state_delete_32 (ZSolveState32 state)
    {
	delete ((Algorithm <int32_t>*) state);
    }
}

