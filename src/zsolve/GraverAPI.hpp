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

#ifndef _4ti2_zsolve__GraverAPI_
#define _4ti2_zsolve__GraverAPI_

#include "4ti2/4ti2xx.h"
#include "zsolve/ZSolveAPI.hpp"

namespace _4ti2_zsolve_ {

template <class T>
class GraverAPI : public ZSolveAPI<T> {
public:
    GraverAPI();

    virtual void write(const char* project);

    virtual _4ti2_matrix* get_matrix(const char* name);

protected:
    virtual void check_consistency();

    // Extract the output after running the algorithm.
    virtual void extract_results(Algorithm <T>* algorithm);
};

template <class T>
GraverAPI<T>::GraverAPI()
{
    ZSolveAPI<T>::free_default = false;
    ZSolveAPI<T>::lower_default = 1;
    ZSolveAPI<T>::upper_default = -1;
}

template <class T>
void
GraverAPI<T>::check_consistency()
{
    ZSolveAPI<T>::check_consistency();

    if (ZSolveAPI<T>::rhs) {
        throw IOException ("No `rhs' allowed for `graver' executable. Use `zsolve' instead!\n");
    }
    if (ZSolveAPI<T>::rel) {
        throw IOException ("No `rel' allowed for `graver' executable. Use `zsolve' instead.");
    }
}

template <class T>
void
GraverAPI<T>::write(const char* project_c_str)
{
    std::string project(project_c_str);   

    if (ZSolveAPI<T>::zhom) { ZSolveAPI<T>::zhom->write((project + ".gra").c_str()); }
    if (ZSolveAPI<T>::zfree && ZSolveAPI<T>::zfree->data.height() > 0) {
        ZSolveAPI<T>::zfree->write((project + ".zfree").c_str());
    }
}

template <class T>
_4ti2_matrix*
GraverAPI<T>::get_matrix(const char* name)
{
    if (!strcmp(name, "gra")) { return ZSolveAPI<T>::zhom; }
    return ZSolveAPI<T>::get_matrix(name);
}

template <class T>
void
GraverAPI<T>::extract_results(Algorithm <T>* algorithm)
{
    delete ZSolveAPI<T>::zhom;
    ZSolveAPI<T>::zhom = new VectorArrayAPI <T> (0, algorithm->get_result_variables ());
    algorithm->extract_graver_results (ZSolveAPI<T>::zhom->data);
}

} // namespace _4ti2_zsolve_

#endif
