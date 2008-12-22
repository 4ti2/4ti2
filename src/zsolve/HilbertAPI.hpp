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

#ifndef _4ti2_zsolve__HilbertAPI_
#define _4ti2_zsolve__HilbertAPI_

#include "4ti2/4ti2xx.h"
#include "zsolve/ZSolveAPI.hpp"

namespace _4ti2_zsolve_ {

template <class T>
class HilbertAPI : public ZSolveAPI<T> {
public:
    HilbertAPI();

    virtual void write(const char* project);

    virtual _4ti2_matrix* get_matrix(const char* name);

protected:
    virtual void check_consistency();

    // Extract the output after running the algorithm.
    virtual void extract_results(Algorithm <T>* algorithm);
};

template <class T>
HilbertAPI<T>::HilbertAPI()
{
    ZSolveAPI<T>::free_default = false;
    ZSolveAPI<T>::lower_default = 0;
    ZSolveAPI<T>::upper_default = -1;
}

template <class T>
void
HilbertAPI<T>::check_consistency()
{
    ZSolveAPI<T>::check_consistency();

    if (ZSolveAPI<T>::rhs) {
        throw IOException ("No `rhs' allowed for `hilbert' executable. Use `zsolve' instead!\n");
    }
    if (ZSolveAPI<T>::rel) {
        throw IOException ("No `rel' allowed for `hilbert' executable. Use `zsolve' instead.");
    }
    if (ZSolveAPI<T>::lb) {
        throw IOException ("No `lb' allowed for `hilbert' executable. Use `zsolve' or `graver' instead.");
    }
    if (ZSolveAPI<T>::sign) {
        for (size_t i = 0; i < ZSolveAPI<T>::sign->data.width(); ++i) {
            if (ZSolveAPI<T>::sign->data[0][i] == 2) {
                throw IOException ("Graver components are not allowed for `hilbert' executable. Use `zsolve' or `graver' instead.");
            }
        }
    }
}

template <class T>
void
HilbertAPI<T>::write(const char* project_c_str)
{
    std::string project(project_c_str);   

    if (ZSolveAPI<T>::zhom) { ZSolveAPI<T>::zhom->write((project + ".hil").c_str()); }
    if (ZSolveAPI<T>::zfree && ZSolveAPI<T>::zfree->data.height() > 0) {
        ZSolveAPI<T>::zfree->write((project + ".zfree").c_str());
    }
}

template <class T>
_4ti2_matrix*
HilbertAPI<T>::get_matrix(const char* name)
{
    if (!strcmp(name, "hil")) { return ZSolveAPI<T>::zhom; }
    return ZSolveAPI<T>::get_matrix(name);
}

template <class T>
void
HilbertAPI<T>::extract_results(Algorithm <T>* algorithm)
{
    delete ZSolveAPI<T>::zhom;
    ZSolveAPI<T>::zhom = new VectorArrayAPI <T> (0, algorithm->get_result_variables ());
    algorithm->extract_hilbert_results (ZSolveAPI<T>::zhom->data);
}

} // namespace _4ti2_zsolve_

#endif
