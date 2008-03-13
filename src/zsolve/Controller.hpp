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

#ifndef _4ti2_zsolve__Controller_
#define _4ti2_zsolve__Controller_

#include "zsolve/LinearSystem.hpp"
#include "zsolve/Lattice.hpp"
#include "zsolve/Timer.h"
#include "zsolve/Algorithm.hpp"

template <typename T> class Algorithm;

template <typename T> class Controller
{
public:
    virtual void log_system (LinearSystem <T> * system) = 0;
    virtual void log_homogenized_system (LinearSystem <T> * system) = 0;
    virtual void log_lattice (Lattice <T> * system) = 0;
    virtual void log_variable_start (size_t variable, size_t vectors) = 0;
    virtual void log_variable_end (size_t variable, size_t vectors) = 0;
    virtual void log_sum_start (size_t variable, const T& sum, size_t vectors) = 0;
    virtual void log_sum_end (size_t variable, const T& sum, size_t vectors) = 0;
    virtual void log_norm_start (size_t variable, const T& sum, const T& norm, size_t vectors) = 0;
    virtual void log_norm_end (size_t variable, const T& sum, const T& norm, size_t vectors) = 0;
    virtual void log_result (size_t inhoms, size_t homs, size_t free) = 0;
    virtual void log_status (size_t variable, const T& sum, const T& max_sum, const T& norm, size_t vectors, int backup_frequency, Timer& timer) = 0;
    virtual void log_resume (size_t variables, size_t variable, const T& sum, const T& norm, size_t vectors) = 0;
    virtual void log_maxnorm (Algorithm <T> * algorithm, bool final) = 0;

    virtual void save_lattice (Lattice <T> * system) = 0;
    virtual void backup_data (Lattice <T> & lattice, size_t current, const T& sum, const T& norm, bool symmetric) = 0;
    virtual void read_backup (std::ifstream& in) = 0;

    virtual ~Controller () {};
};

#endif
