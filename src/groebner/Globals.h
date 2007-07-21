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

#ifndef _4ti2__Globals_
#define _4ti2__Globals_

#include <string>
#include <iostream>

namespace _4ti2_
{

// 4ti2 Global settings.
class Globals
{
public:
    typedef enum {DEFAULT, NORMAL, WEIGHTED, GEBAUER_MOELLER} Algorithm;
    typedef enum {HYBRID, SATURATION, PROJECT_AND_LIFT, MAXMIN} Generation;
    typedef enum {WEIGHT, IP, LP, NONE} Truncation;

    static Algorithm algorithm;
    static Generation generation;
    static std::string context;
    static std::string exec;
    static int output_freq;
    static int auto_reduce_freq;
    static Truncation truncation;
    static bool minimal;
    static bool criteria;
    static int norm;
};

extern std::ostream* out;

} // namespace _4ti2_

#endif
