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

#include "groebner/DataType.h"
#include "groebner/Globals.h"
#include "banner.h"
#include <iostream>

using namespace _4ti2_;

Globals::Algorithm  Globals::algorithm = DEFAULT;
Globals::Generation Globals::generation = HYBRID;
std::string Globals::context;
std::string Globals::exec;
int Globals::output_freq = 1000;
int Globals::auto_reduce_freq = 2500;
Globals::Truncation Globals::truncation = Globals::WEIGHT;
std::ostream* _4ti2_::out = &std::cout;
std::ostream* _4ti2_::err = &std::cerr;
bool Globals::minimal = true;
bool Globals::criteria = false;
int Globals::norm = 1;

void
_4ti2_::print_banner()
{
    *out << FORTY_TWO_BANNER;

#if defined(_4ti2_INT32_) || defined(_4ti2_INT64_)
    *out << "Using " << sizeof(IntegerType)*CHAR_BIT << " bit integers.\n";
# if !defined(HAVE_TRAPV_LONG_LONG)
    *err << "WARNING: Overflow detection is not available on this architecture/compiler.\n"
	      << "WARNING: To guarantee correct results, run 4ti2 with arbitrary precision\n"
	      << "WARNING: by using the option `-parb'\n";
# endif
#elif defined(_4ti2_GMP_)
    *out << "Using arbitrary precision integers.\n";
#endif
}
