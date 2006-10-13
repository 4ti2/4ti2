/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.

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

#include "BasicCompletion.h"
#include "BinomialSetStream.h"
#include "Generation.h"
#include "Globals.h"

#include <iostream>
#include <iomanip>

//#define DEBUG_4ti2(X) X
#include "Debug.h"

using namespace _4ti2_;

BasicCompletion::BasicCompletion()
{
    name = "(F)";
}

BasicCompletion::~BasicCompletion()
{
}

bool
BasicCompletion::algorithm(BinomialSet& bs)
{
    DEBUG_4ti2(*out << "Before auto reduction:\n" << bs << "\n";)
    bs.auto_reduce_once();
    DEBUG_4ti2(*out << "After auto reduction:\n" << bs << "\n";)

    long int num_iterations = 0;
    Index i = 0;
    while(i < bs.get_number())
    {
        if (num_iterations % Globals::output_freq== 0)
        {
            *out << "\r" << Globals::context << name;
            *out << " Size: " << std::setw(6) << bs.get_number();
            *out << ", Index: " << std::setw(6) << i << std::flush;
            DEBUG_4ti2(*out << "\n";)
            DEBUG_4ti2(*out << bs << "\n";)
        }
        gen->generate(bs, i, bs);
        ++num_iterations;
        ++i;
        if (num_iterations % Globals::auto_reduce_freq == 0 && num_iterations != 0)
        {
            bs.auto_reduce_once(i);
        }
    }
    DEBUG_4ti2(*out << "Before minimal:\n" << bs << "\n";)
    bs.minimal();
    DEBUG_4ti2(*out << "After minimal:\n" << bs << "\n";)
    bs.reduced();
    DEBUG_4ti2(*out << "Final:\n" << bs << "\n";)
    return true;
}
