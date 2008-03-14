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

#ifndef _4ti2__Heuristics_
#define _4ti2__Heuristics_

#include <map>
#include "zsolve/BitSet.h"
#include "zsolve/LinearSystem.hpp"
#include "zsolve/Controller.hpp"
#include "zsolve/Norms.hpp"
#include "zsolve/Timer.h"

template <typename T> class Heuristics
{
    static int chooseNextVariableByZeros (Lattice <T> & lattice, BitSet & allowed)
    {
    // count zeros
    int* zeros = new int[lattice.variables ()];
    for (size_t i = 0; i < lattice.variables (); i++)
    {
        zeros[i] = 0;
        if (!allowed[i])
            continue;
        for (size_t j = 0; j < lattice.vectors (); j++)
        {
            if (lattice[j][i] == 0)
                zeros[i] ++;
        }
    }

    // choose one with most zeros
    int best_column = -1;
    for (size_t i = 0; i < lattice.variables (); i++)
    {
        if (!allowed[i])
            continue;
        if (best_column < 0 || (zeros[i] > zeros[best_column]))
            best_column = i;
        }
    delete[] zeros;

    return best_column;
}

static int chooseNextVariableByFile (Lattice <T> & lattice, BitSet & allowed)
{
    static int count = 0;
    static int* permutation = NULL;
    static std::string file = "permutation";

    if (permutation == NULL)
    {
        permutation = new int[lattice.variables ()];
        for (unsigned int i = 0; i < lattice.variables (); i++)
            permutation[i] = i;
    }

    int best_column = -1;

    std::ifstream stream (file.c_str ());
    for (int i = 0; i <= count; i++)
    {
        stream >> best_column;
    }

    if (count < (int) lattice.variables ())
    {
        for (unsigned int i = 0; i < lattice.variables (); i++)
        {
            if (permutation[i] == best_column)
            {
                best_column = i;
                break;
            }
        }
    }
    else
        best_column = -1;

    if (best_column >= 0)
    {
        int tmp = permutation[count];
        permutation[count] = permutation[best_column];
        permutation[best_column] = tmp;
        count++;

        if (!allowed[best_column])
        {
            int count_allowed = 0;
            for (unsigned int i = 0; i < lattice.variables (); i++)
                if (allowed[i])
                    count_allowed++;
            if (count_allowed > 0)
            {
                std::cerr << "Wanted to choose " << best_column << ", but it is not allowed!" << std::endl;
                exit (1);
            }
            else
                return -1;
        }
    }


    return best_column;
}

static int chooseNextVariableByRandom (Lattice <T> & lattice, BitSet & allowed)
{
    int count = 0;
    static bool inited = false;

    if (!inited)
    {
        srand (time (NULL));
        inited = true;
    }

    for (unsigned int i = 0; i < lattice.variables (); i++)
    {
        if (allowed[i])
            count++;
    }

    if (count == 0)
        return -1;

    int random = rand () % count;

    for (unsigned int i = 0; i < lattice.variables (); i++)
    {
        if (allowed[i])
        {
            if (random == 0)
            {
                return i;
            }
            random--;
        }
    }
    return -1;
}

public:
static int chooseNextVariable (Lattice <T> & lattice, BitSet & allowed)
{
    int best_column;

/*    std::cerr << "Allowed columns = [";
    for (unsigned int i = 0; i < lattice.variables (); i++)
        std::cerr << (allowed[i] ? "1 " : "0 ");
    std::cerr << std::endl;
*/
    best_column = chooseNextVariableByZeros (lattice, allowed);
    //best_column = chooseNextVariableByFile (lattice, allowed);
    //best_column = chooseNextVariableByRandom (lattice, allowed);
    
    //std::cerr << "CHOICE: taking " << best_column << " as next column" << std::endl;

    return best_column;
}

};

#endif
