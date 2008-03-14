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

#ifndef _4ti2__Norms_
#define _4ti2__Norms_

#include <map>

template <typename T> class NormPair
{
public:
    T first;
    T second;
    T sum;

    NormPair (T f, T s)
    {
        if (f > s)
        {
            second = f;
            first = s;
        }
        else
        {
            first = f;
            second = s;
        }
        sum = f + s;
    }

    NormPair (const NormPair& other)
    {
        first = other.first;
        second = other.second;
        sum = other.sum;
    }
};

template <typename T> inline bool operator< (const NormPair <T> & a, const NormPair <T> & b)
{
    return a.sum < b.sum || (a.sum == b.sum && a.first < b.first);
}

template <typename T> class Norms
{
protected:
    std::vector <T> norms;
    std::vector <NormPair <T> > pairs;

public:
    Norms ()
    {
        
    }

    void insert (T value)
    {
        for (typename std::vector <T> :: iterator iter = norms.begin (); iter != norms.end (); iter++)
        {
            if (*iter == value)
                return;
            if (*iter > value)
            {
                // insert before
                
            }
        }
    }

    ~Norms ()
    {
        norms.clear ();
        pairs.clear ();
    }
};

#endif
