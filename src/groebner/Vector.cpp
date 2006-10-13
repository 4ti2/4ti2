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

#include "Vector.h"
#include "Euclidean.h"

using namespace _4ti2_;

Vector::Vector()
        : vector(0), size(0)
{
}

Vector::Vector(const Vector& v)
        : size(v.size)
{
    vector = new IntegerType[size];
    *this = v;
}

Vector::Vector(Size s)
        : size(s)
{
    assert(s >= 0);
    vector = new IntegerType[size];
}

Vector::Vector(Size s, IntegerType v)
        : size(s)
{
    assert(s >= 0);
    vector = new IntegerType[size];
    for (Index i = 0; i < size; ++i) { vector[i] = v; }
}

Vector::~Vector()
{
    delete [] vector;
}

void
Vector::permute(const Permutation& p)
{
    assert((Size) p.size() == size);
    Vector temp(*this);
    for (Index i = 0; i < size; ++i)
    {
        vector[i] = temp[p[i]];
    }
}

void
Vector::normalise()
{
    Index i = 0;
    while(i < size && vector[i] == 0) { ++i; }
    if (i == size) return;
    IntegerType gcd = vector[i];
    if (gcd == 1) return;
    ++i;
    while (i < size && vector[i] == 0) { ++i; }
    while (i < size)
    {
        euclidean(gcd, vector[i], gcd);
        if (gcd == 1) return;
        ++i;
        while (i < size && vector[i] == 0) { ++i; }
    }
    if (gcd != 1) { div(gcd); }
}
