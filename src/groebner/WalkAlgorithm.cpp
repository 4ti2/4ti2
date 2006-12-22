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

#include "WalkAlgorithm.h"
#include "Completion.h"
#include "VectorStream.h"
#include "VectorArrayStream.h"
#include "BinomialSet.h"
#include "BinomialSetStream.h"
#include "FlipCompletion.h"
#include "BinomialFactory.h"
#include "Globals.h"

#include <iostream>
#include <iomanip>

//#define DEBUG_4ti2(X) X
#include "Debug.h"

using namespace _4ti2_;

WalkAlgorithm::WalkAlgorithm()
{
    costnew_start = 0;
    costnew_end = 0;
    costold_start = 0;
    costold_end = 0;
}

WalkAlgorithm::~WalkAlgorithm()
{
}

void
WalkAlgorithm::compute(
                Feasible& feasible,
                const VectorArray& costold,
                VectorArray& vs,
                const VectorArray& costnew)
{
    DEBUG_4ti2(*out << "Cost1:\n" << costold << "\n";)
    DEBUG_4ti2(*out << "Cost2:\n" << costnew << "\n";)
    DEBUG_4ti2(*out << "VS:\n" << vs << "\n";)
    assert(costold.get_number() > 0 && costnew.get_number() > 0);
    t.reset();

    // TODO: We might need to add cost vectors to the old and new costs.
    VectorArray costall(costnew);
    costall.insert(costold);

    BinomialFactory factory(feasible, costall);

    costnew_start = Binomial::cost_start;
    costnew_end = Binomial::cost_start+costnew.get_number();
    costold_start = costnew_end;
    costold_end = Binomial::cost_end;

    BinomialSet bs;
    factory.convert(vs, bs, false);

    TermOrder term_order(costnew_start, costnew_end, Binomial::rs_end);

    // TODO: Perhaps the following value should be an command line option.
    const int reduction_freq = 200;

    int iteration = 0;
    DEBUG_4ti2(*out << "\nBefore BINOMIALS:\n" << bs << "\n";)
    int n;
    Binomial b;
    FlipCompletion alg;
    while(!next(bs, term_order, n))
    {
        if (iteration % Globals::output_freq == 0) 
        {
            *out << "\r";
            *out << std::setiosflags(std::ios_base::right);
            *out << "Iteration = " << std::setw(6) << iteration;
            *out << " Size = " << std::setw(6) << bs.get_number();
            *out << " tvalue " << std::setw(6) << std::setprecision(4);
            *out << std::resetiosflags(std::ios_base::right);
            *out << std::setiosflags(std::ios_base::left);
            *out << tvalue(bs[n]) << std::flush;
            *out << std::resetiosflags(std::ios_base::left);
        }
        DEBUG_4ti2(*out << "Iteration " << iteration << "\n";)
        DEBUG_4ti2(*out << "\nBefore BINOMIALS:\n" << bs << "\n";)
        b = bs[n];
        bs.remove(n);
        if (!bs.reducable(b))
        {
            b.flip();
            alg.algorithm(bs, b);
            bs.add(b);
            if (iteration % reduction_freq == 0)
            { bs.minimal(); bs.reduced(); } 
            DEBUG_4ti2(*out << "Size: " << std::setw(6) << bs.get_number() << "\n";)
            DEBUG_4ti2(*out << "After BINOMIALS:\n" << bs << "\n";)
            ++iteration;
        }
    }
    bs.minimal();
    bs.reduced();
    factory.convert(bs, vs);
    vs.sort();
    bs.clear();

    DEBUG_4ti2(*out << "VS:\n" << vs << "\n";)
    *out << "\r" << Globals::context;
    *out << "Iteration = " << std::setw(6) << iteration;
    *out << " Size: " << std::setw(6) << vs.get_number();
    *out << ", Time: " << t << " / ";
    *out << Timer::global << " secs. Done." << std::endl;
}

IntegerType
WalkAlgorithm::compare(const Binomial& b1, const Binomial& b2)
{
    IntegerType result;
    for (int i = costnew_start; i < costnew_end; ++i)
    {
        for (int j = costold_start; j < costold_end; ++j)
        {
            result = b1[j]*b2[i]-b1[i]*b2[j];
            if (result != 0) { return result; }
        }
        for (int j = 0; j < Binomial::rs_end; j++)
        {
            result = -b1[j]*b2[i]+b1[i]*b2[j];
            if (result != 0) { return result; }
        }
    }
    for (int i = 0; i < Binomial::rs_end; i++)
    {
        for (int j = costold_start; j < costold_end; ++j)
        {
            result = -b1[j]*b2[i]+b1[i]*b2[j];
            if (result != 0) { return result; }
        }
        for (int j = 0; j < Binomial::rs_end; j++)
        {
            result = b1[j]*b2[i]-b1[i]*b2[j];
            if (result != 0) { return result; }
        }
    }
    std::cerr << "Software Error: unexpected execution.\n";
    exit(1);
    return 0;
}

RationalType
WalkAlgorithm::tvalue(const Binomial& b)
{
    return (RationalType)b[costold_start]/(b[costold_start]-b[costnew_start]);
}

void
WalkAlgorithm::tvector(Vector& c1, Vector& c2, Vector& v, Vector& tv)
{
    Vector::sub(c2, Vector::dot(c1,v), c1, Vector::dot(c2,v), tv);
}

bool
WalkAlgorithm::next(const BinomialSet& bs, const TermOrder& term_order, int& min)
{
    min = 0;
    while (min < bs.get_number())
    {
        if (TermOrder::direction(term_order, bs[min]) < 0) { break; }
        else { ++min; }
    }
    if (min == bs.get_number()) { return true; }
    DEBUG_4ti2(*out << min << " tvalue = " << tvalue(bs[min]) << "\n";)
    for (int i = min+1; i < bs.get_number(); ++i)
    {
        if (TermOrder::direction(term_order, bs[i]) < 0)
        {
            DEBUG_4ti2(*out << i << " tvalue = " << tvalue(bs[i]) << "\n";)
            if (compare(bs[min], bs[i]) < 0) { min = i; }
        }
    }
    DEBUG_4ti2(*out << min << " min tvalue = " << tvalue(bs[min]) << "\n";)
    return false;
}
