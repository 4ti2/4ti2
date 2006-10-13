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

#include "Markov.h"
#include "BinomialFactory.h"
#include "Generation.h"
#include "Bounded.h"
#include "Globals.h"

#include <iostream>
#include <iomanip>
#include "VectorArrayStream.h"
#include "FeasibleStream.h"
#include "BasicGeneration.h"

#define DEBUG_4ti2(X) //X

using namespace _4ti2_;

Markov::Markov(Generation* _gen)
        :gen(_gen)
{
    if (gen == 0) { gen = new BasicGeneration; }
}

Markov::~Markov()
{
}

void
Markov::compute(
                Feasible& feasible,
                VectorArray& vs)
{
    *out << "Computing Miminal Generation Set ...\n";
    t.reset(); // Reset the stopwatch.

    if (vs.get_number() != 0)
    {
        DEBUG_4ti2(
            *out << "Slow Algorithm\n";
            *out << feasible;
            *out << "VS:\n" << vs << "\n";
            *out << std::endl;
        )

        VectorArray cost(0, feasible.get_dimension());
        BinomialFactory factory(feasible, cost);
        WeightedBinomialSet c;
        factory.convert(vs, c);

        BinomialSet bs;
        algorithm(c, bs);
        factory.convert(bs, vs);
    }

    *out << "\r";
    *out << "  Size: " << std::setw(6) << vs.get_number();
    *out << ", Time: " << t << " / ";
    *out << Timer::global << " secs. Done." << std::endl;
}

void
Markov::compute(
                Feasible& feasible,
                const VectorArray& cost,
                VectorArray& vs)
{
    *out << "Computing Miminal Generation Set (Fast)...\n";
    t.reset(); // Reset the stopwatch.
    if (vs.get_number() != 0)
    {
        DEBUG_4ti2(
            *out << "Fast Algorithm\n";
            *out << "Cost:\n" << cost << "\n";
            *out << feasible;
            *out << "VS:\n" << vs << "\n";
            *out << std::endl;
        )

        BinomialFactory factory(feasible, cost);
        WeightedBinomialSet c;
        factory.convert(vs, c);

        factory.add_weight(feasible.get_grading(), c.max_grade());
        DEBUG_4ti2(*out<<"Max Grade: "<< c.max_grade()<<"\n";)

        BinomialSet bs;
        fast_algorithm(c, bs);
        factory.convert(bs, vs);
    }

    *out << "\r";
    *out << "  Size: " << std::setw(6) << vs.get_number();
    *out << ", Time: " << t << " / ";
    *out << Timer::global << " secs. Done." << std::endl;
}

bool
Markov::algorithm(
                WeightedBinomialSet& v,
                BinomialSet& gens)
{
    Binomial b;
    WeightedBinomialSet s;
    BinomialSet bs;
    Grade current_grade(v.min_grade());

    int num_iterations = 0;
    while(!s.empty() || !v.empty())
    {
        if (s.empty()) { current_grade = v.min_grade(); }
        else if (v.empty()) { current_grade = s.min_grade(); }
        else
        {
            if (s.min_grade() < v.min_grade())
            {
                current_grade = s.min_grade();
            }
            else
            {
                current_grade = v.min_grade();
            }
        }
        while (!s.empty() && s.min_grade() == current_grade)
        {
            ++num_iterations;
            s.next(b);
            bool zero = false;
            bs.reduce(b, zero);
            if (!zero)
            {
                bs.add(b);
                gen->generate(bs, bs.get_number()-1, s);
            }
            if (num_iterations % Globals::output_freq == 0)
            {
                *out << "\r";
                *out << "  Size: " << std::setw(6) << gens.get_number();
                *out << ", Grade: " << std::setw(6) << current_grade;
                *out << ", ToDo: " << std::setw(6) << s.get_size() << std::flush;
            }
        }
        while (!v.empty() && v.min_grade() == current_grade)
        {
            ++num_iterations;
            v.next(b);
            bool zero = false;
            bs.reduce(b, zero);
            if (!zero)
            {
                bs.add(b);
                gens.add(b);
                gen->generate(bs, bs.get_number()-1, s);
            }
            if (num_iterations % Globals::output_freq == 0)
            {
                *out << "\r";
                *out << "  Size: " << std::setw(6) << gens.get_number();
                *out << ", Grade: " << std::setw(6) << current_grade;
                *out << ", ToDo: " << std::setw(6) << s.get_size() << std::flush;
            }
        }
    }
    return true;
}


bool
Markov::fast_algorithm(
                WeightedBinomialSet& v,
                BinomialSet& gens)
{
    Binomial b;
    WeightedBinomialSet s;
    BinomialSet bs;
    Grade current_grade(v.min_grade());

    int num_iterations = 0;
    while(!s.empty() || !v.empty())
    {
        if (s.empty()) { current_grade = v.min_grade(); }
        else if (v.empty()) { current_grade = s.min_grade(); }
        else
        {
            if (s.min_grade() < v.min_grade())
            {
                current_grade = s.min_grade();
            }
            else
            {
                current_grade = v.min_grade();
            }
        }
        while (!s.empty() && s.min_grade() == current_grade)
        {
            ++num_iterations;
            s.next(b);
            bool zero = false;
            bs.reduce(b, zero);
            if (!zero)
            {
                bs.add(b);
                gen->generate(bs, bs.get_number()-1, s);
            }
            if (num_iterations % Globals::output_freq == 0)
            {
                *out << "\r";
                *out << "  Size: " << std::setw(6) << gens.get_number();
                *out << ", Grade: " << std::setw(6) << current_grade;
                *out << ", ToDo: " << std::setw(6) << s.get_size();
                *out << std::flush;
            }
        }
        while (!v.empty() && v.min_grade() == current_grade)
        {
            ++num_iterations;
            v.next(b);
            if (!bs.reducable(b))
            {
                bs.add(b);
                gens.add(b);
                gen->generate(bs, bs.get_number()-1, s);
            }
            if (num_iterations % Globals::output_freq == 0)
            {
                *out << "\r";
                *out << "  Size: " << std::setw(6) << gens.get_number();
                *out << ", Grade: " << std::setw(6) << current_grade;
                *out << ", ToDo: " << std::setw(6) << s.get_size();
                *out << std::flush;
            }
        }
    }
    return true;
}
