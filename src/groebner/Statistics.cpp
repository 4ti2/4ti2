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

#include "Statistics.h"

#include <iomanip>

using namespace _4ti2_;

Statistics::Statistics()
{
#ifdef STATS_4TI2
    num_critical_pairs = 0;
    num_unmarked_pairs = 0;
    num_disjoint_critical_pairs = 0;
    num_graded_critical_pairs = 0;
    num_syzergy_critical_pairs = 0;
    num_non_duplicates = 0;
    num_reduction_steps = 0;
    num_reductions = 0;
    num_reducable_checks = 0;
    size_of_set = 0;
    size_of_set_before_minimal = 0;
#endif
}

void
Statistics::finalise()
{
}

void
Statistics::print(std::ostream& out)
{
#ifdef STATS_4TI2
    static const int WIDTH = 15;
    out << "Statistics for computing test set" << std::endl;
    out << "---------------------------------" << std::endl;
    out << "Size of test set               : "
            << std::setw(WIDTH) << size_of_set << std::endl;
    out << "Size of test set before minimal: "
            << std::setw(WIDTH) << size_of_set_before_minimal << std::endl;
    out << "Number of critical pairs       : "
            << std::setw(WIDTH) << num_critical_pairs << std::endl;
    out << "Number of unmarked pairs       : "
            << std::setw(WIDTH) << num_unmarked_pairs << std::endl;
    out << "Number of disjoint pairs       : "
            << std::setw(WIDTH) << num_disjoint_critical_pairs << std::endl;
    out << "Number of syzergy pairs        : "
            << std::setw(WIDTH) << num_syzergy_critical_pairs << std::endl;
    out << "Number of graded pairs         : "
            << std::setw(WIDTH) << num_graded_critical_pairs << std::endl;
    out << "Number of non duplicate pairs  : "
            << std::setw(WIDTH) << num_non_duplicates << std::endl;
    out << "Number of reductions           : "
            << std::setw(WIDTH) << num_reductions << std::endl;
    out << "Number of reduction steps      : "
            << std::setw(WIDTH) << num_reduction_steps << std::endl;
    out << "Number of reducable checks     : "
            << std::setw(WIDTH) << num_reducable_checks << std::endl;
    out << std::endl;
#endif
}
