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

#ifndef _4ti2__Statistics_
#define _4ti2__Statistics_

#define STATS_4TI2
#ifndef NDEBUG
#define STATS_4TI2
#endif
#include <iostream>

namespace _4ti2_
{

class Statistics
{
public:
    Statistics();

    void finalise();
    void print(std::ostream& out);

    void incr_num_critical_pairs();
    void incr_num_unmarked_pairs();
    void incr_num_disjoint_critical_pairs();
    void incr_num_graded_critical_pairs();
    void incr_num_syzygy_critical_pairs();
    void incr_num_non_duplicates();
    void incr_num_reduction_steps();
    void incr_num_reductions();
    void incr_num_reducable_checks();
    void set_size_of_set(long int size);
    void set_size_of_set_before_minimal(long int size);

protected:

#ifdef STATS_4TI2
    long int num_critical_pairs; // Number of criticals generated.
    long int num_unmarked_pairs;
    long int num_disjoint_critical_pairs;
    long int num_graded_critical_pairs;
    long int num_syzergy_critical_pairs;
    long int num_non_duplicates;
    long int num_reduction_steps;
    long int num_reductions;
    long int num_reducable_checks;
    long int size_of_set;
    long int size_of_set_before_minimal;
#endif
};

inline
void
Statistics::incr_num_critical_pairs()
{
#ifdef STATS_4TI2
    ++num_critical_pairs;
#endif
}

inline
void
Statistics::incr_num_unmarked_pairs()
{
#ifdef STATS_4TI2
    ++num_unmarked_pairs;
#endif
}

inline
void
Statistics::incr_num_disjoint_critical_pairs()
{
#ifdef STATS_4TI2
    ++num_disjoint_critical_pairs;
#endif
}

inline
void
Statistics::incr_num_graded_critical_pairs()
{
#ifdef STATS_4TI2
    ++num_graded_critical_pairs;
#endif
}

inline
void
Statistics::incr_num_syzygy_critical_pairs()
{
#ifdef STATS_4TI2
    ++num_syzergy_critical_pairs;
#endif
}

inline
void
Statistics::incr_num_non_duplicates()
{
#ifdef STATS_4TI2
    ++num_non_duplicates;
#endif
}

inline
void
Statistics::incr_num_reduction_steps()
{
#ifdef STATS_4TI2
    ++num_reduction_steps;
#endif
}

inline
void
Statistics::incr_num_reductions()
{
#ifdef STATS_4TI2
    ++num_reductions;
#endif
}

inline
void
Statistics::incr_num_reducable_checks()
{
#ifdef STATS_4TI2
    ++num_reducable_checks;
#endif
}

inline
void
Statistics::set_size_of_set(long int size)
{
#ifdef STATS_4TI2
    size_of_set = size;
#endif
}

inline
void
Statistics::set_size_of_set_before_minimal(long int size)
{
#ifdef STATS_4TI2
    size_of_set_before_minimal = size;
#endif
}

} // namespace _4ti2_

#endif
