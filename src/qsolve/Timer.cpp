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

#include "Timer.h"

#include <unistd.h>
#include <ctime>
#ifdef _POSIX_TIMERS
#include <sys/times.h>
#endif
#include <iomanip>

using namespace _4ti2_;

Timer Timer::global;

Timer::Timer()
{
    reset();
}

void
Timer::reset()
{
    start_time = get_time();
}

double
Timer::get_elapsed_time() const
{
    return  get_time() - start_time;
}

double
Timer::get_time()
{
#if defined _POSIX_TIMERS
    struct tms t;
    times(&t); 
    return (double) t.tms_utime/sysconf(_SC_CLK_TCK);
#else
    // clock() is the most widely supported function. Unfortunately, its value
    // can wrap around, so we try to use the above functions instead.
    return (double) clock() / CLOCKS_PER_SEC;
#endif
}

std::ostream& _4ti2_::operator<<(std::ostream& out, const Timer& t)
{
    out.precision(2);
    out.flags(std::ios::fixed);
    out.width(5);
    out << t.get_elapsed_time();

    return out;
}
