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

#ifndef _4ti2_groebner__TermOrder_
#define _4ti2_groebner__TermOrder_

namespace _4ti2_
{

class TermOrder
{
public:
    TermOrder(int cost_start, int cost_end, int rs_end);
    static int direction(const TermOrder& term_order, const Binomial& b);
    static bool orientate(const TermOrder& term_order, Binomial& b);
protected:
    int cost_start;
    int cost_end;
    int rs_end;
};

inline
TermOrder::TermOrder(int _cost_start, int _cost_end, int _rs_end)
{
    cost_start = _cost_start;
    cost_end = _cost_end;
    rs_end = _rs_end;
}

inline
bool
TermOrder::orientate(const TermOrder& term_order, Binomial& b)
{
    Index i = term_order.cost_start;
    while (i < term_order.cost_end && b[i] == 0) { ++i; }
    if (i == term_order.cost_end)
    {
        i = 0;
        while (i < term_order.rs_end && b[i] == 0) { ++i; }
        if (i == term_order.rs_end) { return false; } // the binomial is zero.
        else if (b[i] > 0) { b.flip(); }
    }
    else if (b[i] < 0) { b.flip(); }
    return true;
}

inline
int
TermOrder::direction(const TermOrder& term_order, const Binomial& b)
{
    Index i = term_order.cost_start;
    while (i < term_order.cost_end && b[i] == 0) { ++i; }
    if (i == term_order.cost_end)
    {
        i = 0;
        while (i < term_order.rs_end && b[i] == 0) { ++i; }
        if (i == term_order.rs_end) { return 0; } // the binomial is zero.
        else if (b[i] > 0) { return -1; }
    }
    else if (b[i] < 0) { return -1; }
    return 1;
}

} // namespace _4ti2_

#endif

