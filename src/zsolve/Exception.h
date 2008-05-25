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

#ifndef _4ti2_zsolve__Exception_
#define _4ti2_zsolve__Exception_

// IO

namespace _4ti2_zsolve_
{

class PrecisionException
{
 protected:
    int precision_;
 public:
    PrecisionException (int prec)
    {
	precision_ = prec;
    }

    int precision ()
    {
	return precision_;
    }
};


class IOException
{
protected:
    std::string m_message;
    bool m_print;
public:
    IOException (std::string msg = "Unknown reason.", bool p = true)
    {
        m_print = p;
        m_message = msg;
    }

    bool print ()
    {
        return m_print;
    }

    friend std::ostream& operator<< (std::ostream& out, const IOException& exception);
};

std::ostream& operator<< (std::ostream& out, const IOException& exception)
{
    out << "Input error: " << exception.m_message << std::endl;
    return out;
}

} // namespace _4ti2_zsolve_

#endif
