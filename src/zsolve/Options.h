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

#ifndef _4ti2_zsolve__Options_
#define _4ti2_zsolve__Options_

#include <iostream>

namespace _4ti2_zsolve_
{

class Options
{
protected:
    std::string m_project;
    int m_verbosity;
    int m_loglevel;
    int m_backup_frequency;
    bool m_resume;
    bool m_hilbert;
    bool m_graver;
    bool m_maxnorm;
    int m_precision;

public:

    Options (int argc, char ** argv);

    static void print_banner ();
    void print_usage () const;
    void print_precision () const;

    std::string project () const;
    int verbosity () const;
    int loglevel () const;
    int backup_frequency () const;
    bool resume () const;
    bool hilbert () const;
    bool graver () const;
    bool maxnorm () const;
    int precision () const;

    friend std::istream& operator>>(std::istream& in, Options& options);
};

std::istream& operator>>(std::istream& in, Options& options);

} // namespace _4ti2_zsolve_
 
#endif
