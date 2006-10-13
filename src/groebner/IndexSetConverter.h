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

#ifndef _4ti2__IndexSetConverter_
#define _4ti2__IndexSetConverter_

namespace _4ti2_
{

template <class IndexSet1, class IndexSet2>
void convert(const IndexSet1& i1, IndexSet2&i2);

} // namespace _4ti2_

template <class IndexSet1, class IndexSet2>
inline
void
_4ti2_::convert(const IndexSet1& i1, IndexSet2& i2)
{
    i2.zero();
    for (int i = 0; i < i1.get_size(); ++i) { if (i1[i]) { i2.set(i); } }
}

#endif
