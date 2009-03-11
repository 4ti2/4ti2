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

#ifndef _4ti2_zsolve__BitSet_
#define _4ti2_zsolve__BitSet_

#include <cassert>
#include <iostream>

namespace _4ti2_zsolve_
{

typedef uint32_t BlockType;
#define ALL_ZEROS_BLOCK 0
#define ALL_ONES_BLOCK UINT32_MAX
#define BITS_PER_BLOCK 32
#define BYTES_PER_BLOCK 4

class BitSet
{
private:
    BlockType *m_data;
    size_t m_size;
    size_t m_blocks;

    size_t needed_blocks (const size_t size) const;
    BlockType last_block_mask () const;

public:
    BitSet (const size_t size, bool value = false);
    ~BitSet ();

    size_t size () const;
    size_t blocks () const;
    bool operator[] (const size_t index);
    void set (size_t index);
    void unset (size_t index);
    bool get (size_t index) const;
    void one ();
    void zero ();
    void complement ();
    bool is_zero () const;
    bool is_one () const;
    void set_intersection (const BitSet& other);
    void set_union (const BitSet& other);
};

} // namespace _4ti2_zsolve_
 
#endif
