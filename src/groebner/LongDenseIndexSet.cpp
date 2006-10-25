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

#include "LongDenseIndexSet.h"

using namespace _4ti2_;

const Size LongDenseIndexSet::sizeofstorage;
BlockType LongDenseIndexSet::set_masks[];
BlockType LongDenseIndexSet::unset_masks[];
BlockType LongDenseIndexSet::unused_masks[];
const unsigned char LongDenseIndexSet::bit_count[] = {
0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};

void
LongDenseIndexSet::initialise()
{
    static bool initialised = false;
    if (initialised == false)
    {
        BlockType mask = 1;
        for (Index i = 0; i < sizeofstorage; ++i)
        {
            set_masks[i] = mask;
            unset_masks[i] = ~mask;
            mask *= 2;
        }
        initialised = true;
        unused_masks[0] = 0;
        for (Index i = 1; i <= sizeofstorage; ++i)
        {
            unused_masks[i] = unused_masks[i-1] | set_masks[i-1];
        }
    }
}

void
LongDenseIndexSet::resize(Size s)
{
    Size n = get_num_blocks(s);
    if (n == num_blocks)
    {
        size = s;
        unset_unused_bits();
    }
    else if (n < num_blocks)
    {
        BlockType* bs = new BlockType[n];
        for (Index i = 0; i < n; ++i) { bs[i] = blocks[i]; }
        delete [] blocks;
        blocks = bs;
        size = s;
        unset_unused_bits();
    }
    else // n > num_blocks
    {
        BlockType* bs = new BlockType[n];
        for (Index i = 0; i < num_blocks; ++i) { bs[i] = blocks[i]; }
        for (Index i = num_blocks; i < n; ++i) { bs[i] = 0; }
        delete [] blocks;
        blocks = bs;
        size = s;
        unset_unused_bits();
    }
}
