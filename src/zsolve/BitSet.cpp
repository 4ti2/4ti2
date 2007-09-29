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

#include "BitSet.h"

BitSet::BitSet (const size_t size, bool value)
{
    m_size = size;
    m_blocks = needed_blocks (size);
    m_data = new BlockType [m_blocks];
    for (size_t i = 0; i < m_blocks; i++)
        m_data[i] = value ? ALL_ONES_BLOCK : ALL_ZEROS_BLOCK;
}

BitSet::~BitSet ()
{
    delete[] m_data;
}

size_t BitSet::needed_blocks (const size_t size) const
{
    return (size / BITS_PER_BLOCK) + (((size % BITS_PER_BLOCK) == 0) ? 0 : 1 );
}

BlockType BitSet::last_block_mask () const
{
    size_t rest = m_size % BITS_PER_BLOCK;
    if (rest == 0)
        return ALL_ZEROS_BLOCK;
    else if (rest+1 == BITS_PER_BLOCK)
        return ALL_ONES_BLOCK;
    else
        return (1 << (rest+1))-1;
}

size_t BitSet::size () const
{
    return m_size;
}

size_t BitSet::blocks () const
{
    return m_blocks;
}

bool BitSet::operator[] (const size_t index)
{
    return get (index);
}

void BitSet::set (size_t index)
{
    size_t block = index / BITS_PER_BLOCK;
    size_t pos = index % BITS_PER_BLOCK;
    m_data[block] |= (1 << pos);
}

void BitSet::unset (size_t index)
{
    size_t block = index / BITS_PER_BLOCK;
    size_t pos = index % BITS_PER_BLOCK;
    m_data[block] &= ~(1 << pos);
}

bool BitSet::get (size_t index) const
{
    size_t block = index / BITS_PER_BLOCK;
    size_t pos = index % BITS_PER_BLOCK;
    int result = m_data[block] & (1 << pos);
    return result != 0;
}

void BitSet::one ()
{
    for (size_t i = 0; i < m_blocks; i++)
        m_data[i] = ALL_ONES_BLOCK;
}

void BitSet::zero ()
{
    for (size_t i = 0; i < m_blocks; i++)
        m_data[i] = 0;
}

void BitSet::complement ()
{
    for (size_t i = 0; i < m_blocks; i++)
        m_data[i] = ~m_data[i];
}

bool BitSet::is_zero () const
{
    for (size_t i = 0; i < m_blocks-1; i++)
        if (m_data[i] != ALL_ZEROS_BLOCK)
            return false;
    return (m_data[m_blocks-1] & (last_block_mask ())) == ALL_ZEROS_BLOCK;
}

bool BitSet::is_one () const
{
    for (size_t i = 0; i < m_blocks-1; i++)
        if (m_data[i] != ALL_ONES_BLOCK)
            return false;
    return (m_data[m_blocks-1] | ~ (last_block_mask ())) == ALL_ONES_BLOCK;
}

void BitSet::set_intersection (const BitSet& other)
{
    assert (m_size == other.m_size);

    for (size_t i = 0; i < m_blocks; i++)
        m_data[i] &= other.m_data[i];
}

void BitSet::set_union (const BitSet& other)
{
    assert (m_size == other.m_size);

    for (size_t i = 0; i < m_blocks; i++)
        m_data[i] |= other.m_data[i];
}

