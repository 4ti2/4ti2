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

#ifndef _4ti2_groebner__LongDenseIndexSet_
#define _4ti2_groebner__LongDenseIndexSet_

#include <cassert>
#include <climits>

#include "groebner/Size.h"
#include "groebner/Index.h"

#include <inttypes.h>

#include <iostream>

namespace _4ti2_
{
typedef uint64_t BlockType;
#define ALL_ONES_BLOCK UINT64_MAX
#define BITS_PER_BLOCK (CHAR_BIT*sizeof(BlockType))
#define BYTES_PER_BLOCK (sizeof(BlockType))

class LongDenseIndexSet
{
public:
    explicit LongDenseIndexSet(Size _size, bool v = false);
    LongDenseIndexSet(const LongDenseIndexSet&);
    LongDenseIndexSet& operator=(const LongDenseIndexSet&);
    ~LongDenseIndexSet();

    static bool set_subset(
                    const LongDenseIndexSet& b1,
                    const LongDenseIndexSet& b2);
    static bool set_subset(
                    const LongDenseIndexSet& b1,
                    const LongDenseIndexSet& b2,
                    const LongDenseIndexSet& mask);
    static bool set_disjoint(
                    const LongDenseIndexSet& b1,
                    const LongDenseIndexSet& b2);
    static bool set_disjoint(
                    const LongDenseIndexSet& b1,
                    const LongDenseIndexSet& b2,
                    const LongDenseIndexSet& mask);

    static void set_union(
                    const LongDenseIndexSet& b1,
                    const LongDenseIndexSet& b2,
                    LongDenseIndexSet& b3);
    static void set_intersection(
                    const LongDenseIndexSet& b1,
                    const LongDenseIndexSet& b2,
                    LongDenseIndexSet& b3);
    static void set_difference(
                    const LongDenseIndexSet& b1,
                    const LongDenseIndexSet& b2,
                    LongDenseIndexSet& b3);
    static void set_complement(
                    const LongDenseIndexSet& b1,
                    LongDenseIndexSet& b2);

    void set_union(const LongDenseIndexSet& b);
    void set_intersection(const LongDenseIndexSet& b);
    void set_difference(const LongDenseIndexSet& b);
    void set_complement();

    static void extend(const LongDenseIndexSet& b1, LongDenseIndexSet& b2);
    static void shrink(const LongDenseIndexSet& b1, LongDenseIndexSet& b2);

    bool power_of_2() const;
    bool less_than_equal(int s) const;
    bool empty() const;

    Size count() const;

    bool operator[](Index index) const;

    Size get_size() const;
    Size get_num_blocks() const;
    BlockType get_block(int i) const;
    void resize(Size s);

    void set(Index index);
    void unset(Index index);
    void flip(Index index);
    void zero();
    void one();

    static void swap(LongDenseIndexSet& b1, LongDenseIndexSet& b2);

    friend bool operator==(const LongDenseIndexSet& b1, const LongDenseIndexSet& b2);
    friend bool operator<(const LongDenseIndexSet& b1, const LongDenseIndexSet& b2);

    static bool special(const LongDenseIndexSet& b1, const LongDenseIndexSet& b2);

#if 0
    class Indice
    {
    public:
        Indice(int i) { block_num = i/sizeofstorage; bit_num = i%sizeofstorage; }
        friend bool operator!=(Indice i1, Indice i2)
        { return (i1.block_num != i2.block_num) || (i1.bit_num != i2.bit_num); }
    private:
        int block_num;
        int bit_num;
        friend class LongDenseIndexSet;
    };
    bool operator[](Indice i) const { return (blocks[i.block_num] & set_masks[i.bit_num]); }
#endif

protected:
    LongDenseIndexSet();

    static void logical_ior(const LongDenseIndexSet& b1,
                    const LongDenseIndexSet& b2,
                    LongDenseIndexSet& b3);
    static void logical_xor(const LongDenseIndexSet& b1,
                    const LongDenseIndexSet& b2,
                    LongDenseIndexSet& b3);
    static void logical_and(const LongDenseIndexSet& b1,
                    const LongDenseIndexSet& b2,
                    LongDenseIndexSet& b3);
    static void logical_not(const LongDenseIndexSet& b1,
                    LongDenseIndexSet& b2);

    void logical_ior(const LongDenseIndexSet& b);
    void logical_xor(const LongDenseIndexSet& b);
    void logical_and(const LongDenseIndexSet& b);
    void logical_not();

    static Size get_num_blocks(Size _size);

    static BlockType set_masks[BITS_PER_BLOCK];
    static BlockType unset_masks[BITS_PER_BLOCK];
    static BlockType unused_masks[BITS_PER_BLOCK+1];
    static const Size sizeofstorage = BITS_PER_BLOCK;
    static const unsigned char bit_count[256];

    static void initialise();
    void unset_unused_bits();

    BlockType *blocks;
    Size size;
    Size num_blocks;
};

inline
bool
LongDenseIndexSet::power_of_2() const
{
#if 1
    bool c = false;
    for (Index i = 0; i < num_blocks; ++i)
    {
        BlockType tmp = blocks[i];
        while (tmp)
        {
            if (c) { return false; }
            tmp &= tmp-1;
            c = true;
        }
    }
    return true;
#endif

#if 0
    // Find the first non-empty block.
    Index i = 0;
    while (i < num_blocks)
    {
        if (blocks[i]) { break; }
        ++i;
    }
    if (i == num_blocks) { return true; } // There are no non-empty blocks.

    // Check whether there is only one bit set in the block.
    BlockType tmp = blocks[i];
    tmp &= tmp-1;
    if (tmp) { return false; }

    // Check whether there are no other non-empty blocks.
    ++i;
    while (i < num_blocks)
    {
        if (blocks[i]) { return false; }
        ++i;
    }
    return true;
#endif
}

inline
bool
LongDenseIndexSet::less_than_equal(int s) const
{
    int c = 0;
    for (Index i = 0; i < num_blocks; ++i)
    {
        BlockType tmp = blocks[i];
        while (tmp)
        {
            if (c == s) { return false; }
            tmp &= tmp-1;
            ++c;
        }
    }
    return true;
}

inline
bool
LongDenseIndexSet::empty() const
{
    for (Index i = 0; i < num_blocks; ++i)
    {
        if (blocks[i]) { return false; }
    }
    return true;
}

inline
void
LongDenseIndexSet::set(Index index)
{
    assert(index >= 0 && index < size);
    blocks[index / sizeofstorage] |= set_masks[index % sizeofstorage];
}

inline
void
LongDenseIndexSet::unset(Index index)
{
    assert(index >= 0 && index < size);
    blocks[index / sizeofstorage] &= unset_masks[index % sizeofstorage];
}

inline
void
LongDenseIndexSet::flip(Index index)
{
    assert(index >= 0 && index < size);
    blocks[index / sizeofstorage] ^= set_masks[index % sizeofstorage];
}

inline
void
LongDenseIndexSet::zero()
{
    for (Index i = 0; i < num_blocks; ++i) { blocks[i] = 0; }
}

inline
void
LongDenseIndexSet::one()
{
    for (Index i = 0; i < num_blocks; ++i) { blocks[i] = ALL_ONES_BLOCK; }
    unset_unused_bits();
}

inline
void
LongDenseIndexSet::swap(LongDenseIndexSet& b1, LongDenseIndexSet& b2)
{
    assert(b1.size == b2.size);
    BlockType* temp = b1.blocks;
    b1.blocks = b2.blocks;
    b2.blocks = temp;
}

#if 1
inline
bool
LongDenseIndexSet::operator[](Index index) const
{
    assert(index >= 0 && index <= size);
    // TODO: do we need the comparison with 0?
    return (set_masks[index%sizeofstorage] & blocks[index/sizeofstorage]) != 0;
}
#endif

inline
Size
LongDenseIndexSet::get_size() const
{
    return size;
}

inline
Size
LongDenseIndexSet::get_num_blocks() const
{
    return num_blocks;
}

inline
BlockType
LongDenseIndexSet::get_block(int i) const
{
    return blocks[i];
}

inline
bool
operator<(const LongDenseIndexSet& b1, const LongDenseIndexSet& b2)
{
    assert(b1.size == b2.size && b1.num_blocks == b2.num_blocks);
    Index i = 0;
    while (i < b1.num_blocks && b1.blocks[i] == b2.blocks[i]) ++i;
    if (i < b1.num_blocks && b1.blocks[i] < b2.blocks[i]) return true;
    return false;
}

inline
bool
operator==(const LongDenseIndexSet& b1, const LongDenseIndexSet& b2)
{
    assert(b1.size == b2.size);
    assert(b1.num_blocks == b2.num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i)
    {
        if (b1.blocks[i] != b2.blocks[i]) return false;
    }
    return true;
}

inline
void
LongDenseIndexSet::logical_ior(const LongDenseIndexSet& b1,
                 const LongDenseIndexSet& b2,
                 LongDenseIndexSet& b3)
{
    assert(b1.num_blocks == b2.num_blocks && b1.num_blocks == b3.num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i)
    {
        b3.blocks[i] = b1.blocks[i] | b2.blocks[i];
    }
}

inline
void
LongDenseIndexSet::logical_xor(const LongDenseIndexSet& b1,
                 const LongDenseIndexSet& b2,
                 LongDenseIndexSet& b3)
{
    assert(b1.num_blocks == b2.num_blocks && b1.num_blocks == b3.num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i)
    {
        b3.blocks[i] = b1.blocks[i] ^ b2.blocks[i];
    }
}

inline
void
LongDenseIndexSet::logical_and(const LongDenseIndexSet& b1,
                 const LongDenseIndexSet& b2,
                 LongDenseIndexSet& b3)
{
    assert(b1.num_blocks == b2.num_blocks && b1.num_blocks == b3.num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i)
    {
        b3.blocks[i] = b1.blocks[i] & b2.blocks[i];
    }
}

inline
void
LongDenseIndexSet::logical_not(const LongDenseIndexSet& b1,
                 LongDenseIndexSet& b2)
{
    assert(b1.num_blocks == b2.num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i)
    {
        b2.blocks[i] = ~b1.blocks[i];
    }
    b2.unset_unused_bits();
}

inline
void
LongDenseIndexSet::logical_ior(const LongDenseIndexSet& b)
{
    assert(num_blocks == b.num_blocks);
    for (Index i = 0; i < num_blocks; ++i)
    {
        blocks[i] |= b.blocks[i];
    }
}

inline
void
LongDenseIndexSet::logical_xor(const LongDenseIndexSet& b)
{
    assert(num_blocks == b.num_blocks);
    for (Index i = 0; i < num_blocks; ++i)
    {
        blocks[i] ^= b.blocks[i];
    }
}

inline
void
LongDenseIndexSet::logical_and(const LongDenseIndexSet& b)
{
    assert(num_blocks == b.num_blocks);
    for (Index i = 0; i < num_blocks; ++i)
    {
        blocks[i] &= b.blocks[i];
    }
}

inline
void
LongDenseIndexSet::logical_not()
{
    for (Index i = 0; i < num_blocks; ++i)
    {
        blocks[i] = ~blocks[i];
    }
    unset_unused_bits();
}

inline
bool
LongDenseIndexSet::set_subset(const LongDenseIndexSet& b1, const LongDenseIndexSet& b2)
{
    assert(b1.num_blocks == b2.num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i)
    {
        if (b1.blocks[i] != (b1.blocks[i] & b2.blocks[i])) return false;
    }
    return true;
}

inline
bool
LongDenseIndexSet::set_subset(const LongDenseIndexSet& b1, const LongDenseIndexSet& b2, const LongDenseIndexSet& mask)
{
    assert(b1.num_blocks == b2.num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i)
    {
        if ((b1.blocks[i] & mask.blocks[i]) != (b1.blocks[i] & b2.blocks[i] & mask.blocks[i]))
        {
            return false;
        }
    }
    return true;
}

inline
bool
LongDenseIndexSet::set_disjoint(const LongDenseIndexSet& b1, const LongDenseIndexSet& b2)
{
    assert(b1.num_blocks == b2.num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i)
    {
        if ((b1.blocks[i] & b2.blocks[i]) != 0) return false;
    }
    return true;
}

inline
bool
LongDenseIndexSet::set_disjoint(const LongDenseIndexSet& b1, const LongDenseIndexSet& b2, const LongDenseIndexSet& mask)
{
    assert(b1.num_blocks == b2.num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i)
    {
        if ((b1.blocks[i] & b2.blocks[i] & mask.blocks[i]) != 0) return false;
    }
    return true;
}

inline
void
LongDenseIndexSet::set_union(
                const LongDenseIndexSet& b1,
                const LongDenseIndexSet& b2,
                LongDenseIndexSet& b3)
{
    LongDenseIndexSet::logical_ior(b1, b2, b3);
}

inline
void
LongDenseIndexSet::set_intersection(
                const LongDenseIndexSet& b1,
                const LongDenseIndexSet& b2,
                LongDenseIndexSet& b3)
{
    LongDenseIndexSet::logical_and(b1, b2, b3);
}

inline
void
LongDenseIndexSet::set_difference(
                const LongDenseIndexSet& b1,
                const LongDenseIndexSet& b2,
                LongDenseIndexSet& b3)
{
    assert(b1.num_blocks == b2.num_blocks && b1.num_blocks == b3.num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i)
    {
        b3.blocks[i] = b1.blocks[i] & (~b2.blocks[i]);
    }
}

inline
void
LongDenseIndexSet::set_complement(
                const LongDenseIndexSet& b1,
                LongDenseIndexSet& b2)
{
    LongDenseIndexSet::logical_not(b1, b2);
}

inline
void
LongDenseIndexSet::set_union(const LongDenseIndexSet& b)
{
    logical_ior(b);
}

inline
void
LongDenseIndexSet::set_intersection(const LongDenseIndexSet& b)
{
    logical_and(b);
}

inline
void
LongDenseIndexSet::set_difference(const LongDenseIndexSet& b)
{
    assert(num_blocks == b.num_blocks);
    for (Index i = 0; i < num_blocks; ++i)
    {
        blocks[i] &= ~b.blocks[i];
    }
}

inline
void
LongDenseIndexSet::set_complement()
{
    logical_not();
}

inline
void
LongDenseIndexSet::extend(const LongDenseIndexSet& b1, LongDenseIndexSet& b2)
{
    assert(b1.size <= b2.size);
    for (int i = 0; i < b1.num_blocks; ++i)
    {
        b2.blocks[i] = b1.blocks[i];
    }
    for (int i = b1.num_blocks; i < b2.num_blocks; ++i)
    {
        b2.blocks[i] = 0;
    }
}

inline
void
LongDenseIndexSet::shrink(const LongDenseIndexSet& b1, LongDenseIndexSet& b2)
{
    assert(b1.size >= b2.size);
    for (int i = 0; i < b2.num_blocks; ++i)
    {
        b2.blocks[i] = b1.blocks[i];
    }
}

inline
Size
LongDenseIndexSet::count() const
{
    // The following section of code only works for 64 bit block sizes which we can
    // assume (hopefully) because the BlockType is uint64_t.
    // The following code is based upon code obtained from the website
    // http://graphics.stanford.edu/~seander/bithacks.html
    BlockType* block = blocks;
    BlockType* end = blocks + num_blocks;
    Size c = 0;
    while (block != end)
    {
        BlockType const w = *block - ((*block >> 1) & (BlockType) 0x5555555555555555ULL);  // temp
        BlockType const x = (w & (BlockType) 0x3333333333333333ULL) +
            ((w >> 2) & (BlockType) 0x3333333333333333ULL);     // temp
        c += (((x + (x >> 4)) & (BlockType) 0x0F0F0F0F0F0F0F0FULL) * (BlockType) 0x0101010101010101ULL)
                >> 56;
        ++block;
    }
    return c;
#if 0
    // The following code is slower than the above code but it has the advantage
    // that it is independent of the block size. It is not used at present.
    const unsigned char* byte = (const unsigned char*) blocks;
    const unsigned char* end = byte + (num_blocks*BYTES_PER_BLOCK);
    Size c = 0;
    while (byte < end) 
    {
        c += bit_count[*byte];
        ++byte;
    }
    return c;
#endif
}

// TODO: There is a big-endian/little-endian problem here.
inline
bool
LongDenseIndexSet::special(const LongDenseIndexSet& b1, const LongDenseIndexSet& b2)
{
    BlockType* b1ptr = b1.blocks;
    BlockType* b2ptr = b2.blocks;
    BlockType* b1end = b1ptr + b1.num_blocks;
    int count = 0;
    BlockType tmp;
    const unsigned char* byte = (const unsigned char*) &tmp;
    while (b1ptr < b1end)
    {
        tmp = *b1ptr & *b2ptr;
        count += bit_count[*byte];
        count += bit_count[*(byte+1)];
        count += bit_count[*(byte+2)];
        count += bit_count[*(byte+3)];
        ++b1ptr;
        ++b2ptr;
    }
    return count < 2;
}

inline
void
LongDenseIndexSet::unset_unused_bits()
{
    if (size > 0)
    {
        Size unused_index = ((size-1) % sizeofstorage) + 1;
        blocks[num_blocks-1] &= unused_masks[unused_index];
    }
}

// TODO: Change to return (_size-1)/sizeofstorage + 1; Assuming the _size > 0.
inline
Size
LongDenseIndexSet::get_num_blocks(Size _size)
{
    if (_size % sizeofstorage == 0) { return _size/sizeofstorage; }
    else { return _size/sizeofstorage+1; }
}

inline
LongDenseIndexSet::LongDenseIndexSet()
        : blocks(0), size(0), num_blocks(0)
{
    //*out << "LongDenseIndexSet()\n";
}

inline
LongDenseIndexSet::LongDenseIndexSet(Size _size, bool v)
        : size(_size), num_blocks(get_num_blocks(_size))
{
    initialise();
    assert(_size >= 0);
    blocks = new BlockType[num_blocks];

    if (v == false) { zero(); }
    else { one(); } // v == true
}

inline
LongDenseIndexSet::LongDenseIndexSet(const LongDenseIndexSet& b)
        : size(b.size),
          num_blocks(b.num_blocks)
{
    blocks = new BlockType[num_blocks];
    *this = b;
}

inline
LongDenseIndexSet&
LongDenseIndexSet::operator=(const LongDenseIndexSet& b)
{
    assert(size == b.size);
    assert(num_blocks == b.num_blocks);
    for (Size i = 0; i < num_blocks; ++i) { blocks[i] = b.blocks[i]; }
    return *this;
}

inline
LongDenseIndexSet::~LongDenseIndexSet()
{
    delete [] blocks;
}


} // namespace _4ti2_

#undef ALL_ONES_BLOCK
#undef BITS_PER_BLOCK
#undef BYTES_PER_BLOCK

#endif
