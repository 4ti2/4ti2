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

#ifndef _4ti2_groebner__ShortDenseIndexSet_
#define _4ti2_groebner__ShortDenseIndexSet_

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

class ShortDenseIndexSet
{
public:
    explicit ShortDenseIndexSet(Size _size, bool v = false);
    ShortDenseIndexSet(const ShortDenseIndexSet&);
    ShortDenseIndexSet& operator=(const ShortDenseIndexSet&);
    ~ShortDenseIndexSet();

    static const int max_size = BITS_PER_BLOCK;

    static bool set_subset(
                    const ShortDenseIndexSet& b1,
                    const ShortDenseIndexSet& b2);
    static bool set_subset(
                    const ShortDenseIndexSet& b1,
                    const ShortDenseIndexSet& b2,
                    const ShortDenseIndexSet& mask);
    static bool set_disjoint(
                    const ShortDenseIndexSet& b1,
                    const ShortDenseIndexSet& b2);
    static bool set_disjoint(
                    const ShortDenseIndexSet& b1,
                    const ShortDenseIndexSet& b2,
                    const ShortDenseIndexSet& mask);

    static void set_union(
                    const ShortDenseIndexSet& b1,
                    const ShortDenseIndexSet& b2,
                    ShortDenseIndexSet& b3);
    static void set_intersection(
                    const ShortDenseIndexSet& b1,
                    const ShortDenseIndexSet& b2,
                    ShortDenseIndexSet& b3);
    static void set_difference(
                    const ShortDenseIndexSet& b1,
                    const ShortDenseIndexSet& b2,
                    ShortDenseIndexSet& b3);
    static void set_complement(
                    const ShortDenseIndexSet& b1,
                    ShortDenseIndexSet& b2);

    void set_union(const ShortDenseIndexSet& b);
    void set_intersection(const ShortDenseIndexSet& b);
    void set_difference(const ShortDenseIndexSet& b);
    void set_complement();

    static void extend(const ShortDenseIndexSet& b1, ShortDenseIndexSet& b2);
    static void shrink(const ShortDenseIndexSet& b1, ShortDenseIndexSet& b2);

    bool power_of_2() const;
    bool less_than_equal(int s) const;
    bool empty() const;

    Size count() const;

    bool operator[](Index index) const;

    void set(Index index);
    void unset(Index index);
    void flip(Index index);
    void zero();
    void one();

    Size get_size() const;
    void resize(Size s);

    static void swap(ShortDenseIndexSet& b1, ShortDenseIndexSet& b2);

    friend bool operator==(const ShortDenseIndexSet& b1, const ShortDenseIndexSet& b2);
    friend bool operator<(const ShortDenseIndexSet& b1, const ShortDenseIndexSet& b2);

#if 0
    class Indice
    {
    public:
        Indice(int i) { b = set_masks[i]; }
        friend bool operator!=(Indice i1, Indice i2) { return i1.b != i2.b; }
    private:
        BlockType b;
        friend class ShortDenseIndexSet;
    };
    bool operator[](Indice i) const { return (block & i.b); }
#endif

protected:
    ShortDenseIndexSet();

    static void logical_ior(const ShortDenseIndexSet& b1,
                    const ShortDenseIndexSet& b2,
                    ShortDenseIndexSet& b3);
    static void logical_xor(const ShortDenseIndexSet& b1,
                    const ShortDenseIndexSet& b2,
                    ShortDenseIndexSet& b3);
    static void logical_and(const ShortDenseIndexSet& b1,
                    const ShortDenseIndexSet& b2,
                    ShortDenseIndexSet& b3);
    static void logical_not(const ShortDenseIndexSet& b1,
                    ShortDenseIndexSet& b2);

    void logical_ior(const ShortDenseIndexSet& b);
    void logical_xor(const ShortDenseIndexSet& b);
    void logical_and(const ShortDenseIndexSet& b);
    void logical_not();

    static BlockType set_masks[BITS_PER_BLOCK];
    static BlockType unset_masks[BITS_PER_BLOCK];
    static BlockType unused_masks[BITS_PER_BLOCK+1];
    static const Size sizeofstorage = BITS_PER_BLOCK;
    static const unsigned char bit_count[256];

    static void initialise();
    void unset_unused_bits();

    BlockType block;
    Size size;
};

inline
bool
ShortDenseIndexSet::power_of_2() const
{
    BlockType tmp = block;
    tmp &= tmp-1;
    return !(tmp);
}

inline
bool
ShortDenseIndexSet::less_than_equal(int s) const
{
    BlockType tmp = block;
    int c = 0;
    while (tmp)
    {
        if (c == s) { return false; }
        tmp &= tmp-1;
        c++;
    }
    return true;
}

inline
bool
ShortDenseIndexSet::empty() const
{
    return !(block);
}

inline
void
ShortDenseIndexSet::set(Index index)
{
    assert(index >= 0 && index < size);
    block |= set_masks[index];
}

inline
void
ShortDenseIndexSet::unset(Index index)
{
    assert(index >= 0 && index < size);
    block &= unset_masks[index];
}

inline
void
ShortDenseIndexSet::flip(Index index)
{
    assert(index >= 0 && index < size);
    block ^= set_masks[index];
}

inline
void
ShortDenseIndexSet::zero()
{
    block = 0;
}

inline
void
ShortDenseIndexSet::one()
{
    block = ALL_ONES_BLOCK;
    unset_unused_bits();
}

inline
void
ShortDenseIndexSet::swap(ShortDenseIndexSet& b1, ShortDenseIndexSet& b2)
{
    assert(b1.size == b2.size);
    BlockType temp = b1.block;
    b1.block = b2.block;
    b2.block = temp;
}

#if 1
inline
bool
ShortDenseIndexSet::operator[](Index index) const
{
    assert(index >= 0 && index <= size);
    // TODO: do we need the comparison with 0?
    return (set_masks[index] & block) != 0;
}
#endif

inline
Size
ShortDenseIndexSet::get_size() const
{
    return size;
}

inline
bool
operator<(const ShortDenseIndexSet& b1, const ShortDenseIndexSet& b2)
{
    assert(b1.size == b2.size);
    return (b1.block < b2.block);
}

inline
bool
operator==(const ShortDenseIndexSet& b1, const ShortDenseIndexSet& b2)
{
    assert(b1.size == b2.size);
    return (b1.block == b2.block);
}

inline
void
ShortDenseIndexSet::logical_ior(const ShortDenseIndexSet& b1,
                 const ShortDenseIndexSet& b2,
                 ShortDenseIndexSet& b3)
{
    b3.block = b1.block | b2.block;
}

inline
void
ShortDenseIndexSet::logical_xor(const ShortDenseIndexSet& b1,
                 const ShortDenseIndexSet& b2,
                 ShortDenseIndexSet& b3)
{
    b3.block = b1.block ^ b2.block;
}

inline
void
ShortDenseIndexSet::logical_and(const ShortDenseIndexSet& b1,
                 const ShortDenseIndexSet& b2,
                 ShortDenseIndexSet& b3)
{
    b3.block = b1.block & b2.block;
}

inline
void
ShortDenseIndexSet::logical_not(const ShortDenseIndexSet& b1,
                 ShortDenseIndexSet& b2)
{
    b2.block = ~b1.block;
    b2.unset_unused_bits();
}

inline
void
ShortDenseIndexSet::logical_ior(const ShortDenseIndexSet& b)
{
    block |= b.block;
}

inline
void
ShortDenseIndexSet::logical_xor(const ShortDenseIndexSet& b)
{
    block ^= b.block;
}

inline
void
ShortDenseIndexSet::logical_and(const ShortDenseIndexSet& b)
{
    block &= b.block;
}

inline
void
ShortDenseIndexSet::logical_not()
{
    block = ~block;
    unset_unused_bits();
}

inline
bool
ShortDenseIndexSet::set_subset(const ShortDenseIndexSet& b1, const ShortDenseIndexSet& b2)
{
    assert(b1.size == b2.size);
    return (b1.block == (b1.block & b2.block));
}

inline
bool
ShortDenseIndexSet::set_subset(const ShortDenseIndexSet& b1, const ShortDenseIndexSet& b2, const ShortDenseIndexSet& mask)
{
    assert(b1.size == b2.size && b1.size == mask.size);
    return ((b1.block & mask.block) == (b1.block & b2.block & mask.block));
}

inline
bool
ShortDenseIndexSet::set_disjoint(const ShortDenseIndexSet& b1, const ShortDenseIndexSet& b2)
{
    assert(b1.size == b2.size);
    return ((b1.block & b2.block) == 0);
}

inline
bool
ShortDenseIndexSet::set_disjoint(const ShortDenseIndexSet& b1, const ShortDenseIndexSet& b2, const ShortDenseIndexSet& mask)
{
    assert(b1.size == b2.size && b1.size == mask.size);
    return ((b1.block & b2.block & mask.block) == 0);
}

inline
void
ShortDenseIndexSet::set_union(
                const ShortDenseIndexSet& b1,
                const ShortDenseIndexSet& b2,
                ShortDenseIndexSet& b3)
{
    ShortDenseIndexSet::logical_ior(b1, b2, b3);
}

inline
void
ShortDenseIndexSet::set_intersection(
                const ShortDenseIndexSet& b1,
                const ShortDenseIndexSet& b2,
                ShortDenseIndexSet& b3)
{
    ShortDenseIndexSet::logical_and(b1, b2, b3);
}

inline
void
ShortDenseIndexSet::set_difference(
                const ShortDenseIndexSet& b1,
                const ShortDenseIndexSet& b2,
                ShortDenseIndexSet& b3)
{
    b3.block = b1.block & (~b2.block);
}

inline
void
ShortDenseIndexSet::set_complement(
                const ShortDenseIndexSet& b1,
                ShortDenseIndexSet& b2)
{
    ShortDenseIndexSet::logical_not(b1, b2);
}

inline
void
ShortDenseIndexSet::set_union(const ShortDenseIndexSet& b)
{
    logical_ior(b);
}

inline
void
ShortDenseIndexSet::set_intersection(const ShortDenseIndexSet& b)
{
    logical_and(b);
}

inline
void
ShortDenseIndexSet::set_difference(const ShortDenseIndexSet& b)
{
    block &= ~b.block;
}

inline
void
ShortDenseIndexSet::set_complement()
{
    logical_not();
}

inline
void
ShortDenseIndexSet::extend(const ShortDenseIndexSet& b1, ShortDenseIndexSet& b2)
{
    assert(b1.size <= b2.size);
    b2.block = b1.block;
}

inline
void
ShortDenseIndexSet::shrink(const ShortDenseIndexSet& b1, ShortDenseIndexSet& b2)
{
    assert(b1.size >= b2.size);
    b2.block = b1.block;
}

inline
Size
ShortDenseIndexSet::count() const
{
    // The following section of code only works for 64 bit block sizes which we can
    // assume (hopefully) because the BlockType is uint64_t.
    // The following code is based upon code obtained from the website
    // http://graphics.stanford.edu/~seander/bithacks.html
    BlockType const w = block - ((block >> 1) & (BlockType) 0x5555555555555555ULL);  // temp
    BlockType const x = (w & (BlockType) 0x3333333333333333ULL) +
            ((w >> 2) & (BlockType) 0x3333333333333333ULL);     // temp
    Size c = (((x + (x >> 4)) & (BlockType) 0x0F0F0F0F0F0F0F0FULL) * (BlockType) 0x0101010101010101ULL)
            >> 56;
    return c;
#if 0
    // The following code is slower than the above code but it has the advantage
    // that it is independent of the block size. It is not used at present.
    Size c = 0;
    const unsigned char* byte = (const unsigned char*) &block;
    const unsigned char* end = byte + BYTES_PER_BLOCK;
    while (byte < end) 
    {
        c += bit_count[*byte];
        ++byte;
    }
    return c;
#endif
}

inline
void
ShortDenseIndexSet::unset_unused_bits()
{
    block &= unused_masks[size];
}

inline
ShortDenseIndexSet::ShortDenseIndexSet()
        : block(0), size(0) {
}

inline
ShortDenseIndexSet::ShortDenseIndexSet(Size _size, bool v)
        : size(_size)
{
    initialise();
    assert(_size >= 0 && _size <= (Size) BITS_PER_BLOCK);
    if (v == false) { zero(); }
    else { one(); } // v == true
}

inline
ShortDenseIndexSet::ShortDenseIndexSet(const ShortDenseIndexSet& b)
        : size(b.size)
{
    *this = b;
}

inline
ShortDenseIndexSet&
ShortDenseIndexSet::operator=(const ShortDenseIndexSet& b)
{
    assert(size == b.size);
    block = b.block;
    return *this;
}

inline
ShortDenseIndexSet::~ShortDenseIndexSet()
{
}

inline
void
ShortDenseIndexSet::resize(Size s)
{
    assert(s <= (Size) BITS_PER_BLOCK);
    size = s;
    unset_unused_bits();
}

} // namespace _4ti2_

#undef ALL_ONES_BLOCK
#undef BITS_PER_BLOCK
#undef BYTES_PER_BLOCK

#endif
