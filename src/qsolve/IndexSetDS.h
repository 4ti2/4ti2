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

#ifndef _4ti2_groebner__IndexSetDS_
#define _4ti2_groebner__IndexSetDS_

#include <cassert>
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

class IndexSetDS
{
public:
    explicit IndexSetDS(Size _size, bool v = false);
    IndexSetDS(const IndexSetDS&);
    IndexSetDS& operator=(const IndexSetDS&);
    ~IndexSetDS();

    bool set_subset(const IndexSetDS& b2) const;
    bool set_disjoint(const IndexSetDS& b2) const;

    void set_union(const IndexSetDS& b1, const IndexSetDS& b2);
    void set_intersection(const IndexSetDS& b1, const IndexSetDS& b2);
    void set_difference(const IndexSetDS& b1, const IndexSetDS& b2);
    void set_symm_difference(const IndexSetDS& b1, const IndexSetDS& b2);
    void set_complement(const IndexSetDS& b1);

    void set_union(const IndexSetDS& b);
    void set_intersection(const IndexSetDS& b);
    void set_difference(const IndexSetDS& b);
    void set_symm_difference(const IndexSetDS& b);
    void set_complement();

    void assign(const IndexSetDS& b);

    bool empty() const;
    bool full() const;

    bool singleton() const;
    bool singleton_diff(const IndexSetDS& b) const;
    bool singleton_intersection(const IndexSetDS& b) const;
    Size count() const;
    Size count_union(const IndexSetDS& b) const;
    bool count_lte(Size s) const;
    bool count_lte_diff(Size s, const IndexSetDS& b) const;

    bool operator[](Index index) const;

    void set(Index index);
    void unset(Index index);
    void flip(Index index);
    void zero();
    void one();
    void range(Index start, Index end);

    void swap_odd_n_even();

    Size get_size() const;
    void resize(Size s);

    static void swap(IndexSetDS& b1, IndexSetDS& b2);

    friend bool operator==(const IndexSetDS& b1, const IndexSetDS& b2);
    friend bool operator!=(const IndexSetDS& b1, const IndexSetDS& b2);
    friend bool operator<(const IndexSetDS& b1, const IndexSetDS& b2);

    class Iter
    {
    public:
        operator Index() const { return i; }
        Index operator*() const { return i; }
        bool operator!=(Iter it) const { return i != it.i; }
        Iter& operator=(Iter it) { i = it.i; assert(&is == &it.is); return *this; }
        void operator++() { ++i; while (i != is.size && !is[i]) { ++i; } }
    private:
        Iter(Index _i, const IndexSetDS& _is) : i(_i), is(_is)
            {  assert(i <= is.size); if (i != is.size && !is[i]) { operator++(); } }
        Index i;
        const IndexSetDS& is;
        friend class IndexSetDS;
    };

    Iter begin() const { return Iter(0, *this); }
    Iter end() const { return Iter(size, *this); }

    class NIter
    {
    public:
        operator Index() const { return i; }
        Index operator*() const { return i; }
        bool operator!=(NIter it) const { return i != it.i; }
        NIter& operator=(NIter it) { i = it.i; assert(&is == &it.is); return *this; }
        void operator++() { ++i; while (i != is.size && is[i]) { ++i; } }
    private:
        NIter(Index _i, const IndexSetDS& _is) : i(_i), is(_is)
            { assert(i <= is.size); if (i != is.size && is[i]) { operator++(); } }
        Index i;
        const IndexSetDS& is;
        friend class IndexSetDS;
    };

    NIter nbegin() const { return NIter(0, *this); }
    NIter nend() const { return NIter(size, *this); }

    static const int max_size = BITS_PER_BLOCK;

protected:
    IndexSetDS();

    void logical_ior(const IndexSetDS& b1, const IndexSetDS& b2);
    void logical_xor(const IndexSetDS& b1, const IndexSetDS& b2);
    void logical_and(const IndexSetDS& b1, const IndexSetDS& b2);
    void logical_not(const IndexSetDS& b1);

    void logical_ior(const IndexSetDS& b);
    void logical_xor(const IndexSetDS& b);
    void logical_and(const IndexSetDS& b);
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
IndexSetDS::singleton() const
{
    //BlockType tmp = block;
    //tmp &= tmp-1;
    //return !tmp;
    return !(block & (block-1));
}

inline
bool
IndexSetDS::singleton_diff(const IndexSetDS& b) const
{
    BlockType tmp = block & ~b.block;
    return !(tmp & (tmp-1));
}

inline
bool
IndexSetDS::singleton_intersection(const IndexSetDS& b) const
{
    BlockType tmp = block & b.block;
    return !(tmp & (tmp-1));
}

inline
bool
IndexSetDS::count_lte(Size s) const
{
    BlockType tmp = block;
    while (tmp) {
        if (s == 0) { return false; }
        tmp &= tmp-1;
        --s;
    }
    return true;
}

inline
bool
IndexSetDS::count_lte_diff(Size s, const IndexSetDS& b) const
{
    BlockType tmp = block & ~b.block;
    while (tmp) {
        if (s == 0) { return false; }
        tmp &= tmp-1;
        --s;
    }
    return true;
}

inline
bool
IndexSetDS::empty() const
{
    return !(block);
}

inline
bool
IndexSetDS::full() const
{
    return (block == unused_masks[size]);
}

inline
void
IndexSetDS::set(Index index)
{
    assert(index >= 0 && index < size);
    block |= set_masks[index];
}

inline
void
IndexSetDS::unset(Index index)
{
    assert(index >= 0 && index < size);
    block &= unset_masks[index];
}

inline
void
IndexSetDS::flip(Index index)
{
    assert(index >= 0 && index < size);
    block ^= set_masks[index];
}

inline
void
IndexSetDS::zero()
{
    block = 0;
}

inline
void
IndexSetDS::one()
{
    block = ALL_ONES_BLOCK;
    unset_unused_bits();
}

inline
void
IndexSetDS::range(Index start, Index end)
{
    assert(start <= end && end <= size);
    block = unused_masks[end];
    block &= ~unused_masks[start];
}

inline
void
IndexSetDS::swap_odd_n_even()
{
    assert(size%2 == 0);
    block = ((block >> 1) & (BlockType) 0x5555555555555555ULL) |
            ((block << 1) & (BlockType) 0xAAAAAAAAAAAAAAAAULL);
}

inline
void
IndexSetDS::swap(IndexSetDS& b1, IndexSetDS& b2)
{
    assert(b1.size == b2.size);
    BlockType temp = b1.block;
    b1.block = b2.block;
    b2.block = temp;
}

#if 1
inline
bool
IndexSetDS::operator[](Index index) const
{
    assert(index >= 0 && index <= size);
    return (((BlockType) 1 << index) & block) != 0;
    //return (set_masks[index] & block) != 0;
}
#endif

inline
Size
IndexSetDS::get_size() const
{
    return size;
}

inline
bool
operator<(const IndexSetDS& b1, const IndexSetDS& b2)
{
    assert(b1.size == b2.size);
    return (b1.block < b2.block);
}

inline
bool
operator==(const IndexSetDS& b1, const IndexSetDS& b2)
{
    assert(b1.size == b2.size);
    return (b1.block == b2.block);
}

inline
bool
operator!=(const IndexSetDS& b1, const IndexSetDS& b2)
{
    return !operator==(b1, b2);
}

inline
void
IndexSetDS::logical_ior(const IndexSetDS& b1, const IndexSetDS& b2)
{
    block = b1.block | b2.block;
}

inline
void
IndexSetDS::logical_xor(const IndexSetDS& b1, const IndexSetDS& b2)
{
    block = b1.block ^ b2.block;
}

inline
void
IndexSetDS::logical_and(const IndexSetDS& b1, const IndexSetDS& b2)
{
    block = b1.block & b2.block;
}

inline
void
IndexSetDS::logical_not(const IndexSetDS& b1)
{
    block = ~b1.block;
    unset_unused_bits();
}

inline
void
IndexSetDS::logical_ior(const IndexSetDS& b)
{
    block |= b.block;
}

inline
void
IndexSetDS::logical_xor(const IndexSetDS& b)
{
    block ^= b.block;
}

inline
void
IndexSetDS::logical_and(const IndexSetDS& b)
{
    block &= b.block;
}

inline
void
IndexSetDS::logical_not()
{
    block = ~block;
    unset_unused_bits();
}

inline
bool
IndexSetDS::set_subset(const IndexSetDS& b2) const
{
    assert(size == b2.size);
    return (block == (block & b2.block));
}

inline
bool
IndexSetDS::set_disjoint(const IndexSetDS& b2) const
{
    assert(size == b2.size);
    return ((block & b2.block) == 0);
}

inline
void
IndexSetDS::set_union(const IndexSetDS& b1, const IndexSetDS& b2)
{
    logical_ior(b1, b2);
}

inline
void
IndexSetDS::set_intersection(const IndexSetDS& b1, const IndexSetDS& b2)
{
    logical_and(b1, b2);
}

inline
void
IndexSetDS::set_difference(const IndexSetDS& b1, const IndexSetDS& b2)
{
    block = b1.block & (~b2.block);
}

inline
void
IndexSetDS::set_symm_difference(const IndexSetDS& b1, const IndexSetDS& b2)
{
    block = b1.block ^ b2.block;
}

inline
void
IndexSetDS::set_complement(const IndexSetDS& b1)
{
    logical_not(b1);
}

inline
void
IndexSetDS::set_union(const IndexSetDS& b)
{
    logical_ior(b);
}

inline
void
IndexSetDS::set_intersection(const IndexSetDS& b)
{
    logical_and(b);
}

inline
void
IndexSetDS::set_difference(const IndexSetDS& b)
{
    block &= ~b.block;
}

inline
void
IndexSetDS::set_symm_difference(const IndexSetDS& b)
{
    block ^= b.block;
}

inline
void
IndexSetDS::set_complement()
{
    logical_not();
}

inline
void
IndexSetDS::assign(const IndexSetDS& b)
{
    block = b.block;
    unset_unused_bits();
}

inline
Size
IndexSetDS::count() const
{
    // The following section of code only works for 64 bit block sizes which we can
    // assume (hopefully) because the BlockType is uint64_t.
    // The following code is based upon code obtained from the website
    // http://graphics.stanford.edu/~seander/bithacks.html
    BlockType w = block - ((block >> 1) & (BlockType) 0x5555555555555555ULL);  // temp
    w = (w & (BlockType) 0x3333333333333333ULL) +
            ((w >> 2) & (BlockType) 0x3333333333333333ULL);     // temp
    Size c = ((w + (w >> 4) & (BlockType) 0x0F0F0F0F0F0F0F0FULL) * (BlockType) 0x0101010101010101ULL)
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
Size
IndexSetDS::count_union(const IndexSetDS& b) const
{
    // The following section of code only works for 64 bit block sizes which we can
    // assume (hopefully) because the BlockType is uint64_t.
    // The following code is based upon code obtained from the website
    // http://graphics.stanford.edu/~seander/bithacks.html
    BlockType w = block | b.block;
    w = w - ((w >> 1) & (BlockType) 0x5555555555555555ULL);  // temp
    w = (w & (BlockType) 0x3333333333333333ULL) +
            ((w >> 2) & (BlockType) 0x3333333333333333ULL);     // temp
    Size c = ((w + (w >> 4) & (BlockType) 0x0F0F0F0F0F0F0F0FULL) * (BlockType) 0x0101010101010101ULL)
            >> 56;
    return c;
}

inline
void
IndexSetDS::unset_unused_bits()
{
    block &= unused_masks[size];
}

inline
IndexSetDS::IndexSetDS()
        : block(0), size(0) {
}

inline
IndexSetDS::IndexSetDS(Size _size, bool v)
        : size(_size)
{
    initialise();
    assert(_size >= 0 && _size <= (Size) BITS_PER_BLOCK);
    if (v == false) { zero(); }
    else { one(); } // v == true
}

inline
IndexSetDS::IndexSetDS(const IndexSetDS& b)
        : size(b.size)
{
    *this = b;
}

inline
IndexSetDS&
IndexSetDS::operator=(const IndexSetDS& b)
{
    assert(size == b.size);
    block = b.block;
    return *this;
}

inline
IndexSetDS::~IndexSetDS()
{
}

inline
void
IndexSetDS::resize(Size s)
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
