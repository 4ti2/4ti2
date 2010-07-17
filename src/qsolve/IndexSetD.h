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

#ifndef _4ti2_qsolve__IndexSetD_
#define _4ti2_qsolve__IndexSetD_

#include "qsolve/Size.h"
#include "qsolve/Index.h"

#include <cassert>
#include <inttypes.h>
#include <climits>
#include <cstdlib>

#include "IndexSetR.h"

namespace _4ti2_
{
typedef uint64_t BlockType;
#define ALL_ONES_BLOCK UINT64_MAX
#define BITS_PER_BLOCK (CHAR_BIT*sizeof(BlockType))
#define BYTES_PER_BLOCK (sizeof(BlockType))

class IndexSetD
{
public:
    explicit IndexSetD(Size _size, bool v = false);
    IndexSetD(const IndexSetD&);
    IndexSetD& operator=(const IndexSetD&);
    ~IndexSetD();

    bool set_subset(const IndexSetD& b2) const;
    bool set_disjoint(const IndexSetD& b2) const;

    void set_union(const IndexSetD& b1, const IndexSetD& b2);
    void set_intersection(const IndexSetD& b1, const IndexSetD& b2);
    void set_difference(const IndexSetD& b1, const IndexSetD& b2);
    void set_symm_difference(const IndexSetD& b1, const IndexSetD& b2);
    void set_complement(const IndexSetD& b1);

    void set_union(const IndexSetD& b);
    void set_intersection(const IndexSetD& b);
    void set_difference(const IndexSetD& b);
    void set_symm_difference(const IndexSetD& b);
    void set_complement();

    void assign(const IndexSetD& b);

    bool empty() const;
    bool full() const;

    bool singleton() const;
    bool singleton_diff(const IndexSetD& b) const;
    bool singleton_intersection(const IndexSetD& b) const;
    Size count() const;
    Size count_union(const IndexSetD& b) const;
    bool count_lte(Size s) const;
    bool count_lte_diff(Size s, const IndexSetD& b) const;
    bool count_lte_2_diff(const IndexSetD& b) const;

    bool operator[](Index index) const;

    Size get_size() const;
    void resize(Size s);

    void set(Index index);
    void flip(Index index);
    void unset(Index index);
    void zero();
    void one();
    void range(Index start, Index end);

    void swap_odd_n_even();

    static void swap(IndexSetD& b1, IndexSetD& b2);

    friend bool operator==(const IndexSetD& b1, const IndexSetD& b2);
    friend bool operator!=(const IndexSetD& b1, const IndexSetD& b2);

    class Iter
    {
    public:
        operator Index() const { return i; }
        Index operator*() const { return i; }
        bool operator!=(Iter it) const { return i != it.i; }
        Iter& operator=(Iter it) { i = it.i; assert(&is == &it.is); return *this; }
        void operator++() { ++i; while (i != is.size && !is[i]) { ++i; } }
    private:
        Iter(Index _i, const IndexSetD& _is) : i(_i), is(_is)
            { if (i < is.size && !is[i]) { operator++(); } }
        Index i;
        const IndexSetD& is;
        friend class IndexSetD;
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
        NIter(Index _i, const IndexSetD& _is) : i(_i), is(_is)
            { assert(i < is.size); if (i != is.size && is[i]) { operator++(); } }
        Index i;
        const IndexSetD& is;
        friend class IndexSetD;
    };

    NIter nbegin() const { return NIter(0, *this); }
    NIter nend() const { return NIter(size, *this); }

protected:
    IndexSetD();

    Size get_num_blocks() const;
    BlockType get_block(Index i) const;

    void logical_ior(const IndexSetD& b1, const IndexSetD& b2);
    void logical_xor(const IndexSetD& b1, const IndexSetD& b2);
    void logical_and(const IndexSetD& b1, const IndexSetD& b2);
    void logical_not(const IndexSetD& b1);

    void logical_ior(const IndexSetD& b);
    void logical_xor(const IndexSetD& b);
    void logical_and(const IndexSetD& b);
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
IndexSetD::singleton() const
{
    Index i = 0;
    while (i < num_blocks) {
        if (blocks[i]) {
            if (blocks[i] & (blocks[i]-1)) { return false; }
            ++i;
            break;
        }
        ++i;
    }
    while (i < num_blocks) { if (blocks[i]) { return false; } ++i; }
    return true;
}

inline
bool
IndexSetD::singleton_diff(const IndexSetD& b) const
{
    bool c = false;
    BlockType tmp;
    for (Index i = 0; i < num_blocks; ++i) {
        tmp = blocks[i] & ~b.blocks[i];
        if (tmp) {
            if (c) { return false; }
            tmp &= tmp-1;
            if (tmp) { return false; }
            c = true;
        }
    }
    return true;
}

inline
bool
IndexSetD::singleton_intersection(const IndexSetD& b) const
{
    bool c = false;
    BlockType tmp;
    for (Index i = 0; i < num_blocks; ++i) {
        tmp = blocks[i] & b.blocks[i];
        if (tmp) {
            if (c) { return false; }
            tmp &= tmp-1;
            if (tmp) { return false; }
            c = true;
        }
    }
    return true;
}

inline
bool
IndexSetD::count_lte(Size s) const
{
    for (Index i = 0; i < num_blocks; ++i) {
        BlockType tmp = blocks[i];
        while (tmp) {
            if (s == 0) { return false; }
            tmp &= tmp-1;
            --s;
        }
    }
    return true;
}

inline
bool
IndexSetD::count_lte_diff(Size s, const IndexSetD& b) const
{
    for (Index i = 0; i < num_blocks; ++i) {
        BlockType tmp = blocks[i] & ~b.blocks[i];
        while (tmp) {
            if (s == 0) { return false; }
            tmp &= tmp-1;
            --s;
        }
    }
    return true;
}

inline
bool
IndexSetD::count_lte_2_diff(const IndexSetD& b) const
{
    // TODO:
    Size s = 2;
    for (Index i = 0; i < num_blocks; ++i) {
        BlockType tmp = blocks[i] & ~b.blocks[i];
        while (tmp) {
            if (s == 0) { return false; }
            tmp &= tmp-1;
            --s;
        }
    }
    return true;
}

inline
bool
IndexSetD::empty() const
{
    for (Index i = 0; i < num_blocks; ++i) {
        if (blocks[i]) { return false; }
    }
    return true;
}

inline
bool
IndexSetD::full() const
{
    for (Index i = 0; i < num_blocks-1; ++i) {
        if (blocks[i] != ALL_ONES_BLOCK) { return false; }
    }
    if (num_blocks != 0 && blocks[num_blocks-1] != unused_masks[size]) { return false; }
    return true;
}

inline
void
IndexSetD::set(Index index)
{
    assert(index >= 0 && index < size);
    blocks[index / sizeofstorage] |= set_masks[index % sizeofstorage];
}

inline
void
IndexSetD::unset(Index index)
{
    assert(index >= 0 && index < size);
    blocks[index / sizeofstorage] &= unset_masks[index % sizeofstorage];
}

inline
void
IndexSetD::flip(Index index)
{
    assert(index >= 0 && index < size);
    blocks[index / sizeofstorage] ^= set_masks[index % sizeofstorage];
}

inline
void
IndexSetD::zero()
{
    for (Index i = 0; i < num_blocks; ++i) { blocks[i] = 0; }
}

inline
void
IndexSetD::one()
{
    for (Index i = 0; i < num_blocks; ++i) { blocks[i] = ALL_ONES_BLOCK; }
    unset_unused_bits();
}

inline
void
IndexSetD::range(Index start, Index end)
{
    assert(start <= end && end <= size);
    Size s = (start ? (start-1)/sizeofstorage : 0);
    Size e = (end ? (end-1)/sizeofstorage : 0);
    for (Index i = 0; i < s; ++i) { blocks[i] = 0; }
    for (Index i = s; i < e; ++i) { blocks[i] = ALL_ONES_BLOCK; }
    for (Index i = e; i < num_blocks; ++i) { blocks[i] = 0; }
    blocks[s] &= ~unused_masks[(start-1)%sizeofstorage];
    blocks[e] |= unused_masks[(end-1)%sizeofstorage];
}

inline
void
IndexSetD::swap_odd_n_even()
{
    assert(size%2 == 0);
    for (Index i = 0; i < num_blocks; ++i) {
        blocks[i] = ((blocks[i] << 1) & (BlockType) 0x5555555555555555ULL) |
                    ((blocks[i] >> 1) & (BlockType) 0xAAAAAAAAAAAAAAAAULL);
    }
}

inline
void
IndexSetD::swap(IndexSetD& b1, IndexSetD& b2)
{
    assert(b1.size == b2.size);
    BlockType* temp = b1.blocks;
    b1.blocks = b2.blocks;
    b2.blocks = temp;
}

#if 1
inline
bool
IndexSetD::operator[](Index index) const
{
    assert(index >= 0 && index <= size);
    //return (set_masks[index%sizeofstorage] & blocks[index/sizeofstorage]) != 0;
    //return (((BlockType) 1 << (index & (Index) 63)) & blocks[index >> 6]) != 0;
    return (((BlockType) 1 << (index%sizeofstorage)) & blocks[index/sizeofstorage]) != 0;
}
#endif

inline
Size
IndexSetD::get_size() const
{
    return size;
}

inline
Size
IndexSetD::get_num_blocks() const
{
    return num_blocks;
}

inline
BlockType
IndexSetD::get_block(Index i) const
{
    return blocks[i];
}

inline
bool
operator==(const IndexSetD& b1, const IndexSetD& b2)
{
    assert(b1.size == b2.size);
    assert(b1.num_blocks == b2.num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i) {
        if (b1.blocks[i] != b2.blocks[i]) { return false; }
    }
    return true;
}

inline
bool
operator!=(const IndexSetD& b1, const IndexSetD& b2)
{
    return !operator==(b1, b2);
}

inline
void
IndexSetD::logical_ior(const IndexSetD& b1, const IndexSetD& b2)
{
    assert(b1.num_blocks == b2.num_blocks && b1.num_blocks == num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i) {
        blocks[i] = b1.blocks[i] | b2.blocks[i];
    }
}

inline
void
IndexSetD::logical_xor(const IndexSetD& b1, const IndexSetD& b2)
{
    assert(b1.num_blocks == b2.num_blocks && b1.num_blocks == num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i) {
        blocks[i] = b1.blocks[i] ^ b2.blocks[i];
    }
}

inline
void
IndexSetD::logical_and(const IndexSetD& b1, const IndexSetD& b2)
{
    assert(b1.num_blocks == b2.num_blocks && b1.num_blocks == num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i) {
        blocks[i] = b1.blocks[i] & b2.blocks[i];
    }
}

inline
void
IndexSetD::logical_not(const IndexSetD& b1)
{
    assert(b1.num_blocks == num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i) { blocks[i] = ~b1.blocks[i]; }
    unset_unused_bits();
}

inline
void
IndexSetD::logical_ior(const IndexSetD& b)
{
    assert(num_blocks == b.num_blocks);
    for (Index i = 0; i < num_blocks; ++i) {
        blocks[i] |= b.blocks[i];
    }
}

inline
void
IndexSetD::logical_xor(const IndexSetD& b)
{
    assert(num_blocks == b.num_blocks);
    for (Index i = 0; i < num_blocks; ++i) {
        blocks[i] ^= b.blocks[i];
    }
}

inline
void
IndexSetD::logical_and(const IndexSetD& b)
{
    assert(num_blocks == b.num_blocks);
    for (Index i = 0; i < num_blocks; ++i) {
        blocks[i] &= b.blocks[i];
    }
}

inline
void
IndexSetD::logical_not()
{
    for (Index i = 0; i < num_blocks; ++i) {
        blocks[i] = ~blocks[i];
    }
    unset_unused_bits();
}

inline
bool
IndexSetD::set_subset(const IndexSetD& b2) const
{
    assert(num_blocks == b2.num_blocks);
    for (Index i = 0; i < num_blocks; ++i) {
        if (blocks[i] != (blocks[i] & b2.blocks[i])) return false;
    }
    return true;
}

inline
bool
IndexSetD::set_disjoint(const IndexSetD& b2) const
{
    assert(num_blocks == b2.num_blocks);
    for (Index i = 0; i < num_blocks; ++i) {
        if ((blocks[i] & b2.blocks[i]) != 0) { return false; }
    }
    return true;
}

inline
void
IndexSetD::set_union(const IndexSetD& b1, const IndexSetD& b2)
{
    logical_ior(b1, b2);
}

inline
void
IndexSetD::set_intersection(const IndexSetD& b1, const IndexSetD& b2)
{
    logical_and(b1, b2);
}

inline
void
IndexSetD::set_difference(const IndexSetD& b1, const IndexSetD& b2)
{
    assert(b1.num_blocks == b2.num_blocks && b1.num_blocks == num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i) {
        blocks[i] = b1.blocks[i] & (~b2.blocks[i]);
    }
}

inline
void
IndexSetD::set_symm_difference(const IndexSetD& b1, const IndexSetD& b2)
{
    assert(b1.num_blocks == b2.num_blocks && b1.num_blocks == num_blocks);
    for (Index i = 0; i < b1.num_blocks; ++i) {
        blocks[i] = b1.blocks[i] ^ b2.blocks[i];
    }
}

inline
void
IndexSetD::set_complement(const IndexSetD& b1)
{
    logical_not(b1);
}

inline
void
IndexSetD::set_union(const IndexSetD& b)
{
    logical_ior(b);
}

inline
void
IndexSetD::set_intersection(const IndexSetD& b)
{
    logical_and(b);
}

inline
void
IndexSetD::set_difference(const IndexSetD& b)
{
    assert(num_blocks == b.num_blocks);
    for (Index i = 0; i < num_blocks; ++i) {
        blocks[i] &= ~b.blocks[i];
    }
}

inline
void
IndexSetD::set_symm_difference(const IndexSetD& b)
{
    assert(num_blocks == b.num_blocks);
    for (Index i = 0; i < num_blocks; ++i) {
        blocks[i] ^= b.blocks[i];
    }
}

inline
void
IndexSetD::set_complement()
{
    logical_not();
}

inline
void
IndexSetD::assign(const IndexSetD& b)
{
    int n = (num_blocks < b.num_blocks ? num_blocks : b.num_blocks);
    int i = 0;
    while (i < n) { blocks[i] = b.blocks[i]; ++i; }
    while (i < num_blocks) { blocks[i] = 0; ++i; }
    unset_unused_bits();
}

inline
Size
IndexSetD::count() const
{
    // The following section of code only works for 64 bit block sizes which we can
    // assume (hopefully) because the BlockType is uint64_t.
    // The following code is based upon code obtained from the website
    // http://graphics.stanford.edu/~seander/bithacks.html
    BlockType* block = blocks;
    BlockType* end = blocks + num_blocks;
    Size c = 0;
    while (block != end) {
        BlockType w = *block - ((*block >> 1) & (BlockType) 0x5555555555555555ULL);  // temp
        w = (w & (BlockType) 0x3333333333333333ULL) +
            ((w >> 2) & (BlockType) 0x3333333333333333ULL);     // temp
        c += ((w + (w >> 4) & (BlockType) 0x0F0F0F0F0F0F0F0FULL) * (BlockType) 0x0101010101010101ULL)
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
    while (byte < end) {
        c += bit_count[*byte];
        ++byte;
    }
    return c;
#endif
}

inline
Size
IndexSetD::count_union(const IndexSetD& b) const
{
    // The following section of code only works for 64 bit block sizes which we can
    // assume (hopefully) because the BlockType is uint64_t.
    // The following code is based upon code obtained from the website
    // http://graphics.stanford.edu/~seander/bithacks.html
    Size c = 0;
    BlockType w;
    for (Index i = 0; i < num_blocks; ++i) {
        w = blocks[i] | b.blocks[i];
        w = w - ((w >> 1) & (BlockType) 0x5555555555555555ULL);  // temp
        w = (w & (BlockType) 0x3333333333333333ULL) +
            ((w >> 2) & (BlockType) 0x3333333333333333ULL);     // temp
        c += ((w + (w >> 4) & (BlockType) 0x0F0F0F0F0F0F0F0FULL) * (BlockType) 0x0101010101010101ULL)
                >> 56;
    }
    return c;
}

inline
void
IndexSetD::unset_unused_bits()
{
    if (size) {
        Size unused_index = ((size-1) % sizeofstorage) + 1;
        blocks[num_blocks-1] &= unused_masks[unused_index];
    }
}

// TODO: Change to return (_size-1)/sizeofstorage + 1; Assuming the _size > 0.
inline
Size
IndexSetD::get_num_blocks(Size _size)
{
    if (_size % sizeofstorage == 0) { return _size/sizeofstorage; }
    else { return _size/sizeofstorage+1; }
}

inline
IndexSetD::IndexSetD()
        : blocks(0), size(0), num_blocks(0)
{
    //*out << "IndexSetD()\n";
}

inline
IndexSetD::IndexSetD(Size _size, bool v)
        : size(_size), num_blocks(get_num_blocks(_size))
{
    initialise();
    assert(_size >= 0);
    blocks = new BlockType[num_blocks];

    if (v == false) { zero(); }
    else { one(); } // v == true
}

inline
IndexSetD::IndexSetD(const IndexSetD& b)
        : size(b.size),
          num_blocks(b.num_blocks)
{
    blocks = new BlockType[num_blocks];
    *this = b;
}

inline
IndexSetD&
IndexSetD::operator=(const IndexSetD& b)
{
    assert(size == b.size);
    assert(num_blocks == b.num_blocks);
    for (Size i = 0; i < num_blocks; ++i) { blocks[i] = b.blocks[i]; }
    return *this;
}

inline
IndexSetD::~IndexSetD()
{
    delete [] blocks;
}


} // namespace _4ti2_

#undef ALL_ONES_BLOCK
#undef BITS_PER_BLOCK
#undef BYTES_PER_BLOCK

#endif
