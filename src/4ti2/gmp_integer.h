/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2026 4ti2 team.

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

#ifndef _4ti2_GMP_INTEGER_H
#define _4ti2_GMP_INTEGER_H

#include <gmp.h>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ios>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>

// ---------------------------------------------------------------------------
// Portable int64_t <-> mpz_t helpers.
// These never assume sizeof(long) == 8; they work correctly on both
// LP64 (long = 64-bit) and LLP64/Windows (long = 32-bit).
// ---------------------------------------------------------------------------

// Set mpz_t from int64_t using mpz_import.
inline void mpz_set_int64(mpz_ptr z, int64_t v)
{
    // Two's-complement trick: avoids signed overflow when v == INT64_MIN.
    // -(v+1)+1 == -v for all v, but -(v+1) never overflows since v+1 >= INT64_MIN+1.
    uint64_t uv = v < 0 ? static_cast<uint64_t>(-(v + 1)) + 1u
                         : static_cast<uint64_t>(v);
    mpz_import(z, 1, 1, sizeof(uv), 0, 0, &uv);
    if (v < 0) mpz_neg(z, z);
}

// Set mpz_t from uint64_t using mpz_import.
inline void mpz_set_uint64(mpz_ptr z, uint64_t v)
{
    mpz_import(z, 1, 1, sizeof(v), 0, 0, &v);
}

// Check if mpz_t value fits in int64_t.
inline int mpz_fits_int64_p(mpz_srcptr v)
{
    if (mpz_sgn(v) == 0) return 1;
    // mpz_sizeinbase(v, 2) may overestimate by 1: a reported size of N means
    // the real bit-length is in [N-1, N].  We reject reported sizes > 65
    // because real bit-length >= 65 exceeds int64_t range in all cases.
    // A reported size of exactly 65 may correspond to real bit-length 64,
    // which still needs the exact mpz_export check below.
    if (mpz_sizeinbase(v, 2) > 65) return 0;
    // Export absolute value.  At most 65 real bits => at most 2 uint64_t words.
    uint64_t buf[2] = {0, 0};
    size_t count = 0;
    mpz_export(buf, &count, 1, sizeof(uint64_t), 0, 0, v);
    if (count == 0) return 1;
    if (count > 1) return 0;
    // Positive: must be <= INT64_MAX (2^63 - 1).
    // Negative: magnitude must be <= 2^63.
    return mpz_sgn(v) > 0 ? buf[0] <= static_cast<uint64_t>(INT64_MAX)
                           : buf[0] <= static_cast<uint64_t>(INT64_MAX) + 1u;
}

// Extract int64_t from mpz_t.
// Precondition: mpz_fits_int64_p(v) must be true.
// If violated, returns 0 (no buffer overflow, but the result is meaningless).
inline int64_t mpz_get_int64(mpz_srcptr v)
{
    // Use a 2-element buffer so that even if the precondition is violated,
    // mpz_export cannot write past the end of the buffer.
    uint64_t buf[2] = {0, 0};
    size_t count = 0;
    mpz_export(buf, &count, 1, sizeof(uint64_t), 0, 0, v);
    if (count == 0) return 0;
    if (count > 1) return 0;  // precondition violated
    uint64_t limb = buf[0];
    if (mpz_sgn(v) > 0) return static_cast<int64_t>(limb);
    // Two's-complement negation: -(limb-1)-1 handles limb == 2^63 (INT64_MIN).
    return -static_cast<int64_t>(limb - 1) - 1;
}

// Compare mpz_t with int64_t without heap allocation.
// Uses mpz_cmp_si for values that fit in long; falls back to
// sign + mpz_export magnitude comparison otherwise.
inline int mpz_cmp_int64(mpz_srcptr z, int64_t v)
{
    // Fast path: v fits in long (always true on LP64).
    if (v >= static_cast<int64_t>(LONG_MIN) && v <= static_cast<int64_t>(LONG_MAX)) {
        return mpz_cmp_si(z, static_cast<long>(v));
    }
    // Heap-free fallback for LLP64 (Windows) where long is 32-bit.
    int zsign = mpz_sgn(z);
    int vsign = (v > 0) ? 1 : -1; // v != 0 since 0 fits in long
    if (zsign != vsign) return (zsign > vsign) ? 1 : -1;
    // Same sign; compare magnitudes.
    // sizeinbase may overestimate by 1; > 65 guarantees |z| >= 2^65 > 2^63.
    if (mpz_sizeinbase(z, 2) > 65) return zsign;
    uint64_t zbuf[2] = {0, 0};
    size_t zcount = 0;
    mpz_export(zbuf, &zcount, 1, sizeof(uint64_t), 0, 0, z);
    if (zcount > 1) return zsign; // |z| > 64 bits > |v|
    uint64_t zmag = (zcount > 0) ? zbuf[0] : 0;
    uint64_t vmag = (v > 0) ? static_cast<uint64_t>(v)
                             : static_cast<uint64_t>(-(v + 1)) + 1u;
    if (zmag < vmag) return -zsign;
    if (zmag > vmag) return zsign;
    return 0;
}

// Add int64_t to mpz_t without heap allocation when possible.
inline void mpz_add_int64(mpz_ptr z, mpz_srcptr a, int64_t v)
{
    if (v >= 0 && static_cast<uint64_t>(v) <= static_cast<uint64_t>(ULONG_MAX)) {
        mpz_add_ui(z, a, static_cast<unsigned long>(v));
    } else if (v < 0) {
        uint64_t mag = static_cast<uint64_t>(-(v + 1)) + 1u;
        if (mag <= static_cast<uint64_t>(ULONG_MAX))
            mpz_sub_ui(z, a, static_cast<unsigned long>(mag));
        else {
            mpz_t tmp; mpz_init(tmp);
            mpz_set_int64(tmp, v);
            mpz_add(z, a, tmp);
            mpz_clear(tmp);
        }
    } else {
        mpz_t tmp; mpz_init(tmp);
        mpz_set_int64(tmp, v);
        mpz_add(z, a, tmp);
        mpz_clear(tmp);
    }
}

// Subtract int64_t from mpz_t without heap allocation when possible.
inline void mpz_sub_int64(mpz_ptr z, mpz_srcptr a, int64_t v)
{
    if (v >= 0 && static_cast<uint64_t>(v) <= static_cast<uint64_t>(ULONG_MAX)) {
        mpz_sub_ui(z, a, static_cast<unsigned long>(v));
    } else if (v < 0) {
        uint64_t mag = static_cast<uint64_t>(-(v + 1)) + 1u;
        if (mag <= static_cast<uint64_t>(ULONG_MAX))
            mpz_add_ui(z, a, static_cast<unsigned long>(mag));
        else {
            mpz_t tmp; mpz_init(tmp);
            mpz_set_int64(tmp, v);
            mpz_sub(z, a, tmp);
            mpz_clear(tmp);
        }
    } else {
        mpz_t tmp; mpz_init(tmp);
        mpz_set_int64(tmp, v);
        mpz_sub(z, a, tmp);
        mpz_clear(tmp);
    }
}

// Multiply mpz_t by int64_t without heap allocation when possible.
inline void mpz_mul_int64(mpz_ptr z, mpz_srcptr a, int64_t v)
{
    if (v >= static_cast<int64_t>(LONG_MIN) && v <= static_cast<int64_t>(LONG_MAX)) {
        mpz_mul_si(z, a, static_cast<long>(v));
    } else {
        mpz_t tmp; mpz_init(tmp);
        mpz_set_int64(tmp, v);
        mpz_mul(z, a, tmp);
        mpz_clear(tmp);
    }
}

namespace _4ti2_GMP_INTEGER
{

class Integer
{
public:
    // -- Constructors / destructor ------------------------------------------

    Integer() noexcept
    {
        mpz_init(value_);
    }

    // Non-explicit to support widespread template code like `T x = 0;`.
    // NOLINTNEXTLINE(google-explicit-constructor)
    Integer(long long value) noexcept
    {
        mpz_init(value_);
        mpz_set_int64(value_, static_cast<int64_t>(value));
    }

    // Named factory for unsigned 64-bit values.
    // A constructor overload would cause ambiguity with Integer(long long)
    // in direct-initialization of int/long/unsigned int arguments, because
    // int->long long and int->unsigned long long are conversions of equal rank.
    static Integer from_uint64(uint64_t value) noexcept
    {
        Integer result;
        mpz_set_uint64(result.value_, value);
        return result;
    }

    Integer(const Integer& other) noexcept
    {
        mpz_init_set(value_, other.value_);
    }

    Integer(Integer&& other) noexcept
    {
        mpz_init(value_);
        mpz_swap(value_, other.value_);
    }

    Integer& operator=(const Integer& other) noexcept
    {
        // mpz_set handles self-assignment correctly.
        mpz_set(value_, other.value_);
        return *this;
    }

    Integer& operator=(Integer&& other) noexcept
    {
        mpz_swap(value_, other.value_);
        return *this;
    }

    Integer& operator=(long long value) noexcept
    {
        mpz_set_int64(value_, static_cast<int64_t>(value));
        return *this;
    }

    ~Integer() noexcept
    {
        mpz_clear(value_);
    }

    friend void swap(Integer& a, Integer& b) noexcept
    {
        mpz_swap(a.value_, b.value_);
    }

    // -- Accessors ----------------------------------------------------------

    void set_mpz(mpz_srcptr value) noexcept
    {
        mpz_set(value_, value);
    }

    mpz_ptr get_mpz_t() noexcept
    {
        return value_;
    }

    mpz_srcptr get_mpz_t() const noexcept
    {
        return value_;
    }

    double get_d() const noexcept
    {
        return mpz_get_d(value_);
    }

    // WARNING: returns long, which is 32-bit on LLP64/Windows.
    // Prefer get_int64() for portable 64-bit extraction.
    long get_si() const noexcept
    {
        return mpz_get_si(value_);
    }

    int64_t get_int64() const noexcept
    {
        return mpz_get_int64(value_);
    }

    bool fits_int64() const noexcept
    {
        return mpz_fits_int64_p(value_) != 0;
    }

    explicit operator double() const noexcept
    {
        return mpz_get_d(value_);
    }

    // -- Unary --------------------------------------------------------------

    Integer operator-() const noexcept
    {
        Integer result;
        mpz_neg(result.value_, value_);
        return result;
    }

    // -- Compound assignment: Integer ---------------------------------------

    Integer& operator+=(const Integer& rhs) noexcept
    {
        mpz_add(value_, value_, rhs.value_);
        return *this;
    }

    Integer& operator-=(const Integer& rhs) noexcept
    {
        mpz_sub(value_, value_, rhs.value_);
        return *this;
    }

    Integer& operator*=(const Integer& rhs) noexcept
    {
        mpz_mul(value_, value_, rhs.value_);
        return *this;
    }

    Integer& operator/=(const Integer& rhs)
    {
        if (mpz_sgn(rhs.value_) == 0) {
            throw std::domain_error("_4ti2_GMP_INTEGER::Integer: division by zero");
        }
        mpz_tdiv_q(value_, value_, rhs.value_);
        return *this;
    }

    Integer& operator%=(const Integer& rhs)
    {
        if (mpz_sgn(rhs.value_) == 0) {
            throw std::domain_error("_4ti2_GMP_INTEGER::Integer: division by zero");
        }
        mpz_tdiv_r(value_, value_, rhs.value_);
        return *this;
    }

    // -- Compound assignment: long long (heap-free when value fits) ---------

    Integer& operator+=(long long rhs) noexcept
    {
        mpz_add_int64(value_, value_, static_cast<int64_t>(rhs));
        return *this;
    }

    Integer& operator-=(long long rhs) noexcept
    {
        mpz_sub_int64(value_, value_, static_cast<int64_t>(rhs));
        return *this;
    }

    Integer& operator*=(long long rhs) noexcept
    {
        mpz_mul_int64(value_, value_, static_cast<int64_t>(rhs));
        return *this;
    }

    // -- Binary arithmetic: Integer x Integer -------------------------------

    friend Integer operator+(Integer lhs, const Integer& rhs) noexcept
    {
        lhs += rhs;
        return lhs;
    }

    friend Integer operator-(Integer lhs, const Integer& rhs) noexcept
    {
        lhs -= rhs;
        return lhs;
    }

    friend Integer operator*(Integer lhs, const Integer& rhs) noexcept
    {
        lhs *= rhs;
        return lhs;
    }

    friend Integer operator/(Integer lhs, const Integer& rhs)
    {
        lhs /= rhs;
        return lhs;
    }

    friend Integer operator%(Integer lhs, const Integer& rhs)
    {
        lhs %= rhs;
        return lhs;
    }

    // -- Binary arithmetic: Integer x long long (and reverse) ---------------

    friend Integer operator+(Integer lhs, long long rhs) noexcept
    {
        lhs += rhs;
        return lhs;
    }

    friend Integer operator+(long long lhs, Integer rhs) noexcept
    {
        rhs += lhs;
        return rhs;
    }

    friend Integer operator-(Integer lhs, long long rhs) noexcept
    {
        lhs -= rhs;
        return lhs;
    }

    friend Integer operator-(long long lhs, const Integer& rhs) noexcept
    {
        Integer result;
        mpz_set_int64(result.value_, static_cast<int64_t>(lhs));
        mpz_sub(result.value_, result.value_, rhs.value_);
        return result;
    }

    friend Integer operator*(Integer lhs, long long rhs) noexcept
    {
        lhs *= rhs;
        return lhs;
    }

    friend Integer operator*(long long lhs, Integer rhs) noexcept
    {
        rhs *= lhs;
        return rhs;
    }

    // -- Comparisons: Integer x Integer -------------------------------------

    friend bool operator==(const Integer& lhs, const Integer& rhs) noexcept
    {
        return mpz_cmp(lhs.value_, rhs.value_) == 0;
    }

    friend bool operator!=(const Integer& lhs, const Integer& rhs) noexcept
    {
        return mpz_cmp(lhs.value_, rhs.value_) != 0;
    }

    friend bool operator<(const Integer& lhs, const Integer& rhs) noexcept
    {
        return mpz_cmp(lhs.value_, rhs.value_) < 0;
    }

    friend bool operator<=(const Integer& lhs, const Integer& rhs) noexcept
    {
        return mpz_cmp(lhs.value_, rhs.value_) <= 0;
    }

    friend bool operator>(const Integer& lhs, const Integer& rhs) noexcept
    {
        return mpz_cmp(lhs.value_, rhs.value_) > 0;
    }

    friend bool operator>=(const Integer& lhs, const Integer& rhs) noexcept
    {
        return mpz_cmp(lhs.value_, rhs.value_) >= 0;
    }

    // -- Comparisons: Integer x long long (heap-free via mpz_cmp_int64) -----

    friend bool operator==(const Integer& lhs, long long rhs) noexcept
    {
        return mpz_cmp_int64(lhs.value_, static_cast<int64_t>(rhs)) == 0;
    }

    friend bool operator==(long long lhs, const Integer& rhs) noexcept
    {
        return mpz_cmp_int64(rhs.value_, static_cast<int64_t>(lhs)) == 0;
    }

    friend bool operator!=(const Integer& lhs, long long rhs) noexcept
    {
        return mpz_cmp_int64(lhs.value_, static_cast<int64_t>(rhs)) != 0;
    }

    friend bool operator!=(long long lhs, const Integer& rhs) noexcept
    {
        return mpz_cmp_int64(rhs.value_, static_cast<int64_t>(lhs)) != 0;
    }

    friend bool operator<(const Integer& lhs, long long rhs) noexcept
    {
        return mpz_cmp_int64(lhs.value_, static_cast<int64_t>(rhs)) < 0;
    }

    friend bool operator<(long long lhs, const Integer& rhs) noexcept
    {
        return mpz_cmp_int64(rhs.value_, static_cast<int64_t>(lhs)) > 0;
    }

    friend bool operator<=(const Integer& lhs, long long rhs) noexcept
    {
        return mpz_cmp_int64(lhs.value_, static_cast<int64_t>(rhs)) <= 0;
    }

    friend bool operator<=(long long lhs, const Integer& rhs) noexcept
    {
        return mpz_cmp_int64(rhs.value_, static_cast<int64_t>(lhs)) >= 0;
    }

    friend bool operator>(const Integer& lhs, long long rhs) noexcept
    {
        return mpz_cmp_int64(lhs.value_, static_cast<int64_t>(rhs)) > 0;
    }

    friend bool operator>(long long lhs, const Integer& rhs) noexcept
    {
        return mpz_cmp_int64(rhs.value_, static_cast<int64_t>(lhs)) < 0;
    }

    friend bool operator>=(const Integer& lhs, long long rhs) noexcept
    {
        return mpz_cmp_int64(lhs.value_, static_cast<int64_t>(rhs)) >= 0;
    }

    friend bool operator>=(long long lhs, const Integer& rhs) noexcept
    {
        return mpz_cmp_int64(rhs.value_, static_cast<int64_t>(lhs)) <= 0;
    }

    // -- I/O ----------------------------------------------------------------

    friend std::ostream& operator<<(std::ostream& out, const Integer& value)
    {
        // Determine base from stream format flags.
        int base = 10;
        std::ios_base::fmtflags flags = out.flags();
        switch (flags & std::ios_base::basefield) {
        case std::ios_base::hex: base = 16; break;
        case std::ios_base::oct: base = 8;  break;
        default:                 base = 10; break;
        }

        // GMP uses negative base for uppercase hex digits.
        int gmp_base = (base == 16 && (flags & std::ios_base::uppercase))
                        ? -16 : base;

        // Get string representation; copy into std::string immediately
        // and free the GMP buffer, so that no subsequent throw can leak it.
        char* raw = mpz_get_str(nullptr, gmp_base, value.value_);
        std::string text(raw);
        {
            void (*freefunc)(void*, size_t);
            mp_get_memory_functions(nullptr, nullptr, &freefunc);
            freefunc(raw, text.size() + 1);
        }

        // showbase prefix (insert after any leading '-').
        if (flags & std::ios_base::showbase) {
            std::string::size_type pos =
                (!text.empty() && text[0] == '-') ? 1u : 0u;
            if (base == 16) {
                text.insert(pos,
                    (flags & std::ios_base::uppercase) ? "0X" : "0x");
            } else if (base == 8 && text != "0" && text != "-0") {
                text.insert(pos, "0");
            }
        }

        // showpos for non-negative decimal values.
        if ((flags & std::ios_base::showpos) && base == 10
            && mpz_sgn(value.value_) >= 0) {
            text.insert(std::string::size_type(0), "+");
        }

        // The stream's operator<<(const std::string&) handles
        // width, fill, and adjustment.
        out << text;
        return out;
    }

    friend std::istream& operator>>(std::istream& in, Integer& value)
    {
        std::string text;
        in >> text;
        if (!in) {
            return in;
        }
        // Decide base from stream flags.
        int base = 10;
        switch (in.flags() & std::ios_base::basefield) {
        case std::ios_base::hex: base = 16; break;
        case std::ios_base::oct: base = 8;  break;
        default:                 base = 10; break;
        }
        // Parse into a temporary to avoid corrupting value on bad input.
        // GMP documents that mpz_set_str leaves rop in an undefined state
        // on parse failure.
        mpz_t tmp;
        mpz_init(tmp);
        if (mpz_set_str(tmp, text.c_str(), base) != 0) {
            mpz_clear(tmp);
            in.setstate(std::ios::failbit);
        } else {
            mpz_swap(value.value_, tmp);
            mpz_clear(tmp);
        }
        return in;
    }

private:
    mpz_t value_;
};

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

inline Integer gcd(const Integer& a, const Integer& b) noexcept
{
    Integer result;
    mpz_gcd(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return result;
}

// Return the number of bits in the absolute value of `value` (0 for zero).
inline int bit_length(const Integer& value) noexcept
{
    return mpz_sgn(value.get_mpz_t()) == 0
        ? 0
        : static_cast<int>(mpz_sizeinbase(value.get_mpz_t(), 2));
}

// Backward-compatibility alias for bit_length.
inline int calc_precision(const Integer& value) noexcept
{
    return bit_length(value);
}

} // namespace _4ti2_GMP_INTEGER

#endif
