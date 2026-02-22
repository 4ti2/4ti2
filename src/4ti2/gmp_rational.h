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

#ifndef _4ti2_GMP_RATIONAL_H
#define _4ti2_GMP_RATIONAL_H

#include <gmp.h>
#include <cstring>
#include <ostream>
#include <stdexcept>

#include "4ti2/gmp_integer.h"

namespace _4ti2_GMP_RATIONAL
{

// Minimal exact rational type wrapping mpq_t (GMP C API).
// Only the operations actually used via RationalType are provided.
class Rational
{
public:
    Rational() noexcept
    {
        mpq_init(value_);
    }

    // NOLINTNEXTLINE(google-explicit-constructor)
    Rational(long long value) noexcept
    {
        mpq_init(value_);
        mpz_set_int64(mpq_numref(value_), static_cast<int64_t>(value));
        // denominator is already 1 from mpq_init
    }

    // Overloads that prevent int/long ambiguity between long long and double.
    // NOLINTNEXTLINE(google-explicit-constructor)
    Rational(int value) noexcept
    {
        mpq_init(value_);
        mpz_set_si(mpq_numref(value_), static_cast<long>(value));
    }

    // NOLINTNEXTLINE(google-explicit-constructor)
    Rational(long value) noexcept
    {
        mpq_init(value_);
        mpz_set_si(mpq_numref(value_), value);
    }

    // NOLINTNEXTLINE(google-explicit-constructor)
    Rational(double value) noexcept
    {
        mpq_init(value_);
        mpq_set_d(value_, value);
    }

    // Construct from GMP Integer (exact: sets value to integer/1).
    // NOLINTNEXTLINE(google-explicit-constructor)
    Rational(const _4ti2_GMP_INTEGER::Integer& value) noexcept
    {
        mpq_init(value_);
        mpq_set_z(value_, value.get_mpz_t());
    }

    Rational(const Rational& other) noexcept
    {
        mpq_init(value_);
        mpq_set(value_, other.value_);
    }

    Rational(Rational&& other) noexcept
    {
        mpq_init(value_);
        mpq_swap(value_, other.value_);
    }

    Rational& operator=(const Rational& other) noexcept
    {
        mpq_set(value_, other.value_);
        return *this;
    }

    Rational& operator=(Rational&& other) noexcept
    {
        mpq_swap(value_, other.value_);
        return *this;
    }

    Rational& operator=(double value) noexcept
    {
        mpq_set_d(value_, value);
        return *this;
    }

    Rational& operator=(long long value) noexcept
    {
        mpz_set_int64(mpq_numref(value_), static_cast<int64_t>(value));
        mpz_set_ui(mpq_denref(value_), 1);
        return *this;
    }

    Rational& operator=(int value) noexcept
    {
        mpz_set_si(mpq_numref(value_), static_cast<long>(value));
        mpz_set_ui(mpq_denref(value_), 1);
        return *this;
    }

    Rational& operator=(long value) noexcept
    {
        mpz_set_si(mpq_numref(value_), value);
        mpz_set_ui(mpq_denref(value_), 1);
        return *this;
    }

    ~Rational() noexcept
    {
        mpq_clear(value_);
    }

    // -- Accessors ----------------------------------------------------------

    double get_d() const noexcept
    {
        return mpq_get_d(value_);
    }

    mpq_ptr get_mpq_t() noexcept { return value_; }
    mpq_srcptr get_mpq_t() const noexcept { return value_; }

    // -- Compound assignment ------------------------------------------------

    Rational& operator+=(const Rational& rhs) noexcept
    {
        mpq_add(value_, value_, rhs.value_);
        return *this;
    }

    // -- Binary arithmetic --------------------------------------------------

    friend Rational operator+(const Rational& lhs, const Rational& rhs) noexcept
    {
        Rational result;
        mpq_add(result.value_, lhs.value_, rhs.value_);
        return result;
    }

    friend Rational operator*(const Rational& lhs, const Rational& rhs) noexcept
    {
        Rational result;
        mpq_mul(result.value_, lhs.value_, rhs.value_);
        return result;
    }

    friend Rational operator/(const Rational& lhs, const Rational& rhs)
    {
        if (mpq_sgn(rhs.value_) == 0)
        {
            throw std::domain_error(
                "_4ti2_GMP_RATIONAL::Rational: division by zero");
        }
        Rational result;
        mpq_div(result.value_, lhs.value_, rhs.value_);
        return result;
    }

    // -- Comparisons --------------------------------------------------------

    friend bool operator>(const Rational& lhs, const Rational& rhs) noexcept
    {
        return mpq_cmp(lhs.value_, rhs.value_) > 0;
    }

    friend bool operator<(const Rational& lhs, const Rational& rhs) noexcept
    {
        return mpq_cmp(lhs.value_, rhs.value_) < 0;
    }

    friend bool operator>=(const Rational& lhs, const Rational& rhs) noexcept
    {
        return mpq_cmp(lhs.value_, rhs.value_) >= 0;
    }

    friend bool operator<=(const Rational& lhs, const Rational& rhs) noexcept
    {
        return mpq_cmp(lhs.value_, rhs.value_) <= 0;
    }

    friend bool operator==(const Rational& lhs, const Rational& rhs) noexcept
    {
        return mpq_equal(lhs.value_, rhs.value_) != 0;
    }

    friend bool operator!=(const Rational& lhs, const Rational& rhs) noexcept
    {
        return mpq_equal(lhs.value_, rhs.value_) == 0;
    }

    // -- I/O ----------------------------------------------------------------

    friend std::ostream& operator<<(std::ostream& out, const Rational& value)
    {
        char* raw = mpq_get_str(nullptr, 10, value.value_);
        std::string text(raw);
        {
            void (*freefunc)(void*, size_t);
            mp_get_memory_functions(nullptr, nullptr, &freefunc);
            freefunc(raw, text.size() + 1);
        }
        out << text;
        return out;
    }

private:
    mpq_t value_;
};

} // namespace _4ti2_GMP_RATIONAL

#endif
