#include <fstream>
#include "groebner/VectorArrayAPI.h"
#include "groebner/VectorArrayStream.h"

using namespace _4ti2_;

VectorArrayAPI::VectorArrayAPI(int num_rows, int num_cols)
    : data(num_rows, num_cols, 0)
{
}

VectorArrayAPI::~VectorArrayAPI()
{
}

int
VectorArrayAPI::get_num_rows() const
{
    return data.get_number();
}

int
VectorArrayAPI::get_num_cols() const
{
    return data.get_size();
}

void
VectorArrayAPI::write(const char* filename) const
{
    std::ofstream out(filename);
    write(out);
}

void
VectorArrayAPI::write(std::ostream& out) const
{
    output(out, data);
}

void
VectorArrayAPI::read(std::istream& in)
{
    in >> data;
}

void
VectorArrayAPI::set_entry_int32_t(int r, int c, const int32_t& value)
{
    convert(value, data[r][c]);
}

void
VectorArrayAPI::get_entry_int32_t(int r, int c, int32_t& value) const
{
    convert(data[r][c], value);
}

void
VectorArrayAPI::set_entry_int64_t(int r, int c, const int64_t& value)
{
    convert(value, data[r][c]);
}

void
VectorArrayAPI::get_entry_int64_t(int r, int c, int64_t& value) const
{
    convert(data[r][c], value);
}

#ifdef _4ti2_GMP_
void
VectorArrayAPI::set_entry_mpz_class(int r, int c, const mpz_class& value)
{
    convert(value, data[r][c]);
}

void
VectorArrayAPI::get_entry_mpz_class(int r, int c, mpz_class& value) const
{
    convert(data[r][c], value);
}
#endif
