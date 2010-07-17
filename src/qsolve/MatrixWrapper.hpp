#include "qsolve/VectorArrayStream.h"
#include "qsolve/TypeConversion.h"

namespace _4ti2_ {

template <class MatrixT>
MatrixWrapper<MatrixT>::MatrixWrapper()
{
}

template <class MatrixT>
MatrixWrapper<MatrixT>::MatrixWrapper(Size m, Size n)
    : MatrixT(m,n)
{
}

template <class MatrixT>
void
MatrixWrapper<MatrixT>::resize(Size m, Size n)
{
    MatrixT::init(m, n);
}

template <class MatrixT>
MatrixWrapper<MatrixT>::~MatrixWrapper()
{
}

template <class MatrixT>
int
MatrixWrapper<MatrixT>::get_num_rows() const
{
    return MatrixT::get_number();
}

template <class MatrixT>
int
MatrixWrapper<MatrixT>::get_num_cols() const
{
    return MatrixT::get_size();
}

template <class MatrixT>
void
MatrixWrapper<MatrixT>::write(std::ostream& out) const
{
    output(out, *this);
}

template <class MatrixT>
void
MatrixWrapper<MatrixT>::read(std::istream& in)
{
    in >> *this;
}

template <class MatrixT>
void
MatrixWrapper<MatrixT>::assign(const _4ti2_matrix& _m)
{
    resize(_m.get_num_rows(), _m.get_num_cols());
    // TODO: Handle other transfer types.
    const MatrixT& m = dynamic_cast<const MatrixT&>(_m);
    MatrixT::operator=(m);
}

template <class MatrixT>
void
MatrixWrapper<MatrixT>::swap(_4ti2_matrix& _m)
{
    // TODO: Handle other swap types.
    MatrixT& m = dynamic_cast<MatrixT&>(_m);
    MatrixT::swap(m);
}

#define _MatrixWrapper_define_(TYPE) \
template <class MatrixT> void MatrixWrapper<MatrixT>::set_entry(int r, int c, const TYPE& v) \
{ type_conversion(v, MatrixT::vectors[r][c]); } \
template <class MatrixT> void MatrixWrapper<MatrixT>::get_entry(int r, int c, TYPE& v) const \
{ type_conversion(MatrixT::vectors[r][c], v); } \

_MatrixWrapper_define_(int32_t)
_MatrixWrapper_define_(int64_t)
#ifdef _4ti2_GMP_
_MatrixWrapper_define_(mpz_class)
#endif

template <class MatrixT>
void
MatrixWrapper<MatrixT>::set_entry(int r, int c, std::istream& in)
{
    in >> MatrixT::operator()(r,c);
}

template <class MatrixT>
void
MatrixWrapper<MatrixT>::get_entry(int r, int c, std::ostream& out) const
{
    out << MatrixT::operator()(r,c);
}

} // namespace _4ti2_
