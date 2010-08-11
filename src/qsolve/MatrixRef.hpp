#include "qsolve/VectorArrayStream.h"
#include "qsolve/TypeConversion.h"

namespace _4ti2_ {

template <class MatrixT>
MatrixRef<MatrixT>::MatrixRef(MatrixT& _m)
    : data(_m)
{
}

template <class MatrixT>
void
MatrixRef<MatrixT>::resize(Size m, Size n)
{
    data.init(m, n);
}

template <class MatrixT>
MatrixRef<MatrixT>::~MatrixRef()
{
}

template <class MatrixT>
int
MatrixRef<MatrixT>::get_num_rows() const
{
    return data.get_number();
}

template <class MatrixT>
int
MatrixRef<MatrixT>::get_num_cols() const
{
    return data.get_size();
}

template <class MatrixT>
void
MatrixRef<MatrixT>::write(std::ostream& out) const
{
    output(out, data);
}

template <class MatrixT>
void
MatrixRef<MatrixT>::read(std::istream& in)
{
    in >> data;
}

template <class MatrixT>
void
MatrixRef<MatrixT>::assign(const _4ti2_matrix& m)
{
    resize(m.get_num_rows(), m.get_num_cols());
    // TODO: Handle other transfer types.
    const MatrixRef<MatrixT>& vs = dynamic_cast<const MatrixRef<MatrixT>&>(m);
    data = vs.data;
}

template <class MatrixT>
void
MatrixRef<MatrixT>::swap(_4ti2_matrix& m)
{
    // TODO: Handle other swap types.
    MatrixRef<MatrixT>& vs = dynamic_cast<MatrixRef<MatrixT>&>(m);
    data.swap(vs.data);
}

#define _MatrixRef_define_(TYPE) \
template <class MatrixT> void MatrixRef<MatrixT>::set_entry(int r, int c, const TYPE& v) \
{ type_conversion(v, data(r,c)); } \
template <class MatrixT> void MatrixRef<MatrixT>::get_entry(int r, int c, TYPE& v) const \
{ type_conversion(data(r,c), v); } \

_MatrixRef_define_(int32_t)
_MatrixRef_define_(int64_t)
#ifdef _4ti2_GMP_
_MatrixRef_define_(mpz_class)
#endif

template <class MatrixT>
void
MatrixRef<MatrixT>::set_entry(int r, int c, std::istream& in)
{
    in >> data(r,c);
}

template <class MatrixT>
void
MatrixRef<MatrixT>::get_entry(int r, int c, std::ostream& out) const
{
    out << data(r,c);
}

} // namespace _4ti2_
