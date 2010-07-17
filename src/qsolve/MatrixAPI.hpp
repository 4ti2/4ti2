#include "qsolve/VectorArrayStream.h"
#include "qsolve/TypeConversion.h"

namespace _4ti2_ {

template <class MatrixT>
MatrixAPI<MatrixT>::MatrixAPI()
{
}

template <class MatrixT>
MatrixAPI<MatrixT>::MatrixAPI(Size m, Size n)
{
    data.init(m, n);
}

template <class MatrixT>
void
MatrixAPI<MatrixT>::resize(Size m, Size n)
{
    data.init(m, n);
}

template <class MatrixT>
MatrixAPI<MatrixT>::~MatrixAPI()
{
}

template <class MatrixT>
int
MatrixAPI<MatrixT>::get_num_rows() const
{
    return data.get_number();
}

template <class MatrixT>
int
MatrixAPI<MatrixT>::get_num_cols() const
{
    return data.get_size();
}

template <class MatrixT>
void
MatrixAPI<MatrixT>::write(std::ostream& out) const
{
    output(out, data);
}

template <class MatrixT>
void
MatrixAPI<MatrixT>::read(std::istream& in)
{
    in >> data;
}

template <class MatrixT>
void
MatrixAPI<MatrixT>::assign(const _4ti2_matrix& m)
{
    resize(m.get_num_rows(), m.get_num_cols());
    // TODO: Handle other transfer types.
    const MatrixAPI<MatrixT>& vs = dynamic_cast<const MatrixAPI<MatrixT>&>(m);
    data = vs.data;
}

template <class MatrixT>
void
MatrixAPI<MatrixT>::swap(_4ti2_matrix& m)
{
    // TODO: Handle other swap types.
    MatrixAPI<MatrixT>& vs = dynamic_cast<MatrixAPI<MatrixT>&>(m);
    data.swap(vs.data);
}

#define _MatrixAPI_define_(TYPE) \
template <class MatrixT> void MatrixAPI<MatrixT>::set_entry(int r, int c, const TYPE& v) \
{ type_conversion(v, data(r,c)); } \
template <class MatrixT> void MatrixAPI<MatrixT>::get_entry(int r, int c, TYPE& v) const \
{ type_conversion(data(r,c), v); } \

_MatrixAPI_define_(int32_t)
_MatrixAPI_define_(int64_t)
#ifdef _4ti2_GMP_
_MatrixAPI_define_(mpz_class)
#endif

template <class MatrixT>
void
MatrixAPI<MatrixT>::set_entry(int r, int c, std::istream& in)
{
    in >> data(r,c);
}

template <class MatrixT>
void
MatrixAPI<MatrixT>::get_entry(int r, int c, std::ostream& out) const
{
    out << data(r,c);
}

} // namespace _4ti2_
