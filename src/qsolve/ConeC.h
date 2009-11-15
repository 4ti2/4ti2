#ifndef _4ti2_qsolve__ConeC_
#define _4ti2_qsolve__ConeC_

#include "qsolve/Size.h"
#include "qsolve/Globals.h"
#include "qsolve/Debug.h"
#include "qsolve/HermiteAlgorithm.h"
#include "qsolve/VectorArray.h"
#include "qsolve/DataType.h"
#include "qsolve/IndexSetR.h"

namespace _4ti2_ {

// A class representing a cone in the form Ax>=0, x>=0.
template <class T, template <class> class M = VectorArrayT>
class ConeC {
public:
    ConeC(const M<T>&);

    // The number of variables.
    Size num_vars() const;
    // The number of constraints.
    Size num_cons() const;

    // The transpose constraint matrix A.
    const M<T>& constraint_matrix() const;

    // Checks whether the given support determines a d dimensional face of the cone.
    template <class IS>
    bool is_d_dimensional_face(const IS& supp, int d);
 
private:
    M<T> matrix;
};

// Constructor.
template <class T, template <class> class M> inline
ConeC<T,M>::ConeC(const M<T>& m) : matrix(m)
{
}

template <class T, template <class> class M> inline
const M<T>&
ConeC<T,M>::constraint_matrix() const
{
    return matrix;
}

template <class T, template <class> class M> inline
Size
ConeC<T,M>::num_vars() const
{
    return matrix.get_size();
}

template <class T, template <class> class M> inline
Size
ConeC<T,M>::num_cons() const
{
    return matrix.get_number();
}

// Checks whether the given support determines an d dimensional face of the cone.
template <class T, template <class> class M>
template <class IS>
bool
ConeC<T,M>::is_d_dimensional_face(const IS& supp, int d)
{
    Index n = num_vars();
    Index m = num_cons();
    // We insert only the constraints that correspond to zero slack entries and
    // only the columns that correspond to the non-zero support entries.
    DEBUG_4ti2(*out << "CONSTRAINTS:\n" << matrix << "\n";)
    IS var_supp(n);
    IS::shrink(supp, var_supp); 
    Index supp_size = var_supp.count();
    M<T> matrix(0, supp_size);
    VectorT<T> temp_vec(supp_size);
    DEBUG_4ti2(*out << "SUPP:\n" << supp << "\n";)
    DEBUG_4ti2(*out << "VAR SUPP:\n" << var_supp << "\n";)
    for (Index i = 0; i < m; ++i) { 
        if (!supp[i+n]) {
            temp_vec.assign(matrix[i], var_supp, IndexSetR(0,supp_size));
            matrix.insert(temp_vec);
        }
    }
    DEBUG_4ti2(*out << "REDUCED MATRIX:\n" << matrix << "\n";)

    // Compute the rank of the matrix.
    Index rank = upper_triangle(matrix);
    DEBUG_4ti2(*out << "M Rank is " << rank << "\n";)
    if (rank == supp_size-d) { return true; }
    return false;
}

} // namespace _4ti2_

#endif
