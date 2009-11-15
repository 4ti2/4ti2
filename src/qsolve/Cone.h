#ifndef _4ti2_qsolve__ConeT_
#define _4ti2_qsolve__ConeT_

#include "qsolve/Size.h"
#include "qsolve/Globals.h"
#include "qsolve/Debug.h"
#include "qsolve/HermiteAlgorithm.h"
#include "qsolve/Vector.h"
#include "qsolve/VectorArray.h"
#include "qsolve/DataType.h"
#include "qsolve/IndexSetR.h"
#include "4ti2/4ti2.h"

namespace _4ti2_ {

// A class representing a cone in the form Ax>=0, x>=0.
template <class T, template <class> class M = VectorArrayT>
class ConeT {
public:
    typedef _4ti2_constraint ConstraintType;

    ConeT();
    ConeT(Size m, Size n);
    void init(Size m, Size n);
    ConeT(const M<T>&);
    ConeT(const M<T>&, const std::vector<ConstraintType>& cts);

    // The number of variables.
    Size num_vars() const;
    // The number of constraints.
    Size num_cons() const;

    // The constraint matrix A.
    const M<T>& get_matrix() const;
    M<T>& get_matrix();
    // The types of constraints (vars then rows).
    const std::vector<ConstraintType>& get_constraint_types() const;
    // Get the type of constraint.
    ConstraintType get_constraint_type(Index i) const;
    // Set the type of constraint.
    void set_constraint_type(Index i, ConstraintType t);

    // Adds the set of constraints matching the given constraint type to the
    // index set is.
    template <class IS>
    void constraint_set(ConstraintType t, IS& is) const;

    // Returns the slacks of the ray r for the constraint i.
    void get_slack(const VectorR<T>& r, Index i, T& d) const;
    // Returns the slacks of the ray r.
    void get_slacks(const VectorR<T>& r, VectorR<T>& s) const;
    // Gets the support of the ray r.
    template <class IS>
    void get_support(const VectorR<T>& r, IS& supp) const;
    // Returns a vector of slacks of the rays rs for the constraint i.
    void get_slacks(const VectorArrayT<T>& rs, Index i, VectorR<T>& s) const;

    // Counts the polarity of the slacks for the constraint i and the given set of rays.
    void slack_count(const VectorArrayT<T>& rs, Index i, Size& pos, Size& neg, Size& zero) const;

    // Checks whether the given support determines a d dimensional face of the cone.
    template <class IS>
    bool is_d_dimensional_face(const IS& supp, int d);

private:
    M<T> matrix;
    std::vector<ConstraintType> types;
};

// Constructor.
template <class T, template <class> class M> inline
ConeT<T,M>::ConeT() 
    : matrix(), types()
{
}

// Constructor.
template <class T, template <class> class M> inline
ConeT<T,M>::ConeT(Size m, Size n) 
    : matrix(m, n), types(n+m, _4ti2_LB)
{
}

// Constructor.
template <class T, template <class> class M> inline
void
ConeT<T,M>::init(Size m, Size n) 
{
    matrix.init(m, n);
    types.resize(n+m, _4ti2_LB);
}

// Constructor.
template <class T, template <class> class M> inline
ConeT<T,M>::ConeT(const M<T>& m) 
    : matrix(m), types(m.get_number()+m.get_size(), _4ti2_LB)
{
}

// Constructor.
template <class T, template <class> class M> inline
ConeT<T,M>::ConeT(const M<T>& m, const std::vector<ConstraintType>& cts)
    : matrix(m), types(cts)
{
}

template <class T, template <class> class M> inline
const M<T>&
ConeT<T,M>::get_matrix() const
{
    return matrix;
}

template <class T, template <class> class M> inline
M<T>&
ConeT<T,M>::get_matrix()
{
    return matrix;
}

template <class T, template <class> class M> inline
const std::vector<typename ConeT<T,M>::ConstraintType>&
ConeT<T,M>::get_constraint_types() const
{
    return types;
}

template <class T, template <class> class M> inline
typename ConeT<T,M>::ConstraintType
ConeT<T,M>::get_constraint_type(Index i) const
{
    return types[i];
}

template <class T, template <class> class M> inline
void
ConeT<T,M>::set_constraint_type(Index i, ConstraintType t)
{
    types[i] = t;
}

template <class T, template <class> class M>
template <class IS>
inline
void
ConeT<T,M>::constraint_set(ConstraintType t, IS& is) const
{
    assert(is.size() == types.size());
    for (Index i = 0; i < (Index) types.size(); ++i) {
        if (types[i] == t) { is.set(i); }
    }
}

template <class T, template <class> class M> inline
Size
ConeT<T,M>::num_vars() const
{
    return matrix.get_size();
}

template <class T, template <class> class M> inline
Size
ConeT<T,M>::num_cons() const
{
    return matrix.get_number();
}

template <class T, template <class> class M> inline
void
ConeT<T,M>::get_slack(const VectorR<T>& v, Index i, T& d) const
{
    if (i < matrix.get_size()) { d = v[i]; }
    else { v.dot(matrix[i-matrix.get_size()], d); }
}

template <class T, template <class> class M> inline
void
ConeT<T,M>::get_slacks(const VectorR<T>& r, VectorR<T>& s) const
{
    assert(s.get_size() == num_vars()+num_cons());
    for (Index i = 0; i < num_vars(); ++i) { s[i] = r[i]; }
    for (Index i = num_vars(); i < num_vars()+num_cons(); ++i) {
        r.dot(matrix[i-num_vars()], s[i]);
    }
}

// Should we set supports for free variables and free constraints?
template <class T, template <class> class M>
template <class IS>
void
ConeT<T,M>::get_support(const VectorR<T>& r, IS& supp) const
{
    assert(s.get_size() == num_vars()+num_cons());
    supp.zero();
    for (Index i = 0; i < num_vars(); ++i) { if (r[i]!=0) { supp.set(i); } }
    for (Index i = num_vars(); i < num_vars()+num_cons(); ++i) {
        if (r.dot(matrix[i-num_vars()])) { supp.set(i); }
    }
}

template <class T, template <class> class M> inline
void
ConeT<T,M>::get_slacks(const VectorArrayT<T>& rs, Index i, VectorR<T>& slacks) const
{
    assert(slacks.get_size() == rs.get_number());
    assert(rs.get_size() == matrix.get_size());
    if (i < matrix.get_size()) {
        for (Index j = 0; j < rs.get_number(); ++j) {
            slacks.set(j,rs[j][i]);
        }
    }
    else {
        const VectorR<T>& v = matrix[i-matrix.get_size()];
        T d;
        for (Index j = 0; j < rs.get_number(); ++j) {
            v.dot(rs[j], d);
            slacks.set(j,d);
        }
    }
}

// Checks whether the given support determines a two dimensional face of the cone.
template <class T, template <class> class M>
template <class IS>
bool
ConeT<T,M>::is_d_dimensional_face(const IS& supp, int d)
{
    assert(supp.size() == num_vars()+num_cons());
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
    //*out << "REDUCED MATRIX:\n" << m << "\n";
    // Compute the rank of the matrix.
    Index rank = upper_triangle(matrix);
    DEBUG_4ti2(*out << "M Rank is " << rank << "\n";)
    if (rank == supp_size-d) { return true; }
    return false;
}

} // namespace _4ti2_

#endif
