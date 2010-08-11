#ifndef _4ti2_qsolve__ConeAPI_
#define _4ti2_qsolve__ConeAPI_

#include "qsolve/Size.h"
#include "qsolve/Globals.h"
#include "qsolve/Debug.h"
#include "qsolve/HermiteAlgorithm.h"
#include "qsolve/Vector.h"
#include "qsolve/VectorArray.h"
#include "qsolve/DataType.h"
#include "qsolve/IndexSetR.h"
#include "qsolve/IndexSetD.h"
#include "qsolve/IndexSetDS.h"
#include "qsolve/MatrixWrapper.h"
#include "qsolve/MatrixRef.h"
#include "4ti2/4ti2.h"

namespace _4ti2_ {

class ConeAPI {
public:
    virtual ~ConeAPI();
    virtual void resize(Size m, Size n) = 0;

    // The number of variables.
    virtual Size num_vars() const = 0;
    // The number of constraints.
    virtual Size num_cons() const = 0;

    // Determine whether the given support gives a d-dimensional face of the cone.
    virtual Size is_d_dimensional_face(const IndexSetD& supp) const = 0;
    virtual Size is_d_dimensional_face(const IndexSetDS& supp) const = 0;

    // The constraint matrix A.
    const _4ti2_matrix& get_matrix() const;
    _4ti2_matrix& get_matrix();

    // The set of constraints matching the given constraint type.
    template <class IndexSet>
    void get_constraint_set(_4ti2_constraint t, IndexSet& is) const;

    // Get and set the type of constraint.
    const std::vector<_4ti2_constraint>& get_constraint_types() const;
    _4ti2_constraint get_constraint_type(Index i) const;
    void set_constraint_type(Index i, _4ti2_constraint t);

protected:
    ConeAPI();
    ConeAPI(Size m, Size n);
    ConeAPI(const std::vector<_4ti2_constraint>& cts);
    _4ti2_matrix* matrix_api;
    std::vector<_4ti2_constraint> types;
};

template <class T>
class ConeT : public ConeAPI {
    typedef VectorArrayT<T> MatrixT;
public:
    ConeT();
    virtual ~ConeT();
    ConeT(Size m, Size n);
    virtual void resize(Size m, Size n);

    // The number of variables.
    virtual Size num_vars() const;
    // The number of constraints.
    virtual Size num_cons() const;

    // Determine whether the given support gives a d-dimensional face of the cone.
    virtual Size is_d_dimensional_face(const IndexSetD& supp) const;
    virtual Size is_d_dimensional_face(const IndexSetDS& supp) const;

    void canonize(ConeT<T>& cone, MatrixT& subspace, MatrixT& map) const;

    // The constraint matrix A.
    virtual const MatrixT& get_matrix() const;
    virtual MatrixT& get_matrix();

    // Returns the slacks of the ray r for the constraint i.
    void get_slack(const VectorR<T>& r, Index i, T& d) const;
    // Returns the slacks of the ray r.
    void get_slacks(const VectorR<T>& r, VectorR<T>& s) const;
    // Gets the support of the ray r.
    template <class IS>
    void get_support(const VectorR<T>& r, IS& supp) const;
    // Returns a vector of slacks of the rays rs for the constraint i.
    void get_slacks(const MatrixT& rs, Index i, VectorR<T>& s) const;

    // Counts the polarity of the slacks for the constraint i and the given set of rays.
    void slack_count(const MatrixT& rs, Index i, Size& pos, Size& neg, Size& zero) const;

    // Checks whether the given support determines a d dimensional face of the cone.
    template <class IS>
    Size is_d_dimensional_faceT(const IS& supp) const;

private:
    MatrixT matrix; // Constraints.
};

template <class T, class IndexSet>
class ConeDDT {
public:
    ConeDDT(ConeT<T>& primal, ConeT<T>& dual);

    const ConeT<T>& get_primal();
    const ConeT<T>& get_dual();
    const std::vector<IndexSet>& get_supps();

    // Get the slack of ray r and constraint c.
    void get_slack(Index r, Index c, T& d) const;
private:
    ConeT<T>& primal;
    std::vector<IndexSet> supps;
    ConeT<T>& dual;
};

// Constructor.
inline
ConeAPI::ConeAPI()
    : matrix_api(0)
{
}

// Constructor.
inline
ConeAPI::~ConeAPI()
{
    delete matrix_api;
}

// Constructor.
inline
ConeAPI::ConeAPI(Size m, Size n) 
    : matrix_api(0), types(n+m, _4ti2_LB)
{
}

// Constructor.
inline
ConeAPI::ConeAPI(const std::vector<_4ti2_constraint>& cts)
    : matrix_api(0), types(cts)
{
}

inline
const std::vector<_4ti2_constraint>&
ConeAPI::get_constraint_types() const
{
    return types;
}

inline
_4ti2_constraint
ConeAPI::get_constraint_type(Index i) const
{
    return types[i];
}

inline
void
ConeAPI::set_constraint_type(Index i, _4ti2_constraint t)
{
    types[i] = t;
}

template <class IndexSet> inline
void
ConeAPI::get_constraint_set(_4ti2_constraint t, IndexSet& is) const
{
    assert(is.get_size() == (Size) types.size());
    for (Index i = 0; i < (Index) types.size(); ++i) {
        if (types[i] == t) { is.set(i); }
    }
}

inline
const _4ti2_matrix&
ConeAPI::get_matrix() const
{
    return *matrix_api;
}

inline
_4ti2_matrix&
ConeAPI::get_matrix()
{
    return *matrix_api;
}

// Constructor.
template <class T> inline
ConeT<T>::ConeT() 
    : ConeAPI(), matrix()
{
    ConeAPI::matrix_api = new MatrixRef<MatrixT>(matrix);
}

// Constructor.
template <class T> inline
ConeT<T>::ConeT(Size m, Size n) 
    : ConeAPI(m, n), matrix(m, n)
{
    ConeAPI::matrix_api = new MatrixRef<MatrixT>(matrix);
}

// Constructor.
template <class T> inline
void
ConeT<T>::resize(Size m, Size n) 
{
    matrix.init(m, n);
    ConeAPI::types.resize(n+m, _4ti2_LB);
    ConeAPI::matrix_api = new MatrixRef<MatrixT>(matrix);
}

// Destructor.
template <class T> inline
ConeT<T>::~ConeT() 
{
}

template <class T> inline
const typename ConeT<T>::MatrixT&
ConeT<T>::get_matrix() const
{
    return matrix;
}

template <class T> inline
typename ConeT<T>::MatrixT&
ConeT<T>::get_matrix()
{
    return matrix;
}

template <class T> inline
Size
ConeT<T>::num_vars() const
{
    return matrix.get_size();
}

template <class T> inline
Size
ConeT<T>::num_cons() const
{
    return matrix.get_number();
}

template <class T> inline
void
ConeT<T>::get_slack(const VectorR<T>& v, Index i, T& d) const
{
    if (i < matrix.get_size()) { d = v[i]; }
    else { v.dot(matrix[i-matrix.get_size()], d); }
}

template <class T> inline
void
ConeT<T>::get_slacks(const VectorR<T>& r, VectorR<T>& s) const
{
    assert(s.get_size() == num_vars()+num_cons());
    for (Index i = 0; i < num_vars(); ++i) { s[i] = r[i]; }
    for (Index i = num_vars(); i < num_vars()+num_cons(); ++i) {
        r.dot(matrix[i-num_vars()], s[i]);
    }
}

// Should we set supports for free variables and free constraints?
template <class T>
template <class IS>
void
ConeT<T>::get_support(const VectorR<T>& r, IS& supp) const
{
    assert(supp.get_size() == num_vars()+num_cons());
    supp.zero();
    for (Index i = 0; i < num_vars(); ++i) { if (r[i]!=0) { supp.set(i); } }
    for (Index i = num_vars(); i < num_vars()+num_cons(); ++i) {
        if (r.dot(matrix[i-num_vars()])) { supp.set(i); }
    }
}

template <class T> inline
void
ConeT<T>::get_slacks(const MatrixT& rs, Index i, VectorR<T>& slacks) const
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

template <class T>
Size
ConeT<T>::is_d_dimensional_face(const IndexSetD& supp) const
{
    return is_d_dimensional_faceT(supp);
}

template <class T>
Size
ConeT<T>::is_d_dimensional_face(const IndexSetDS& supp) const
{
    return is_d_dimensional_faceT(supp);
}

// Count how many zero, positive and negative entries there are in a column.
template <class T>
void
ConeT<T>::slack_count(
        const MatrixT& rays, Index next_con,
        Size& pos_count, Size& neg_count, Size& zero_count) const
{
    zero_count = 0; pos_count = 0; neg_count = 0;
    T slack;
    for (Index r = 0; r < rays.get_number(); ++r) {
        get_slack(rays[r], next_con, slack);
        if (slack == 0) { ++zero_count; }
        else if (slack > 0) { ++pos_count; }
        else { ++neg_count; }
    }
}


// Checks whether the given support determines a two dimensional face of the cone.
template <class T>
template <class IndexSet>
Size
ConeT<T>::is_d_dimensional_faceT(const IndexSet& supp) const
{
    // We insert only the constraints that correspond to zero slack entries and
    // only the columns that correspond to the non-zero support entries.
    Index n = num_vars();
    Index m = num_cons();
    assert(supp.get_size() == n+m);
    IndexSet var_supp(n);
    IndexSet con_supp(m);
    for (Index i = 0; i < n; ++i) { if (supp[i]) { var_supp.set(i); } }
    for (Index i = n; i < n+m; ++i) { if (!supp[i]) { con_supp.set(i-n); } }

#if 1
    MatrixT temp(con_supp.count(), var_supp.count());
    temp.assign(matrix, con_supp, var_supp);

    // Compute the rank of the matrix.
    Index rank = upper_triangle(temp);
    return (temp.get_size()-rank);
#else
    MatrixT temp(var_supp.count(), con_supp.count());
    temp.assign_trans(matrix, con_supp, var_supp);

    // Compute the rank of the matrix.
    Index rank = upper_triangle(temp);
    return (temp.get_number()-rank);
#endif
}

template <class T>
void
ConeT<T>::canonize(ConeT<T>& proj_cone, MatrixT& subspace, MatrixT& map) const
{
    DEBUG_4ti2(*out << "MATRIX:\n" << cone.get_matrix() << "\n";)

    Size n = num_vars();
    Size m = num_cons();
    Size num_cons = n+m;

    IndexSetD full_rs(num_cons,0);
    get_constraint_set(_4ti2_LB, full_rs);
    IndexSetD full_cir(num_cons,0);
    get_constraint_set(_4ti2_DB, full_cir);
    IndexSetD full_eq(num_cons,0);
    get_constraint_set(_4ti2_EQ, full_eq);

    DEBUG_4ti2(*out << "RS:\n" << full_rs << "\n";)
    DEBUG_4ti2(*out << "CIR:\n" << full_cir << "\n";)
    DEBUG_4ti2(*out << "EQ:\n" << full_eq << "\n";)

    MatrixT trans(n, num_cons);
    // Add an identity matrix at the beginning of the transpose.
    trans.assignT(0, IndexSetR(0,n), IndexSetR(0,n));
    for (Index i = 0; i < n; ++i) { trans[i][i] = 1; }
    // Add transpose after the identity matrix.
    trans.assign_trans(get_matrix(), IndexSetR(0,m), IndexSetR(0,n), IndexSetR(0,n), IndexSetR(n,n+m));
    DEBUG_4ti2(*out << "TRANS:\n" << trans << "\n";)

    // Process the equality constraints.
    Index eq_row = upper_triangle(trans, 0, n, full_eq.begin(), full_eq.end());
    DEBUG_4ti2(*out << "EQ TRANS:\n" << trans << "\n";)
    trans.remove(0, eq_row);

    // Process the inequality constraints.
    IndexSetD rs_pivots(num_cons, false);
    Index rs_dim = diagonal(trans, 0, trans.get_number(), full_rs.begin(), full_rs.end(), rs_pivots);
    DEBUG_4ti2(*out << "RS PIVOTS:\n" << rs_pivots << "\n";)
    DEBUG_4ti2(*out << "RS TRANS:\n" << trans << "\n";)

    // Process the circuit constraints.
    IndexSetD cir_pivots(num_cons, false);
    Index cir_dim=diagonal(trans, rs_dim, trans.get_number(), full_cir.begin(), full_cir.end(), cir_pivots);
    DEBUG_4ti2(*out << "CIR PIVOTS:\n" << cir_pivots << "\n";)
    DEBUG_4ti2(*out << "CIR TRANS:\n" << trans << "\n";)

    // Extract a linear subspace basis.
    // TODO: Is the subspace basis necessarily independent?
    subspace.init(trans.get_number()-cir_dim, n);
    subspace.assign(trans, IndexSetR(cir_dim, trans.get_number()), IndexSetR(0,n));
    DEBUG_4ti2(*out << "SUB BASIS:\n" << subspace << "\n";)

    // Construct the projected constraint matrix.
    full_rs.set_difference(rs_pivots);
    full_cir.set_difference(cir_pivots);
    proj_cone.resize(full_rs.count()+full_cir.count(), cir_dim);
    MatrixT& proj_matrix = proj_cone.get_matrix();
    // First add rs constraints.
    proj_matrix.assign_trans(trans, IndexSetR(0,cir_dim), full_rs, 
                    IndexSetR(0,full_rs.count()), IndexSetR(0,cir_dim));
    // Then add cir constraints.
    proj_matrix.assign_trans(trans, IndexSetR(0, cir_dim), full_cir, 
                    IndexSetR(full_rs.count(),full_rs.count()+full_cir.count()), IndexSetR(0,cir_dim));
    proj_matrix.normalise();
    DEBUG_4ti2(*out << "PROJ MATRIX:\n" << proj_matrix << "\n";)

    // Construct projected cone and projected initial rays.
    for (Index i = 0; i < rs_dim; ++i) { proj_cone.set_constraint_type(i, _4ti2_LB); }
    for (Index i = rs_dim; i < cir_dim; ++i) { proj_cone.set_constraint_type(i, _4ti2_DB); }
    for (Index i = cir_dim+full_rs.count(); i < cir_dim+full_rs.count()+full_cir.count(); ++i) {
        proj_cone.set_constraint_type(i, _4ti2_DB);
    }
    DEBUG_4ti2(*out << "CONE CONS:\n" << proj_cone.get_constraint_types() << "\n";)

    // The map to lift the rays and circuits back into the original space.
    map.init(n, cir_dim);
    map.assign_trans(trans, IndexSetR(0,cir_dim), IndexSetR(0,n));
    DEBUG_4ti2(*out << "MAP:\n" << map << std::endl;)

    return;
}

} // namespace _4ti2_

#endif
