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
#include "4ti2/4ti2.h"

namespace _4ti2_ {

class ConeAPI {
public:
    // The number of variables.
    virtual Size num_vars() const = 0;
    // The number of constraints.
    virtual Size num_cons() const = 0;

    // Determine whether the given support gives a d-dimensional face of the cone.
    virtual bool is_d_dimensional_face(const IndexSetD& supp, int d) = 0;
    virtual bool is_d_dimensional_face(const IndexSetDS& supp, int d) = 0;

    // The set of constraints matching the given constraint type.
    virtual void get_constraint_set(_4ti2_constraint t, IndexSetD& is) const = 0;
    virtual void get_constraint_set(_4ti2_constraint t, IndexSetDS& is) const = 0;

    // Get and set the type of constraint.
    virtual _4ti2_constraint get_constraint_type(Index i) const = 0;
    virtual void set_constraint_type(Index i, _4ti2_constraint t) = 0;

protected:
    ConeAPI() {}
    virtual ~ConeAPI() {}
};

template <class T>
class ConeT : public ConeAPI {
public:
    ConeT();
    virtual ~ConeT();
    ConeT(Size m, Size n);
    void resize(Size m, Size n);
    ConeT(const VectorArrayT<T>&);
    ConeT(const VectorArrayT<T>&, const std::vector<_4ti2_constraint>& cts);

    // The number of variables.
    virtual Size num_vars() const;
    // The number of constraints.
    virtual Size num_cons() const;

    // Determine whether the given support gives a d-dimensional face of the cone.
    virtual bool is_d_dimensional_face(const IndexSetD& supp, int d);
    virtual bool is_d_dimensional_face(const IndexSetDS& supp, int d);

    // The set of constraints matching the given constraint type.
    virtual void get_constraint_set(_4ti2_constraint t, IndexSetD& is) const;
    virtual void get_constraint_set(_4ti2_constraint t, IndexSetDS& is) const;

    void canonize(ConeT<T>& cone, VectorArrayT<T>& subspace, VectorArrayT<T>& map) const;

    // The constraint matrix A.
    const VectorArrayT<T>& get_matrix() const;
    VectorArrayT<T>& get_matrix();
    // The types of constraints (vars then rows).
    const std::vector<_4ti2_constraint>& get_constraint_types() const;

    // Get the type of constraint.
    virtual _4ti2_constraint get_constraint_type(Index i) const;
    // Set the type of constraint.
    virtual void set_constraint_type(Index i, _4ti2_constraint t);

    // Adds the set of constraints matching the given constraint type to the
    // index set is.
    template <class IS>
    void constraint_set(_4ti2_constraint t, IS& is) const;

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

private:
    // Checks whether the given support determines a d dimensional face of the cone.
    template <class IS>
    bool is_d_dimensional_faceT(const IS& supp, int d);

    VectorArrayT<T> matrix; // Constraints.
    VectorArrayT<T> generators; // Generators.
    std::vector<_4ti2_constraint> types;
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
template <class T> inline
ConeT<T>::ConeT() 
    : matrix(), types()
{
}

// Constructor.
template <class T> inline
ConeT<T>::~ConeT() 
{
}

// Constructor.
template <class T> inline
ConeT<T>::ConeT(Size m, Size n) 
    : matrix(m, n), types(n+m, _4ti2_LB)
{
}

// Constructor.
template <class T> inline
void
ConeT<T>::resize(Size m, Size n) 
{
    matrix.init(m, n);
    types.resize(n+m, _4ti2_LB);
}

// Constructor.
template <class T> inline
ConeT<T>::ConeT(const VectorArrayT<T>& m) 
    : matrix(m), types(m.get_number()+m.get_size(), _4ti2_LB)
{
}

// Constructor.
template <class T> inline
ConeT<T>::ConeT(const VectorArrayT<T>& m, const std::vector<_4ti2_constraint>& cts)
    : matrix(m), types(cts)
{
}

template <class T> inline
const VectorArrayT<T>&
ConeT<T>::get_matrix() const
{
    return matrix;
}

template <class T> inline
VectorArrayT<T>&
ConeT<T>::get_matrix()
{
    return matrix;
}

template <class T> inline
const std::vector<_4ti2_constraint>&
ConeT<T>::get_constraint_types() const
{
    return types;
}

template <class T> inline
_4ti2_constraint
ConeT<T>::get_constraint_type(Index i) const
{
    return types[i];
}

template <class T> inline
void
ConeT<T>::set_constraint_type(Index i, _4ti2_constraint t)
{
    types[i] = t;
}

template <class T> inline
void
ConeT<T>::get_constraint_set(_4ti2_constraint t, IndexSetD& is) const
{
    constraint_set(t, is);
}

template <class T> inline
void
ConeT<T>::get_constraint_set(_4ti2_constraint t, IndexSetDS& is) const
{
    constraint_set(t, is);
}

template <class T>
template <class IS>
inline
void
ConeT<T>::constraint_set(_4ti2_constraint t, IS& is) const
{
    assert(is.get_size() == (Size) types.size());
    for (Index i = 0; i < (Index) types.size(); ++i) {
        if (types[i] == t) { is.set(i); }
    }
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
ConeT<T>::get_slacks(const VectorArrayT<T>& rs, Index i, VectorR<T>& slacks) const
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
bool 
ConeT<T>::is_d_dimensional_face(const IndexSetD& supp, int d)
{
    return is_d_dimensional_faceT(supp, d);
}

template <class T>
bool 
ConeT<T>::is_d_dimensional_face(const IndexSetDS& supp, int d)
{
    return is_d_dimensional_faceT(supp, d);
}

// Count how many zero, positive and negative entries there are in a column.
template <class T>
void
ConeT<T>::slack_count(
        const VectorArrayT<T>& rays, Index next_con,
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


#if 0
// Counts the polarity of the slacks for the ith constraint for the extreme rays of the cone.
virtual void slack_count(Index i, Size& pos, Size& neg, Size& zero) const;

// Count how many zero, positive and negative entries there are in a column.
template <class T>
void
ConeT<T>::slack_count(Index next_con, Size& pos_count, Size& neg_count, Size& zero_count) const
{
    slack_count(rays, next_con, pos_count, neg_count, zero_count);
}
#endif

// Checks whether the given support determines a two dimensional face of the cone.
template <class T>
template <class IS>
bool
ConeT<T>::is_d_dimensional_faceT(const IS& supp, int d)
{
#if 0 // TODO
    assert(supp.size() == num_vars()+num_cons());
    Index n = num_vars();
    Index m = num_cons();
    // We insert only the constraints that correspond to zero slack entries and
    // only the columns that correspond to the non-zero support entries.
    DEBUG_4ti2(*out << "CONSTRAINTS:\n" << matrix << "\n";)
    IS var_supp(n);
    IS::shrink(supp, var_supp); 
    Index supp_size = var_supp.count();
    VectorArrayT<T> matrix(0, supp_size);
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
#endif
    return false;
}

template <class T>
void
ConeT<T>::canonize(ConeT<T>& proj_cone, VectorArrayT<T>& subspace, VectorArrayT<T>& map) const
{
    DEBUG_4ti2(*out << "MATRIX:\n" << cone.get_matrix() << "\n";)
    //DEBUG_4ti2(*out << "RELS:\n" << rels << "\n";)
    //DEBUG_4ti2(*out << "SIGN:\n" << sign << "\n";)

    Size n = num_vars();
    Size m = num_cons();
    Size num_cons = n+m;

    IndexSetD full_rs(num_cons,0);
    constraint_set(_4ti2_LB, full_rs);
    IndexSetD full_cir(num_cons,0);
    constraint_set(_4ti2_DB, full_cir);
    IndexSetD full_eq(num_cons,0);
    constraint_set(_4ti2_EQ, full_eq);

    DEBUG_4ti2(*out << "RS:\n" << full_rs << "\n";)
    DEBUG_4ti2(*out << "CIR:\n" << full_cir << "\n";)
    DEBUG_4ti2(*out << "EQ:\n" << full_eq << "\n";)

    VectorArrayT<T> trans(n, num_cons);
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
    VectorArrayT<T>& proj_matrix = proj_cone.get_matrix();
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
