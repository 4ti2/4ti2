#ifndef _4ti2_qsolve__ConeC_
#define _4ti2_qsolve__ConeC_

#include "qsolve/Size.h"
#include "qsolve/Globals.h"
#include "qsolve/Debug.h"
#include "qsolve/HermiteAlgorithm.h"
#include "qsolve/Matrix.h"
#include "qsolve/DataType.h"
#include "qsolve/IndexSetR.h"
#include "qsolve/Cone.h"

#undef DEBUG_4ti2
#define DEBUG_4ti2(X) //X

namespace _4ti2_ {

class ConeCAPI {
};

// A class representing a cone in the canonical form Ax>=0, x>=0.
template <class T>
class ConeC : public ConeCAPI {
public:
    ConeC() {}

    // The number of variables.
    Size num_vars() const;
    // The number of constraints.
    Size num_cons() const;

    // The transpose constraint matrix A.
    const MatrixT<T>& get_matrix() const;

    template <class IS>
    void project_cone(const ConeT<T>& trans, const IS& zero_supp, std::vector<Index>& con_map);

    template <class IndexSet>
    bool is_one_dimensional_face(const IndexSet& vars, const IndexSet& cons);

private:
    // The constraint matrix transposed.
    MatrixT<T> matrix;

    // Temporary variable for computation only.
    MatrixT<T> temp; 
};

template <class T> inline
const MatrixT<T>&
ConeC<T>::get_matrix() const
{
    return matrix;
}

template <class T> inline
Size
ConeC<T>::num_vars() const
{
    return matrix.get_m();
}

template <class T> inline
Size
ConeC<T>::num_cons() const
{
    return matrix.get_n();
}

template <class T> template <class IndexSet>
void
ConeC<T>::project_cone(const ConeT<T>& cone, const IndexSet& rel, std::vector<Index>& con_map)
{
    Index n = cone.num_vars();
    Index m = cone.num_cons();
    DEBUG_4ti2(*out << "Constraint Matrix:\n" << cone.get_matrix() << "\n";)
    DEBUG_4ti2(*out << "Zero Support:\n" << rel << "\n";)

    IndexSet vars(n, 0);
    IndexSet cons(m, 0);
    con_map.reserve(n+m);
    for (Index i = 0; i < n; ++i) {
        if (rel[i]) { con_map.push_back(i); }
        else { vars.set(i); con_map.push_back(-1); }
    }
    for (Index i = n; i < m+n; ++i) {
        if (rel[i]) { 
            cons.set(i-n);
            con_map.push_back(i);
        }
    }
    DEBUG_4ti2(*out << "Constraint Map:\n" << con_map << "\n";)
    DEBUG_4ti2(*out << "Con projection:\n" << cons << "\n";)
    DEBUG_4ti2(*out << "Var relaxation:\n" << vars << "\n";)

    matrix.init(n, cons.count());
    matrix.assign_trans(cone.get_matrix(), cons, IndexSetR(0,n));
    DEBUG_4ti2(*out << "NEW TRANS:\n"; matrix.print();)

    typename MatrixT<T>::Pivots pivots;
    matrix.row_diagonalise(vars, IndexSetR(0,matrix.get_n()), &pivots);
    DEBUG_4ti2(*out << "Triangle:\n"; matrix.print();)
    vars.set_complement();
    matrix.row_diagonalise(pivots, vars);
    matrix.row_normalise();
    DEBUG_4ti2(*out << "Diagonal:\n"; matrix.print();)

    for (Index i = 0; i < (Index) pivots.size(); ++i) {
        con_map[pivots[i].first] = con_map[pivots[i].second+n];
        con_map[pivots[i].second+n] = -1;
    } 
    DEBUG_4ti2(*out << "Constraint Map:\n" << con_map << "\n";)
}

// Checks whether the given support determines a one dimensional face of the cone.
template <class T> template <class IndexSet>
bool
ConeC<T>::is_one_dimensional_face(const IndexSet& vars, const IndexSet& cons)
{
    Size proj_n = vars.count();
    Size proj_m = cons.count();
    if (proj_n <= 1) { return true; }
    temp.init(proj_n, proj_m);
    temp.assign(matrix, vars, cons);
    DEBUG_4ti2(*out << "MATRIX:\n"; temp.print();)
    Size rank = temp.row_triangulate();
    DEBUG_4ti2(*out << "Rank: " << rank << "\n";)
    if (rank == temp.get_m()-1) { return true; }

    return false;
}

} // namespace _4ti2_

#undef DEBUG_4ti2
#endif
