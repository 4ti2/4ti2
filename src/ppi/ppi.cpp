// Bestimmung der Primitiven Partitionsidentit"aten (PPI) nach [Urbaniak];
// Verwaltung der Testvektoren mit Range-Trees (BB-alpha based). 

// $Id$

#include <stdio.h>
#include <bool.h>
#include <set>
#include <vector>
#include <iostream.h>
#include <iomanip.h>
#include "bbalpha.h"

#define VVVV 1
#undef TALKATIVE

typedef set<Vector, less<Vector> > SimpleVectorSet;
static vector<Vector *> VectorRepository;

float alpha = 0.15;

#if defined(USE_BB)
class VectorSet {
  BBTree *tree;
  VectorSet(const VectorSet &);
  operator=(const VectorSet &);
public:
  VectorSet(int level) : tree(new BBTree(level, alpha)) {}
  ~VectorSet() { delete tree; }
  typedef BBTree::iterator iterator;
  iterator begin() { return tree->begin(); }
  iterator end() { return tree->end(); }
  int size() const { return tree->root ? tree->root->numLeaves : 0; }
  void insert(Vector v) { 
    Vector *vv = new Vector(v);
    VectorRepository.push_back(vv);
    tree->Insert(Leaf(vv)); 
  }
  bool OrthogonalRangeSearch(Leaf min, Leaf max, 
			     Report report)
  { return tree->OrthogonalRangeSearch(min, max, report); }
  operator SimpleVectorSet();
};

VectorSet::operator SimpleVectorSet() 
{
  SimpleVectorSet S;
  iterator i;
  for (i = begin(); i!=end(); ++i)
    if (!S.insert(*i).second) cerr << "duplicate" << endl;
  return S;
}

#else

// Implementation with digital trees (FIXME: name?)

class Node {
public:
  bool isleaf;
  Node(bool leaf) : isleaf(leaf) {}
};

class InnerNode : public Node {
public:
  vector<Node*> Children;
  InnerNode(int size) : Node(false), Children(size) {}
};

class LeafNode : public Node {
public:
  Leaf leaf;
  LeafNode() : Node(true) {}
};  

class DigitalTree {
  InnerNode *root;
  bool DoSearch(Leaf min, Leaf max, Report report, int level, Node *node);
public:
  int Dimension;
  DigitalTree(int dimension) : 
    Dimension(dimension), 
    root(new InnerNode(2*dimension+1)) {}
  ~DigitalTree() { /* FIXME: */ }
  void Insert(Leaf leaf);
  bool OrthogonalRangeSearch(Leaf min, Leaf max, 
			     Report report)
    { return DoSearch(min, max, report, Dimension-1, root); }
};

bool DigitalTree::DoSearch(Leaf min, Leaf max, 
			   Report report, int level, Node *node)
{
  if (node->isleaf) {
    Vector::pointer vi, mini, maxi;
    int pos = level;
    int count = level+1;
    for (vi = &(((Vector&)(((LeafNode*)node)->leaf))[pos]),
	   mini = &(((Vector&)(min))[pos]),
	   maxi = &(((Vector&)(max))[pos]);
	 count; vi--, mini--, maxi--, count--)
      if ((*vi)<(*mini) || (*vi)>(*maxi)) return true;
    return report(((LeafNode*)node)->leaf);
  }
  else {
    vector<Node*>::pointer ci;
    int count = ((Vector&)max)[level] - ((Vector&)min)[level] + 1;
    for (ci = &(((InnerNode*)node)->Children[Dimension + ((Vector&)min)[level]]); 
	 count; ci++, count--)
      if (*ci && !DoSearch(min, max, report, level-1, *ci)) 
	return false;
    return true;
  }
}

void DigitalTree::Insert(Leaf leaf)
{
  InnerNode *n = root;
  int level;
  for (level = Dimension-1; level>=0; level--) {
    int pos = Dimension + ((Vector&)leaf)[level];
    if (n->Children[pos]) {
      if (n->Children[pos]->isleaf) {
	// Build a chain of inners up to tie 
	LeafNode *ol = (LeafNode*) n->Children[pos];
	do {
	  InnerNode *i = new InnerNode(2 * Dimension + 1);
	  n->Children[pos] = i;
	  n = i, level--, pos = Dimension + ((Vector&)leaf)[level];
	} while (level > 0 
		 && (((Vector&)(ol->leaf))[level] 
		     == ((Vector&)leaf)[level]));
	// Put the old and new leaves at the end
	LeafNode *l = new LeafNode;
	l->leaf = leaf;
	n->Children[Dimension + ((Vector&)leaf)[level]] = l;
	n->Children[Dimension + ((Vector&)(ol->leaf))[level]] = ol;
	// we are done
	return;
      }
      else { // Is inner node
	n = (InnerNode*) n->Children[pos];
	// go on
      }
    }
    else { // Is empty slot
      LeafNode *l = new LeafNode;
      l->leaf = leaf;
      n->Children[pos] = l;
      // we are done
      return; 
    }
  }
}

class VectorSet {
  DigitalTree *tree;
  VectorSet(const VectorSet &);
  operator=(const VectorSet &);
  
public:
  VectorSet(int level) : tree(new DigitalTree(level+1)) {}
  ~VectorSet() { delete tree; }
  typedef BBTree::iterator iterator;
  void insert(Vector v) { 
    Vector *vv = new Vector(v);
    VectorRepository.push_back(vv);
    tree->Insert(Leaf(vv)); 
  }
  bool OrthogonalRangeSearch(Leaf min, Leaf max, 
			     Report report)
    { return tree->OrthogonalRangeSearch(min, max, report); }
  operator SimpleVectorSet();
};

static SimpleVectorSet *inserthackset;

static bool inserthack(const Leaf &y)
{
  inserthackset->insert((Vector&)y);
  return true;
}

VectorSet::operator SimpleVectorSet()
  return Set;
{
  inserthackset = &Set;
  Vector min(tree->Dimension), max(tree->Dimension);
  for (int i = 1; i<=tree->Dimension; i++)
    min(i) = -tree->Dimension, max(i) = tree->Dimension;
  OrthogonalRangeSearch(Leaf(&min), Leaf(&max), inserthack);
}

#endif

ostream &operator<<(ostream &s, const Vector &z)
{
  for (int i = 0; i<z.size(); i++)
    s << setw(4) << (int)(z[i]);
  return s;
}

void ClearVectorRepository()
{
  for (vector<Vector *>::iterator i = VectorRepository.begin();
       i!=VectorRepository.end();
       ++i)
    delete *i;
  VectorRepository = vector<Vector *>();
}

int HilbertDivide(Vector z, Vector y)
{
  // Find maximal integer f with f*y<=z in Hilbert-base sense.
  int maxfactor = INT_MAX;
  for (int i = 0; i<z.size(); i++) {
    if (y[i] > 0) {
      if (y[i] > z[i]) return 0;
      // here is z[i]>=y[i]>0.
      maxfactor = maxfactor <? (z[i] / y[i]);
    }
    else if (y[i] < 0) {
      if (y[i] < z[i]) return 0;
      // here is z[i]<=y[i]<0.
      maxfactor = maxfactor <? (z[i] / y[i]);
    }
  }
  return maxfactor;
}

void writeppi(ostream &c, Vector z, int n)
{
  bool first = true;
  for (int i = 0; i<n; i++) 
    if (z[i] > 0) {
      if (first) first = false;
      else c << " + ";
      for (int j = 1; j < z[i]; j++)
	c << i+1 << " + ";
      c << i+1;
    }
  c << "\t= ";
  first = true;
  for (int i = 0; i<n; i++) 
    if (z[i] < 0) {
      if (first) first = false;
      else c << " + ";
      for (int j = 1; j < -z[i]; j++)
	c << i+1 << " + ";
      c << i+1;
    }
  c << endl;
}

static int count;
static int ppicount;

/* rangereport parameters */
/* FIXME: Use a class instead */
static Vector rangez, rangemin, rangemax;
static int LastNonzeroPos;
#ifdef POSTCHECK
static Vector *RangeException;
#endif

// returns false iff reduced to zero or new range shall be set up.
static bool RangeReportWithInversion(const Leaf &y)
{
  static int count = 0;
#ifdef POSTCHECK
  if (&(Vector&)(y) == RangeException) return true;
#endif
  int maxfactor = HilbertDivide(rangez, y);
  if (maxfactor) {
    //    cerr << "Vector: " << rangez << "Reducer: " << Vector(y) << endl;
#if 0
    count = 0;
#endif
    int i;
    for (i = LastNonzeroPos; 
	 i && !(rangez(i) -= maxfactor * Vector(y)(i)); i--);
    if (i < LastNonzeroPos) { // we have more trailing zeros
      LastNonzeroPos = i;
      if (!i) return false; // reduced to zero
      if (rangez(i)<0) { // we must take the negative
	rangez(i) = -rangez(i);
	for (i--; i; i--) 
	  rangez(i) = maxfactor * Vector(y)(i) - rangez(i);
	return false; // we need a new range
      }
      else {
	for (i--; i; i--) 
	  rangez(i) -= maxfactor * Vector(y)(i);
	return false; // we want a new range
      }
    }
    for (i--; i; i--) 
      rangez(i) -= maxfactor * Vector(y)(i);
    return LastNonzeroPos;
  }
  else {
#if 0 // a heuristic decision whether a new range shall be set up.
    return (++count) % 20;
#endif
#if 1 // always use a single range
    return true;
#else // use a new range for every search.
    return false;
#endif
  }
}

static bool RangeReport(const Leaf &y)
{
#ifdef POSTCHECK
  if (&(Vector&)(y) == RangeException) return true;
#endif
  int maxfactor = HilbertDivide(rangez, y);
  if (maxfactor) {
    int i;
    for (i = LastNonzeroPos; 
	 i && !(rangez(i) -= maxfactor * Vector(y)(i)); i--)
#if defined(DYNAMIC_BOUND)
      rangemax(i) = rangemin(i) = 0;
#else
      ;
#endif
    LastNonzeroPos = i;
#if defined(DYNAMIC_BOUND)
    if (rangez(i)<0) rangemin(i) = rangez(i);
    else rangemax(i) = rangez(i);
#endif
    if (i) {
      for (i--; i; i--) {
	rangez(i) -= maxfactor * Vector(y)(i);
#if defined(DYNAMIC_BOUND)
	if (rangez(i)<0) rangemin(i) = rangez(i);
	else rangemax(i) = rangez(i);
#endif
      }
    }
    return LastNonzeroPos;
  }
  else 
    return true;
}

static bool False(const Leaf &y)
{
#ifdef POSTCHECK
  if (&(Vector&)(y) == RangeException) return true;
#endif
  return false;
}

bool HilbertReduce(Vector &z, VectorSet &S)
{
  int i;
  rangemin = Vector(z.size()), rangemax = Vector(z.size());

  for (i = z.size(); i && !z(i); i--);
  LastNonzeroPos = i;
  rangez = z;

  //
  // Step 1:  Try to reduce the vector to zero, beginning from the
  // last positions, taking care of inversions.
  //

#if 0
  // first look if we find the vector itself in S
  rangemin = rangemax = rangez;
  if (!S.OrthogonalRangeSearch(Leaf(&rangemin), Leaf(&rangemax), &False))
    return false;
#endif
  // positive search
  do {
    for (i = rangez.size(); i!=LastNonzeroPos; i--)
      rangemin(i) = rangemax(i) = 0;
    // at LastNonzeroPos, ensure that we strictly reduce
    rangemin(LastNonzeroPos) = 1, 
      rangemax(LastNonzeroPos) = rangez(LastNonzeroPos);
    // use full range at the remaining positions
    for (i = LastNonzeroPos - 1; i; i--) {
      if (rangez(i) >= 0) rangemin(i) = 0, rangemax(i) = rangez(i);
      else rangemin(i) = rangez(i), rangemax(i) = 0;
    }
  } while (!S.OrthogonalRangeSearch(Leaf(&rangemin), Leaf(&rangemax),
				    &RangeReportWithInversion) 
	   && LastNonzeroPos);
  if (!LastNonzeroPos) return false;

#if 0
  // look if we find the vector itself in S
  rangemin = rangemax = rangez;
  if (!S.OrthogonalRangeSearch(Leaf(&rangemin), Leaf(&rangemax), &False))
    return false;
#endif
  
  // 
  // Step 2:  Reduce the vector as far as possible. No inversion
  // needed beyond this point.
  //

  // zero bounds of trailing zeros and last nonzero pos
  for (i = rangez.size(); i >= LastNonzeroPos; i--) 
    rangemin(i) = rangemax(i) = 0;

  // positive search
  do {
    for (i = LastNonzeroPos - 1; i; i--)
      if (rangez(i) >= 0) rangemin(i) = 0, rangemax(i) = rangez(i);
      else rangemin(i) = rangez(i), rangemax(i) = 0;
  } while (!S.OrthogonalRangeSearch(Leaf(&rangemin), Leaf(&rangemax),
				    &RangeReport));

  // negative search
  rangez = -rangez;
  do {
    for (i = LastNonzeroPos - 1; i; i--)
      if (rangez(i) >= 0) rangemin(i) = 0, rangemax(i) = rangez(i);
      else rangemin(i) = rangez(i), rangemax(i) = 0;
  } while (!S.OrthogonalRangeSearch(Leaf(&rangemin), Leaf(&rangemax),
				    &RangeReport));

  // we have found a new irreducible vector
  z = -rangez;
  return true;
}

/////////////////

// This reduce procedure bails out when value at LastNonzeroPos
// reduced. 
bool HilbertReduceVariant(Vector &z, VectorSet &S)
{
  int i;
  rangemin = Vector(z.size()), rangemax = Vector(z.size());

  for (i = z.size(); i && !z(i); i--);
  LastNonzeroPos = i;
  rangez = z;

  //
  // Step 1:  Try to reduce the value at LastNonzeroPos.
  //

  for (i = rangez.size(); i!=LastNonzeroPos; i--)
    rangemin(i) = rangemax(i) = 0;
  // at LastNonzeroPos, ensure that we strictly reduce
  rangemin(LastNonzeroPos) = 1, 
    rangemax(LastNonzeroPos) = rangez(LastNonzeroPos);
  // use full range at the remaining positions
  for (i = LastNonzeroPos - 1; i; i--) {
    if (rangez(i) >= 0) rangemin(i) = 0, rangemax(i) = rangez(i);
    else rangemin(i) = rangez(i), rangemax(i) = 0;
  }
  if (!S.OrthogonalRangeSearch(Leaf(&rangemin), Leaf(&rangemax), &False))
    return false;

  // 
  // Step 2:  Reduce the vector as far as possible. No inversion
  // needed beyond this point.
  //

  // zero bounds of trailing zeros and last nonzero pos
  for (i = rangez.size(); i >= LastNonzeroPos; i--) 
    rangemin(i) = rangemax(i) = 0;

  // positive search
  do {
    for (i = LastNonzeroPos - 1; i; i--)
      if (rangez(i) >= 0) rangemin(i) = 0, rangemax(i) = rangez(i);
      else rangemin(i) = rangez(i), rangemax(i) = 0;
  } while (!S.OrthogonalRangeSearch(Leaf(&rangemin), Leaf(&rangemax),
				    &RangeReport));

  // negative search
  rangez = -rangez;
  do {
    for (i = LastNonzeroPos - 1; i; i--)
      if (rangez(i) >= 0) rangemin(i) = 0, rangemax(i) = rangez(i);
      else rangemin(i) = rangez(i), rangemax(i) = 0;
  } while (!S.OrthogonalRangeSearch(Leaf(&rangemin), Leaf(&rangemax),
				    &RangeReport));

  // we have found a new irreducible vector
  z = -rangez;
  return true;
}

bool IsReducible(Vector &z, VectorSet &S)
{
  int i;
  rangemin = Vector(z.size()), rangemax = Vector(z.size());
  rangez = z;

  for (i = z.size(); i && !z(i); i--) rangemin(i) = rangemax(i) = 0;
  LastNonzeroPos = i;

  // Let x = z(LastNonzeroPos).

  // positive search is enough because there must be one summand that
  // is positive at LastNonzeroPos.
  for (i--; i; i--)
    if (rangez(i) >= 0) rangemin(i) = 0, rangemax(i) = rangez(i);
    else rangemin(i) = rangez(i), rangemax(i) = 0;

#if 0
  // first check for decomposition `x = 0 + x'
  
  rangemin(LastNonzeroPos) = rangemax(LastNonzeroPos) = z(LastNonzeroPos);
  if (!S.OrthogonalRangeSearch(Leaf(&rangemin), Leaf(&rangemax),
			       &False)) return true;

  // then check for proper decomposition

  rangemin(LastNonzeroPos) = 1, 
    rangemax(LastNonzeroPos) = z(LastNonzeroPos) / 2;
  return (!S.OrthogonalRangeSearch(Leaf(&rangemin), Leaf(&rangemax),
				   &False));

#else
  rangemin(LastNonzeroPos) = 1, 
    rangemax(LastNonzeroPos) = z(LastNonzeroPos);
  return (!S.OrthogonalRangeSearch(Leaf(&rangemin), Leaf(&rangemax),
				   &False));
#endif
}

//
// Specialized code for generating all primitive partition identities
//

void reportx(Vector z)
{
  count++, ppicount++;
#ifdef TALKATIVE
  cout << z << "\t"; writeppi(cout, z, z.size()); 
#endif
}

// Reduces v with repect to P and inserts a nonzero remainder into
// both P and Q.
bool ReduceAndInsert(Vector &v, VectorSet &P, SimpleVectorSet &Q)
{
  // reduce 
  if (HilbertReduceVariant(v, P)) {
    P.insert(v);
    Q.insert(v);
    reportx(v);
    return true;
  }
  return false;
}

SimpleVectorSet ExtendPPI(const SimpleVectorSet &Pn, int n)
{
  VectorSet P(n);
  SimpleVectorSet Pold;
  SimpleVectorSet Pnew;

  // Implementation of Algorithms 4.3.8, 4.3.11 from [Urbaniak] 

  // (1) Create the range-searchable set P of (n+1)-vectors from the
  // set Pn of n-vectors.
  cerr << "# Vectors copied from n = " << n << ": " << endl;
  for (SimpleVectorSet::iterator i = Pn.begin(); i!=Pn.end(); ++i) {
    Vector v(n+1);
    for (int j = 1; j <= n; j++) v(j) = (*i)(j);
    v(n+1) = 0;
    P.insert(v); Pold.insert(v); reportx(v);
  }

  // (2) Add `(n+1) = a + b' identities.
  cerr << "# Vectors of type " << n+1 << " = a + b:" << endl;
  {
    Vector v(n+1);
    for (int j = 1; j <= n; j++) v(j) = 0;
    v(n+1) = +1;
    for (int p = 1; p<=(n+1)/2; p++) {
      // takes care of case: p = n+1-p.
      v(p)--, v(n+1-p)--;
      Pnew.insert(v); 
#if defined(VVVV)
      P.insert(v); reportx(v);
#endif
      v(p)++, v(n+1-p)++;
    }
  }

  // (3) Build all other primitive identities with exactly one
  // component of (n+1).
  cerr << "# Vectors of P" << 1 << "(" << n+1 << "):" << endl;
  for (SimpleVectorSet::iterator i = Pold.begin(); i!=Pold.end(); ++i) {
    Vector v = *i;
    for (int j = 1; j<=(n+1)/2; j++) {
      int k = (n+1) - j;
      if (v(j) != 0 || v(k) != 0) { // otherwise, w reducible by v
	if (v(j) >= 0 || v(k) >= 0) { // otherwise, w reducible by v
	  Vector w = v;
	  w(n+1)++, w(j)--, w(k)--;
#if defined(VVVV)
	  ReduceAndInsert(w, P, Pnew);
#else
	  if (HilbertReduce(w, P)) Pnew.insert(w);
#endif
	}
	if (v(j) <= 0 || v(k) <= 0) { // otherwise, w reducible by v
	  Vector w = -v;
	  w(n+1)++, w(j)--, w(k)--;
#if defined(VVVV)
	  ReduceAndInsert(w, P, Pnew);
#else
	  if (HilbertReduce(w, P)) Pnew.insert(w);
#endif
	}
      }
    }
  }
#if !defined(VVVV)
  // insert them into the searchable vector set
  for (SimpleVectorSet::iterator i = Pnew.begin(); i!=Pnew.end(); ++i) {
    P.insert(*i);
    reportx(*i);
  }
#endif
  //
  Pold = Pnew;
  Pnew = SimpleVectorSet();

  // (4) Build all other primitive identities with exactly (t+1)
  // components of (n+1).
  for (int t = 1; t<n; t++) {
    cerr << "# Vectors of P" << t+1 << "(" << n+1 << "):" << endl;
    for (SimpleVectorSet::iterator i = Pold.begin(); i!=Pold.end(); ++i) {
      Vector v = *i;
      for (int j = 1; j<=(n+1)/2; j++) {
	int k = (n+1) - j;
	if (v(j) != 0 || v(k) != 0) { // otherwise, w reducible by v
	  if (v(j) >= 0 || v(k) >= 0) { // otherwise, w reducible by v
	    Vector w = v;
	    w(n+1)++, w(j)--, w(k)--;
	    ReduceAndInsert(w, P, Pnew);
	  }
	}
      }
    }
    Pold = Pnew;
    Pnew = SimpleVectorSet();
  }

#ifdef POSTCHECK
  cerr << "# Postcheck..." << endl;
  SimpleVectorSet S;
  VectorSet::iterator i;
  for (i = P.begin(); i!=P.end(); ++i) {
    RangeException = &((Vector&)(*i));
    Vector v = *i;
    if (IsReducible(v, P)) {
      cerr << "Haha... a reducible vector: " << *i << endl;
      ppicount--;
    }
    else {
      if (!S.insert(*i).second) {
	cerr << "Hoho... a duplicate: " << *i << endl;
      }
    }
  }
  RangeException = 0;
  return S;
#else
  return P; 
#endif

}

int main(int argc, char *argv[])
{
  // PPI n = 5.
  int n = 0;
  if (argc >= 2) {
    sscanf(argv[1], "%d", &n);
    if (argc == 3)       
      sscanf(argv[2], "%f", &alpha);
    //sscanf(argv[2], "%d", &BBTree::FewLeavesBound);
  }
  if (!n) n = 5;

  // Setup PPI set for n=2
  SimpleVectorSet V;
  Vector v(2); v(1) = -2, v(2) = +1; V.insert(v);
  for (int i = 2; i<n; i++) {
    ppicount = 0;
    cerr << "### Extending to n = " << i+1 << endl;
    V = ExtendPPI(V, i);
    ClearVectorRepository();
    cerr << "### This makes " << ppicount 
	 << " PPI up to sign, " << V.size() << endl;
  }

}

/* $Log$
 * Revision 1.20  1999/03/08 20:58:23  mkoeppe
 * Ok, 11 takes 53sec.
 *
 * Revision 1.19  1999/03/08 16:59:16  mkoeppe
 * Some optimization.
 *
 * Revision 1.18  1999/03/07 18:59:41  mkoeppe
 * Added some special-purpose reduction code, for use in iterative
 * computation (4.3.11) and post-check for reducibility.
 *
 * Revision 1.17  1999/03/07 16:02:19  mkoeppe
 * n=14 still need POST_CHECK, but everything else quite ok now.
 *
 * Revision 1.16  1999/03/06 16:59:20  mkoeppe
 * Added reducibility checks.
 *
 * Revision 1.15  1999/03/05 17:39:27  mkoeppe
 * Remove obsolete check. Gives back memory.
 *
 * Revision 1.14  1999/03/05 17:09:51  mkoeppe
 * Now the last nonzero component is always positive. Took some more care
 * in the positive/negative reduction code. Hopefully this will fix n=14.
 *
 * Revision 1.13  1999/03/05 14:36:17  mkoeppe
 * Minor changes.
 *
 * Revision 1.12  1999/03/05 12:33:27  mkoeppe
 * Unified the search modes; they perform about equally.
 *
 * Revision 1.11  1999/03/04 23:53:02  mkoeppe
 * Initial implementation of R. Urbaniak's PPI algorithms.
 * */
