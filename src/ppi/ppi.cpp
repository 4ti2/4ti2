// Berechnung einer Hilbertbasis; Modifikation eines Verfahrens aus
// der Vorlesung Optimierung II (Weismantel). Hier angewendet auf
// Primitive Partitionsidentit"aten (PPI); das Verfahren ist aber
// allgemein implementiert. "

// $Id$

// Verwaltung der Testvektoren mit Range-Trees (BB-alpha based). 

#include <stdio.h>
#include <bool.h>
#include <set>
#include <vector>
#include <iostream.h>
#include <iomanip.h>
#include "bbalpha.h"

typedef set<Vector, less<Vector> > SimpleVectorSet;
static vector<Vector *> VectorRepository;

float alpha = 0.15;

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
			     bool (*report)(const Leaf &)) 
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

ostream &operator<<(ostream &s, const Vector &z)
{
  for (int i = 0; i<z.size(); i++)
    s << setw(4) << z[i];
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
static Vector rangez;
static int LastNonzeroPos;
static bool DidInvert;
static Vector *RangeException;

// returns false iff reduced to zero or new range shall be set up.
static bool rangereport(const Leaf &y)
{
  static int count = 0;
  if (&(Vector&)(y) == RangeException) return true;
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
	DidInvert = true;
	rangez(i) = -rangez(i);
	for (i--; i; i--) 
	  rangez(i) = maxfactor * Vector(y)(i) - rangez(i);
	return false; // we need a new range
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

bool HilbertReduce(Vector &z, VectorSet &S)
{
  int i;
  Vector min(z.size()), max(z.size());

  for (i = z.size(); i && !z(i); i--);
  LastNonzeroPos = i;

  do {
    // positive search
    rangez = z;
    do {
      for (i = rangez.size(); i; i--)
	if (rangez(i) >= 0) min(i) = 0, max(i) = rangez(i);
      else min(i) = rangez(i), max(i) = 0;
    } while (!S.OrthogonalRangeSearch(Leaf(&min), Leaf(&max),
				      &rangereport) 
	     && LastNonzeroPos);
    z = rangez;

    DidInvert = false; // if we inverted while in negative search, we
    // must visit positive search again. FIXME: If so, we need not go
    // to negative search again unless another inversion while in
    // positive search takes place.
    if (LastNonzeroPos) { // negative search
      rangez = -z;
      do {
	for (i = rangez.size(); i && rangez(i) == 0; i--) // zero trailing zeros
	  min(i) = max(i) = 0;
	min(i) = max(i) = 0; // zero last positive component
	for (i--; i; i--) // handle all other components
	  if (rangez(i) >= 0) min(i) = 0, max(i) = rangez(i);
	  else min(i) = rangez(i), max(i) = 0;
      } while (!S.OrthogonalRangeSearch(Leaf(&min), Leaf(&max),
					&rangereport)
	       && LastNonzeroPos);
      z = -rangez;
    }
  } while (DidInvert);

  return LastNonzeroPos;
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
  if (HilbertReduce(v, P)) {
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
      P.insert(v); Pnew.insert(v); reportx(v);
      v(p)++, v(n+1-p)++;
    }
  }

  // (3) Build all other primitive identities with exactly (t+1)
  // components of (n+1).
  for (int t = 0; t<n; t++) {
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
	  if (!v(n+1)) {
	    // if there is no (n+1) component yet, 
	    // we can invert the identity.
	    if (v(j) <= 0 || v(k) <= 0) { // otherwise, w reducible by v
	      Vector w = -v;
	      w(n+1)++, w(j)--, w(k)--;
	      ReduceAndInsert(w, P, Pnew);
	    }
	  }
	}
      }
    }
    Pold = Pnew;
    Pnew = SimpleVectorSet();
  }

  // This check should be obsolete but.

  for (VectorSet::iterator i = P.begin(); i!=P.end(); ++i) {
    RangeException = &((Vector&)(*i));
    Vector v = *i;
    if (!HilbertReduce(v, P)) 
      cerr << "Haha... a reducible vector: " << *i << " -> " << v << endl;
  }    
  RangeException = 0;

  return P;
}

int main(int argc, char *argv[])
{
  // PPI n = 5.
  int n = 0;
  if (argc >= 2) {
    sscanf(argv[1], "%d", &n);
    if (argc == 3)       
      //sscanf(argv[2], "%f", &alpha);
      sscanf(argv[2], "%d", &BBTree::FewLeavesBound);
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
	 << " PPI up to sign." << endl;
  }

}

/* $Log$
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
