// Berechnung einer Hilbertbasis; Modifikation eines Verfahrens aus
// der Vorlesung Optimierung II (Weismantel). Hier angewendet auf
// Primitive Partitionsidentit"aten (PPI); das Verfahren ist aber
// allgemein implementiert. "

// $Id$

// Verwaltung der Testvektoren mit Range-Trees (BB-alpha based). 

// Ich speichere jetzt (nach R. Urbaniak) nur noch `modulo
// Vorzeichen', dh, erster Nichtnulleintrag ist immer positiv.
// Laufzeit vorher: 7 -> user 1m50.480s, 7 SINGLE_RANGE -> user
// 1m42.930s. Laufzeit jetzt: 7 SINGLE_RANGE -> user 0m55.900s

#include <stdio.h>
#include <bool.h>
#include <set>
#include <vector>
#include <iostream.h>
#include <iomanip.h>
#include "bbalpha.h"

typedef vector<int> Vector;
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

void writeppi(ostream &c, Vector z)
{
  bool first = true;
  for (int i = 0; i<z.size()-1; i++) 
    if (z[i] > 0) {
      if (first) first = false;
      else c << " + ";
      for (int j = 1; j < z[i]; j++)
	c << i+1 << " + ";
      c << i+1;
    }
  c << "\t= ";
  first = true;
  for (int i = 0; i<z.size()-1; i++) 
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

void report(Vector z)
{
  count++;
  if (z[z.size()-1]) cout << z << endl;
  else { ppicount++; cout << z << "\t"; writeppi(cout, z); }
}    

#ifdef SINGLE_RANGE

/* rangereport parameter */
static Vector rangez;

static bool rangereport(const Leaf &y)
{
  bool nonzeroz = false;
  int maxfactor = HilbertDivide(rangez, y);
//   if (!maxfactor) cerr << "Bad reducer" << endl;
  for (int i = 0; i<rangez.size(); i++) 
    nonzeroz |= !!(rangez[i] -= maxfactor * Vector(y)[i]);
  return nonzeroz;
}

bool HilbertReduce(Vector &z, VectorSet &S)
{
  int i;
  Vector min(z.size()), max(z.size());
  
  // positive search

  for (i = 0; i<z.size(); i++)
    if (z[i] >= 0) min[i] = 0, max[i] = z[i];
    else min[i] = z[i], max[i] = 0;
  rangez = z;
  bool nonzeroz = (S.OrthogonalRangeSearch(Leaf(&min), Leaf(&max),
					   &rangereport));
  z = rangez;

  if (nonzeroz) { // negative search

    for (i = 0; i<z.size(); i++) rangez[i] = -z[i];
    for (i = 0; i<rangez.size() && rangez[i] == 0; i++) // zero leading zeros
      min[i] = max[i] = 0;
    min[i] = max[i] = 0; // zero first positive component
    for (i = 0; i<rangez.size(); i++)
      if (rangez[i] >= 0) min[i] = 0, max[i] = rangez[i];
      else min[i] = rangez[i], max[i] = 0;
    nonzeroz = (S.OrthogonalRangeSearch(Leaf(&min), Leaf(&max),
					&rangereport));
    for (i = 0; i<z.size(); i++) z[i] = -rangez[i];
  }

  return nonzeroz;
}

#else

#error This revision only supports SINGLE_RANGE

static Vector rangeresult;
static Vector rangez;
static int maxfactor;
static Vector rangemax, rangemin;

static bool rangereport(const Leaf &y)
{
  rangeresult = y;
  maxfactor = HilbertDivide(rangez, rangeresult);
  if (!maxfactor) { 
    /* This check is superfluous, but it doesn't cost much, and
       RangeSearch used to work badly; so rather keep this just in
       case it gets broken again. */
    cerr << "Bad reducer: " << rangeresult 
	 << " min " << rangemin << " max " << rangemax << endl;
    return true;
  }
  return false; /* no further */
}

bool HilbertReduce(Vector &z, VectorSet &S)
{
  bool nonzeroz;
  do {
    int i;
    Vector min(z.size()), max(z.size());
    for (i = 0; i<z.size(); i++)
      if (z[i] >= 0) min[i] = 0, max[i] = z[i];
      else min[i] = z[i], max[i] = 0;
    rangez = z, rangemin = min, rangemax = max;
    if (S.OrthogonalRangeSearch(Leaf(&min), Leaf(&max),
				&rangereport)) {
      /* no reducer found */
      return true;
    }
    else {
      nonzeroz = false;
      for (i = 0; i<z.size(); i++) 
	nonzeroz |= !!(z[i] -= maxfactor * Vector(rangeresult)[i]);
    }
  } while (nonzeroz);
  return false;
}

#endif

void /*VectorSet*/ HilbertBase(VectorSet &T, SimpleVectorSet freshT)
{
  SimpleVectorSet oldT;
  SimpleVectorSet oldFreshT;
  while (freshT.size()) {
    cerr << "New cycle: " << oldT.size() << " " << T.size() << endl;
    oldT = T;
    swap(oldFreshT, freshT);
    freshT = SimpleVectorSet();
    SimpleVectorSet::iterator iv, iw;
    bool vFresh;
    for (iv = oldT.begin(); iv != oldT.end(); ++iv) 
      for (iw = (((vFresh = oldFreshT.find(*iv) != oldFreshT.end())) 
		 ? oldFreshT.upper_bound(*iv) 
		 : oldFreshT.begin()); 
	   // if v is fresh, take only fresh w lexi-greater than v.
	   // otherwise, take all fresh w.
	   iw != oldFreshT.end(); ++iw) {
	Vector z((*iv).size());
	int i;
	// sum of v and w: first nonzero component will be positive
	for (i = 0; i < z.size(); i++)
	  z[i] = (*iv)[i] + (*iw)[i];
	if (!HilbertDivide(z, *iv)) {
	  if (HilbertReduce(z, T)) {
	    T.insert(z);
	    freshT.insert(z);
	    report(z);
	  }
	}
	// difference of v an w: first nonzero component may be negative
	for (i = 0; 
	     i < z.size() && !(z[i] = (*iv)[i] - (*iw)[i]); i++);
	if (i < z.size()) { // v != w 
	  if (z[i] < 0) { // first nonzero component was negative; change
	    z[i] = -z[i];
	    for (i++; i < z.size(); i++) 
	      z[i] = (*iw)[i] - (*iv)[i];
	  }
	  else { // first nonzero component was positive; keep
	    for (i++; i < z.size(); i++) 
	      z[i] = (*iv)[i] - (*iw)[i];
	  }
	  if (!HilbertDivide(z, *iv)) { /* FIXME: this first check can
					   be strengthened */
	    if (HilbertReduce(z, T)) {
	      T.insert(z);
	      freshT.insert(z);
	      report(z);
	    }
	  }
	}
      }
  }
  return /*T*/;
}

void /*VectorSet*/ HilbertBase(VectorSet &T)
{
  /*return*/ HilbertBase(T, T);
}

#if 1
int main(int argc, char *argv[])
{
  // PPI n = 5.
  int n = 0;
  if (argc >= 2) {
    sscanf(argv[1], "%d", &n);
    if (argc == 3)       
      sscanf(argv[2], "%f", &alpha);
  }
  if (!n) n = 5;
  Vector a(n);
  for (int i = 0; i<n; i++) a[i] = i + 1;

  VectorSet H0(n);
  for (int i = 0; i<n; i++) {
    Vector e(n + 1);
    for (int j = 0; j<n; j++) e[j] = 0;
    e[i] = 1; e[n] = a[i];
    H0.insert(e);
    report(e);
  }
  /*VectorSet H =*/ HilbertBase(H0);
  cerr << "This makes " << 2*ppicount << " PPI of " 
       << 2*count << " test vectors." << endl;
}
#endif
