// Berechnung einer Hilbertbasis; Modifikation eines Verfahrens aus
// der Vorlesung Optimierung II (Weismantel). Hier angewendet auf
// Primitive Partitionsidentit"aten (PPI); das Verfahren ist aber
// allgemein implementiert. "

// Verwaltung der Testvektoren mit Range-Trees (BB-alpha based). 
// Laufzeiten: n=5 3.5s, n=6 26s, n=7 4m43s, n=8 36m, n=9 308m, n=10 1900m

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
  for (i = 0; i<z.size(); i++)
    if (z[i] >= 0) min[i] = 0, max[i] = z[i];
    else min[i] = z[i], max[i] = 0;
  rangez = z;
  bool nonzeroz = (S.OrthogonalRangeSearch(Leaf(&min), Leaf(&max),
					   &rangereport));
  z = rangez;
  return nonzeroz;
}

#else

static Vector rangeresult;
static Vector rangez;
static int maxfactor;
static Vector rangemax, rangemin;

static bool rangereport(const Leaf &y)
{
  rangeresult = y;
  maxfactor = HilbertDivide(rangez, rangeresult);
  if (!maxfactor) { /* FIXME: This check should should be superfluous, but 
		       RangeSearch does not work properly */
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
	bool nonzeroz = false;
	for (int i = 0; i < z.size(); i++)
	  nonzeroz |= !!(z[i] = (*iv)[i] + (*iw)[i]);
	if (nonzeroz && !HilbertDivide(z, *iv)) {
	  if (HilbertReduce(z, T)) {
	    T.insert(z);
	    freshT.insert(z);
	    report(z);
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
    e[i] = -1; e[n] = -a[i];
    H0.insert(e);
    report(e);
  }
  /*VectorSet H =*/ HilbertBase(H0);
  cerr << "This makes " << ppicount << " PPI of " 
       << count << " test vectors." << endl;
}
#endif
