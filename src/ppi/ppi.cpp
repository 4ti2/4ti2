#include <stdio.h>
#include <bool.h>
#include <set>
#include <vector>
#include <iostream.h>
#include <iomanip.h>

typedef vector<int> Vector;
typedef set<Vector, less<Vector> > VectorSet;

ostream &operator<<(ostream &s, const Vector &z)
{
  for (int i = 0; i<z.size(); i++)
    cout << setw(4) << z[i];
  return cout;
}

int HilbertDivide(Vector z, Vector y, Vector a)
{
  // Find maximal integer f with f*y<=z in Hilbert-base sense.
  int az = 0;
  for (int i = 0; i<a.size(); i++) az += a[i] * z[i];
  int ay = 0;
  int maxfactor = INT_MAX;
  for (int i = 0; i<a.size(); i++) {
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
    ay += a[i] * y[i];
  }
  if (ay > 0) {
    if (ay > az) return 0;
    // here is az>=ay>0.
    maxfactor = maxfactor <? (az / ay);
  }
  else if (ay < 0) {
    if (ay < az) return 0;
    // here is az<=ay<0.
    maxfactor = maxfactor <? (az / ay);
  }
  return maxfactor;
}

void writeppi(ostream &c, Vector z)
{
  bool first = true;
  for (int i = 0; i<z.size(); i++) 
    if (z[i] > 0) {
      if (first) first = false;
      else c << " + ";
      for (int j = 1; j < z[i]; j++)
	c << i+1 << " + ";
      c << i+1;
    }
  c << "\t= ";
  first = true;
  for (int i = 0; i<z.size(); i++) 
    if (z[i] < 0) {
      if (first) first = false;
      else c << " + ";
      for (int j = 1; j < -z[i]; j++)
	c << i+1 << " + ";
      c << i+1;
    }
  c << endl;
}

void report(Vector z, Vector a)
{
  int az = 0;
  for (int i = 0; i<a.size(); i++) az += z[i]*a[i];
  if (az) cout << z << "\t" << az << endl;
  else { cout << z << "\t\t"; writeppi(cout, z); }
}    

VectorSet HilbertBase(VectorSet T, Vector a, VectorSet freshT)
{
  VectorSet oldT, oldFreshT;
  while (oldT.size() != T.size()) {
    oldT = T;
    swap(oldFreshT, freshT);
    freshT = VectorSet();
    VectorSet::iterator iv, iw;
    bool vFresh;
    for (iv = oldT.begin(); iv != oldT.end(); ++iv) 
      for (iw = (((vFresh = oldFreshT.find(*iv) != oldFreshT.end())) 
		 ? oldFreshT.upper_bound(*iv) 
		 : oldFreshT.begin()); 
	   // if v is fresh, take only fresh w lexi-greater than v.
	   // otherwise, take all fresh w.
	   iw != oldFreshT.end(); ++iw) {
	Vector z(a.size());
	for (int i = 0; i<a.size(); i++)
	  z[i] = (*iv)[i] + (*iw)[i];
	VectorSet::iterator iy;
	bool nonzeroz = false;
	for (int i = 0; i<a.size(); i++)
	  if (nonzeroz |= !!z[i]) break;
	if (nonzeroz) {
	  int good = 0;
	  int tried = 0;
	  // T ist lexikographisch sortiert. Betrachte ich die erste
	  // Nichtnull-Koordinate von z, so kann ich den Bereich der
	  // zul"assigen y in T eingrenzen. "
	  int firstnonzero;
	  for (firstnonzero = 0; !z[firstnonzero]; firstnonzero++);
	  VectorSet::iterator begin, end;
	  if (z[firstnonzero] > 0) {
	    Vector v(a.size());
	    int i;
	    for (i = 0; i<=firstnonzero; i++) v[i] = 0;
	    if (firstnonzero+1 < a.size()) {
	      if (z[firstnonzero+1] >= 0) v[firstnonzero+1] = 0;
	      else v[firstnonzero+1] = z[firstnonzero+1];
	      for (i = firstnonzero+2; i<a.size(); i++) 
		v[i] = -INT_MAX;
	    }
	    begin = T.lower_bound(v);
	    v[firstnonzero] = z[firstnonzero];
	    if (firstnonzero+1 < a.size()) {
	      if (z[firstnonzero+1] < 0) v[firstnonzero+1] = 0;
	      else v[firstnonzero+1] = z[firstnonzero+1];
	      for (i = firstnonzero+2; i<a.size(); i++) 
		v[i] = INT_MAX;
	    }
	    end = T.upper_bound(v);
	  }
	  else {
	    Vector v(a.size());
	    int i;
	    for (i = 0; i<=firstnonzero; i++) v[i] = 0;
	    if (firstnonzero+1 < a.size()) {
	      if (z[firstnonzero+1] < 0) v[firstnonzero+1] = 0;
	      else v[firstnonzero+1] = z[firstnonzero+1];
	      for (i = firstnonzero+2; i<a.size(); i++) 
		v[i] = INT_MAX;
	    }
	    end = T.upper_bound(v);
	    v[firstnonzero] = z[firstnonzero];
	    if (firstnonzero+1 < a.size()) {
	      if (z[firstnonzero+1] >= 0) v[firstnonzero+1] = 0;
	      else v[firstnonzero+1] = z[firstnonzero+1];
	      for (i = firstnonzero+2; i<a.size(); i++) 
		v[i] = -INT_MAX;
	    }
	    begin = T.lower_bound(v);
	  }
	  // Reduktion von z durch y aus T
	  for (iy = begin; iy != end; ++iy) {
	    Vector y = *iy;
	    tried++;
	    int maxfactor = HilbertDivide(z, y, a);
	    if (maxfactor) {
	      good++;
	      nonzeroz = false;
	      for (int i = 0; i<a.size(); i++) 
		nonzeroz |= !!(z[i] -= maxfactor * y[i]);
	      if (!nonzeroz) break;
	    }
	  }
// 	  cout << "Good " << good << " of " << tried << " size " <<
// 	    T.size() << endl;
	  if (nonzeroz) {
	    T.insert(z);
	    freshT.insert(z);
	    report(z, a);
	  }
	}
      }
  }
  return T;
}

VectorSet HilbertBase(VectorSet T, Vector a)
{
  return HilbertBase(T, a, T);
}

#if 0
int main(int argc, char *argv[])
{
  /* inductive version; might be used for a cacheing solution. */
  /* this takes approx the same time as the normal version */
  int n = 0;
  if (argc == 2) sscanf(argv[1], "%d", &n);
  if (!n) n = 5;
  Vector a(n);
  for (int i = 0; i<n; i++) a[i] = i + 1;
  VectorSet H;
  for (int i = 0; i<n; i++) {
    cerr << "Inductive step " << i << endl;
    VectorSet Hfresh;
    Vector e(n);
    for (int j = 0; j<n; j++) e[j] = 0;
    e[i] = 1;
    H.insert(e), Hfresh.insert(e);
    report(e, a);
    e[i] = -1;
    H.insert(e), Hfresh.insert(e);
    report(e, a);
    H = HilbertBase(H, a, Hfresh);
  }
}
#endif

#if 1
int main(int argc, char *argv[])
{
  // PPI n = 5.
  int n = 0;
  if (argc == 2) sscanf(argv[1], "%d", &n);
  if (!n) n = 5;
  Vector a(n);
  for (int i = 0; i<n; i++) a[i] = i + 1;
  VectorSet H0;
  for (int i = 0; i<n; i++) {
    Vector e(n);
    for (int j = 0; j<n; j++) e[j] = 0;
    e[i] = 1;
    H0.insert(e);
    report(e, a);
    e[i] = -1;
    H0.insert(e);
    report(e, a);
  }
  VectorSet H = HilbertBase(H0, a);
}
#endif
