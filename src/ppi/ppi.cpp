// Bestimmung der Primitiven Partitionsidentit"aten (PPI) nach [Urbaniak];
// Verwaltung der Testvektoren mit `Digital trees'.

// $Id$

#define _NOTHREADS
// STL would use slow mutexes 

#include <stdio.h>
#include <bool.h>
#include <set>
#include <vector>
#include <iostream.h>
#include <iomanip.h>
#include "vec.h"

#define TALKATIVE 1
#define HASH
#undef WITH_STATS
#undef ERASE_SOURCES_OF_IRREDUCIBLES
#define WITH_ADVANCE
#undef BACKWARD_LEVEL

#if defined(HASH) // Use a hash table for SimpleVectorSet

#include <hash_set>

struct hash<Vector> 
{
  size_t operator()(const Vector &v) const {
    // FIXME: Is this a good hash function? Ask Knuth.
    signed long h = 0;
    for (Vector::const_iterator i = v.begin(); i!=v.end(); ++i)
      h = 5*h + *i;
    return size_t(h);
  }
};

typedef hash_set<Vector> SimpleVectorSet;

#else // Use a lexicographical balanced binary tree for SimpleVectorSet

typedef set<Vector, less<Vector> > SimpleVectorSet;

#endif

Vector Vector::operator-() const 
  return v; 
{
  v = Vector(size());
  for (int i = 0; i<size(); i++) v[i] = -(*this)[i];
}

// Implementation with digital trees (FIXME: name?)

class Node {
public:
  signed char isleaf;
  unsigned char advance;
  Node(char leaf, unsigned char adv = 1) : 
    isleaf(leaf), advance(adv) {}
};

class InnerNode : public Node {
public:
  signed char Delta;
  signed char Size;
  Node* Children[1];
  InnerNode(int which) : Node(false), Delta(-which), 
    Size(1) { Children[0] = 0; }
};

class LeafNode : public Node {
public:
  Vector vec;
  LeafNode() : Node(true) {}
};  

class DigitalTree {
public:
  int Dimension;
private:
  InnerNode *root;
  InnerNode *AdvanceNodes;
  bool DoSearch(const Vector &Min, const Vector &Max, int level, 
		Node *node);
  void Destroy(Node *node);
  void DoFinish(InnerNode *node);
  Node *&Put(InnerNode *&inner, int where, Node* n);
  Node *&Get(InnerNode *inner, int where) const { 
    static Node *null = 0;
    int i = where + inner->Delta;
    if (i < 0 || i>=inner->Size) return null;
    return inner->Children[i];
  }
public:
  DigitalTree(int dimension) : 
    Dimension(dimension), 
    root(new InnerNode(0)),
    AdvanceNodes(0) {}
  ~DigitalTree() { 
    if (root) Destroy(root);
    /*if (AdvanceNodes) delete AdvanceNodes;*/ /* FIXME: really kill 'em */}
  void Finish();
  bool Insert(const Vector &v);
  bool OrthogonalRangeSearch(const Vector & min, const Vector & max)
#if defined(BACKWARD_LEVEL)
    { return DoSearch(min, max, 0, root); }
#else
    { return DoSearch(min, max, Dimension-1, root); }
#endif
};

Node *&DigitalTree::Put(InnerNode *&inner, int where, Node* n)
{
  int i = where + inner->Delta;
  int j;
  if (i < 0) {
    // Must extend the vector to the left
    inner = (InnerNode *) 
      realloc(inner, sizeof(InnerNode) 
	      + (inner->Size - i - 1) * sizeof(Node*));
    for (j = inner->Size-1; j>=0; j--) 
      inner->Children[j-i] = inner->Children[j];
    for (j = 1; j< -i; j++) inner->Children[j] = 0;
    inner->Children[0] = n;
    inner->Delta -= i;
    inner->Size -= i;
    return inner->Children[0];
  }
  else if (i >= inner->Size) {
    // Must extend the vector to the right
    inner = (InnerNode *)
      realloc(inner, sizeof(InnerNode) + i * sizeof(Node*));
    for (j = inner->Size; j<i; j++) inner->Children[j] = 0;
    inner->Size = i + 1;
    /* FALLTHRU */
  }
  inner->Children[i] = n;
  return inner->Children[i];
}

void DigitalTree::Destroy(Node *node)
{
  if (node->isleaf == true) delete (LeafNode*)node;
  else if (node->isleaf == false) {
    int count = ((InnerNode*)node)->Size;
    for (Node **i = ((InnerNode*)node)->Children;
	 count; ++i, count--)
      if (*i) Destroy(*i);
    delete (InnerNode*)node;
  }
}

bool DigitalTree::DoSearch(const Vector & Min, const Vector & Max, 
			   int level, Node *node)
{
  vector<Node*>::pointer ci;
  int mi = max(Min[level] + ((InnerNode*)node)->Delta, 0);
  int ma = min(Max[level] + ((InnerNode*)node)->Delta,
	       (int)((InnerNode*)node)->Size - 1);
  int count = ma - mi + 1;
  int advance;
  if (count<=0) return true;
  for (ci = &(((InnerNode*)node)->Children[mi]);
       count>0; ci+=advance, count-=advance) {
#if !defined(WITH_ADVANCE)
    if (*ci) {
#endif
      if ((*ci)->isleaf == true) {
	Vector::const_pointer vi, mini, maxi;
#if defined(BACKWARD_LEVEL)
	int pos = level + 1;
	int count = Dimension - level - 1;
	for (vi = &((LeafNode*)(*ci))->vec[pos],
	       mini = &Min[pos],
	       maxi = &Max[pos];
	     count; vi++, mini++, maxi++, count--)
	  if ((*vi)<(*mini) || (*vi)>(*maxi)) goto next;
#else
	int pos = level - 1;
	int count = level;
	for (vi = &((LeafNode*)(*ci))->vec[pos],
	       mini = &Min[pos],
	       maxi = &Max[pos];
	     count; vi--, mini--, maxi--, count--)
	  if ((*vi)<(*mini) || (*vi)>(*maxi)) goto next;
#endif
	return false; /*report(((LeafNode*)node)->leaf);*/
      next:
	void(0);
      }
      else if ((*ci)->isleaf == false) {
#if defined(BACKWARD_LEVEL)
	if (!DoSearch(Min, Max, level+1, *ci)) 
	  return false;
#else
	if (!DoSearch(Min, Max, level-1, *ci)) 
	  return false;
#endif
      }
      advance = (*ci)->advance;
#if !defined(WITH_ADVANCE)
    }
    else advance = 1;
#endif
  }
  return true;
}

// Calculates the `advance' values
void DigitalTree::DoFinish(InnerNode *node)
{
  Node **i, **j;
  i = node->Children + node->Size;
  for (j = node->Children + node->Size - 1; j!=node->Children-1; j--) {
    if (*j) {
      if ((*j)->isleaf == false) DoFinish((InnerNode*)(*j));
      (*j)->advance = i - j;
      i = j;
    }
    else (*j) = Get(AdvanceNodes, i - j);
  }
}

void DigitalTree::Finish()
{ 
#if defined(WITH_ADVANCE)
  // Allocate shared `advance nodes' to put where nulls live
  AdvanceNodes = new InnerNode(2*Dimension+1);
  for (int i = 1; i<=2*Dimension+1; i++)
    Put(AdvanceNodes, i, new Node(-1, i));
  // Do the actual computation
  DoFinish(root); 
#endif
}

bool DigitalTree::Insert(const Vector &v)
{
  InnerNode **n = &root;
  int level;
#if defined(BACKWARD_LEVEL)
  for (level = 0; level<Dimension; level++) 
#else
  for (level = Dimension-1; level>=0; level--) 
#endif
    {
      int pos = v[level];
      if (Get(*n, pos)) {
	if (Get(*n, pos)->isleaf) {
	  LeafNode *ol = (LeafNode*) Get(*n, pos);
	  // Build a chain of inners up to tie 
	  do {
#if defined(BACKWARD_LEVEL)
	    level++;
#else
	    level--;
#endif
	    int newpos = v[level];
	    InnerNode *i = new InnerNode(newpos);
	    n = (InnerNode **) &Put(*n, pos, i);
	    pos = newpos;
	  } while (level > 0 
		   && (ol->vec[level] == v[level]));
	  // Put the old and new leaves at the end
	  LeafNode *l = new LeafNode;
	  l->vec = v;
	  Put(*n, v[level], l);
	  Put(*n, ol->vec[level], ol);
	  // we are done
	  return true;
	}
	else { // Is inner node
	  n = (InnerNode**) &Get(*n, pos);
	  // go on
	}
      }
      else { // Is empty slot
	LeafNode *l = new LeafNode;
	l->vec = v;
	Put(*n, pos, l);
	// we are done
	return true; 
      }
    }
  return false;
}

bool isHoppi(const Vector &z)
{
  int cnt = 0;
  for (int i = 1; i<=z.size(); i++)
    cnt += z(i);
  return (!cnt);
}

class VectorSet {
  DigitalTree *tree;
  DigitalTree *hoppitree;
  VectorSet(const VectorSet &);
  VectorSet &operator=(const VectorSet &);
  
public:
  VectorSet(int level) : tree(new DigitalTree(level+1)),
    hoppitree(new DigitalTree(level+1)) {}
  ~VectorSet() { delete tree; delete hoppitree; }
  void Finish() { tree->Finish(); hoppitree->Finish(); }
  bool insert(const Vector &v) {
    if (isHoppi(v)) {
      Vector w = v;
      Vector::iterator i;
      // normalize hoppi
      for (i = w.begin(); i!=w.end() && !*i; ++i);
      int cnt = i - w.begin();
      if (cnt) {
	w.erase(w.begin(), i);
	w.insert(w.end(), cnt, Vector::value_type(0));
      }	
      return hoppitree->Insert(w);
    }
    else return tree->Insert(v); 
  }
  bool OrthogonalRangeSearch(const Vector & min, const Vector & max)
    { 
      if (!tree->OrthogonalRangeSearch(min, max)) return false;
      Vector minn = min;
      Vector maxx = max;
      int i;
      for (i = minn.size(); !minn(i) && !maxx(i); i--);
      for (; i; i--) { 
	if (!hoppitree->OrthogonalRangeSearch(minn, maxx)) return false;
	minn.erase(minn.begin()); minn.push_back(0);
	maxx.erase(maxx.begin()); maxx.push_back(0);
      }
      return true;
    }
};

ostream &operator<<(ostream &s, const Vector &z)
{
  for (int i = 0; i<z.size(); i++)
    s << setw(4) << (int)(z[i]);
  return s;
}

void ClearVectorRepository()
{
}

int HilbertDivide(Vector z, Vector y)
{
  // Find maximal integer f with f*y<=z in Hilbert-base sense.
  int maxfactor = INT_MAX;
  for (int i = 0; i<z.size(); i++) {
    if (y[i] > 0) {
      if (y[i] > z[i]) return 0;
      // here is z[i]>=y[i]>0.
      maxfactor = min(maxfactor, z[i] / y[i]);
    }
    else if (y[i] < 0) {
      if (y[i] < z[i]) return 0;
      // here is z[i]<=y[i]<0.
      maxfactor = min(maxfactor, z[i] / y[i]);
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

static int ppicount, hoppicount;
static int dupcount, redcount, hitcount[2], failcount[2];

/* rangereport parameters */
static Vector rangemin, rangemax;
static int LastNonzeroPos;

/////////////////

//
// Specialized code for generating all primitive partition identities
//

void reportx(Vector z)
{
  ppicount++;
  if (isHoppi(z)) hoppicount++;
#if (TALKATIVE>=1)
  cout << z << "\t"; writeppi(cout, z, z.size()); 
#endif
}

inline void SetupAttribute(Vector &v, int n, bool WithNegative)
{
  // Vector::Attribute is used here as follows: A bit is clear if
  // the corresponding position does not need to be touched any more
  // (Source Erasing). We take two bits per position because of
  // positive/negative stuff.
  // NOTE: This imposes the arbitrary limit n <= 33.
  unsigned long attrib = 0;
  unsigned long mask = 1;
  if (WithNegative) {
    for (int j = 1; j<=n/2; j++, mask<<=1) {
      if (v(j) > 0 || v(n+1-j) > 0) attrib |= mask;
      if (v(j) < 0 || v(n+1-j) < 0) attrib |= (mask<<16);
    }
  }
  else {
    for (int j = 1; j<=n/2; j++, mask<<=1)
      if (v(j) > 0 || v(n+1-j) > 0) attrib |= mask;
  }
  v.Attribute = attrib;
}

inline void EraseSource(const Vector &w, int n, SimpleVectorSet &Pold,
			bool irreducible, bool WithNegative)
{
#if !defined(ERASE_SOURCES_OF_IRREDUCIBLES)
  if (!irreducible) return; // they often fail, which is expensive
#endif
  Vector u = w;
  u(n+1)--;
  for (int i = 1; i<=n/2; i++) {
    if (u(i) < 0 || u(n+1-i) < 0) {
      u(i)++, u(n+1-i)++;
      SimpleVectorSet::iterator ui = Pold.find(u);
      if (ui != Pold.end()) {
#if defined(WITH_STATS)
	hitcount[irreducible]++;
#endif
	const_cast<Vector&>(*ui).Attribute &= ~(1<<(i-1));
      }
#if defined(WITH_STATS)
      else failcount[irreducible]++;
#endif
      u(i)--, u(n+1-i)--;
    } 
  }
  if (WithNegative) {
    u = -u;
    for (int i = 1; i<=n/2; i++) {
      if (u(i) > 0 || u(n+1-i) > 0) {
	u(i)--, u(n+1-i)--;
	SimpleVectorSet::iterator ui = Pold.find(u); 
	if (ui != Pold.end()) {
#if defined(WITH_STATS)
	  hitcount[irreducible]++;
#endif
	  const_cast<Vector&>(*ui).Attribute &= ~((1<<16)<<(i-1));
	}
#if defined(WITH_STATS)
	else failcount[irreducible]++;
#endif
	u(i)++, u(n+1-i)++;
      } 
    }    
  }
}

inline void RaisePPI(const Vector &v, int j, int k, int n, 
		     VectorSet &P, SimpleVectorSet &Pold,
		     SimpleVectorSet &Pnew, bool WithNegative)
{
  // assert(j<k);
  Vector w = v;
  w(n+1)++, w(j)--, w(k)--;
  // reducibility check
  if (w(j) >= 0 && w(k) >= 0) { // w is irreducible 
    EraseSource(w, n, Pold, true, WithNegative);
    SetupAttribute(w, n, false);
    if (Pnew.insert(w).second) reportx(w);
  }
  else { // may be reducible, must try to reduce 
    int i;
#if defined(WITH_STATS)
    redcount++;
#endif
    rangemin = Vector(n+1), rangemax = Vector(n+1);
    rangemin(n+1) = rangemax(n+1) = 0;
    // use full range at the remaining positions
    for (i = n; i; i--) {
      if (w(i) >= 0) rangemin(i) = 0, rangemax(i) = w(i);
      else rangemin(i) = w(i), rangemax(i) = 0;
    }
    // FIXME: maybe special-case k or j = n (only one sign is possible
    // in this case.)
    
    // at the position that increased (abs value), demand this value
    if (w(k) < 0) rangemin(k) = w(k), rangemax(k) = w(k);
    else rangemin(j) = w(j), rangemax(j) = w(j);
    // search!
    if (P.OrthogonalRangeSearch(rangemin, rangemax)) {
      // try negative search
      for (i = n; i; i--) {
	if (w(i) <= 0) rangemin(i) = 0, rangemax(i) = -w(i);
	else rangemin(i) = -w(i), rangemax(i) = 0;
      }
      if (w(k) < 0) rangemin(k) = -w(k), rangemax(k) = -w(k);
      else rangemin(j) = -w(j), rangemax(j) = -w(j);
      if (P.OrthogonalRangeSearch(rangemin, rangemax)) {
	// didn't find reducer, but vector may already be known
	SetupAttribute(w, n, false);
	if (Pnew.insert(w).second) {
	  EraseSource(w, n, Pold, true, WithNegative);
	  reportx(w);
	}
	else {
#if defined(WITH_STATS)
	  dupcount++;
#endif
#if (TALKATIVE>=2)
	  cerr << "Duplicate: " << w << endl;
#endif
	}
	return;
      }
    }
    EraseSource(w, n, Pold, false, WithNegative);
  }
}

void ExtendPPI(SimpleVectorSet &Pn, int n)
{
  VectorSet P(n);
  SimpleVectorSet *Pold = new SimpleVectorSet;
  SimpleVectorSet *Pnew = new SimpleVectorSet;

  // Implementation of Algorithms 4.3.8, 4.3.11 from [Urbaniak] 

  // (1) Create the range-searchable set P of (n+1)-vectors from the
  // set Pn of n-vectors.
  cerr << "# Vectors copied from n = " << n << ": " << endl;
  for (SimpleVectorSet::iterator i = Pn.begin(); i!=Pn.end(); ++i) {
    Vector v(n+1);
    for (int j = 1; j <= n; j++) v(j) = (*i)(j);
    v(n+1) = 0;
    SetupAttribute(v, n, true); 
    P.insert(v); Pold->insert(v); reportx(v);
    if (isHoppi(v) && v(n)!=0) {
      // perform hoppi shift
      cerr << "H";
      v.pop_back();
      v.insert(v.begin(), 0);
      SetupAttribute(v, n, false);
      Pnew->insert(v); reportx(v);
    }
  }
  cerr << "Erasing sources of shifted hoppis..." << endl;
  for (SimpleVectorSet::iterator i = Pnew->begin(); 
       i!=Pnew->end(); ++i) 
    EraseSource(*i, n, *Pold, true, false);
  // destroy Pn to save memory
  Pn = SimpleVectorSet();
  // finish P to make search faster
  cerr << "Finishing..." << endl;
  P.Finish();

  // (2) Add `(n+1) = a + b' identities.
  cerr << "# Vectors of type " << n+1 << " = a + b:" << endl;
  {
    Vector v(n+1);
    for (int j = 1; j <= n; j++) v(j) = 0;
    v(n+1) = +1;
    for (int p = 1; p<=(n+1)/2; p++) {
      // takes care of case: p = n+1-p.
      v(p)--, v(n+1-p)--;
      SetupAttribute(v, n, false);
      Pnew->insert(v); 
      reportx(v);
      v(p)++, v(n+1-p)++;
    }
  }

  // (3) Build all other primitive identities with exactly one
  // component of (n+1).
  cerr << "# Vectors of P" << 1 << "(" << n+1 << "):" << endl;
  for (SimpleVectorSet::iterator i = Pold->begin(); i!=Pold->end(); ++i) {
    Vector v = *i;
    if (v.Attribute) {
#if (TALKATIVE>=2)
      cerr << "Raising " << v << endl;
#endif
      for (int j = 1; j<=n/2; j++) { // yes n/2 is enough
	int k = (n+1) - j;
	if (v.Attribute & (1<<(j-1)))
	  RaisePPI(v, j, k, n, P, *Pold, *Pnew, true);
	if (v.Attribute & ((1<<16)<<(j-1)))
	  RaisePPI(-v, j, k, n, P, *Pold, *Pnew, true);
      }
    }
  }
  //
  Pn.insert(Pold->begin(), Pold->end());
  delete Pold;
  Pold = Pnew;
  Pnew = new SimpleVectorSet;

  // (4) Build all other primitive identities with exactly (t+1)
  // components of (n+1).
  for (int t = 1; t<n; t++) {
    cerr << "# Vectors of P" << t+1 << "(" << n+1 << "):" << endl;
    for (SimpleVectorSet::iterator i = Pold->begin(); i!=Pold->end(); ++i) {
      Vector v = *i;
      if (v.Attribute) {
#if (TALKATIVE>=2)
	cerr << "Raising " << v << endl;
#endif
	for (int j = 1; j<=n/2; j++) {
	  if (v.Attribute & (1<<(j-1))) { // source erasing
	    int k = (n+1) - j;
	    RaisePPI(v, j, k, n, P, *Pold, *Pnew, false);
	  }
	}
      }
    }
    Pn.insert(Pold->begin(), Pold->end());
    delete Pold;
    Pold = Pnew;
    Pnew = new SimpleVectorSet;
  }

  Pn.insert(Pold->begin(), Pold->end());
  delete Pold;
  delete Pnew;
}

int main(int argc, char *argv[])
{
  // PPI n = 5.
  int n = 0;
  if (argc >= 2) {
    sscanf(argv[1], "%d", &n);
  }
  if (!n) n = 5;

  // Setup PPI set for n=2
  SimpleVectorSet V;
  Vector v(2); v(1) = -2, v(2) = +1; V.insert(v);
  for (int i = 2; i<n; i++) {
#if defined(WITH_STATS)
    hitcount[0] = failcount[0] = hitcount[1] = failcount[1] 
      = redcount = dupcount = 0;
#endif
    hoppicount = ppicount = 0;
    cerr << "### Extending to n = " << i+1 << endl;
    ExtendPPI(V, i);
    cerr << "### This makes " << ppicount 
	 << " PPIs up to sign, including " 
	 << hoppicount << " HoPPIs"
#if defined(WITH_STATS)
	 << ", with " << redcount << " reduce ops, "
	 << hitcount[0] << " red-erase hits, "
	 << failcount[0] << " red-erase fails, "
	 << hitcount[1] << " irr-erase hits, "
	 << failcount[1] << " irr-erase fails, and "
	 << dupcount << " duplicates. " 
#endif
	 << endl;
    cerr << "### Clearing vector repository..." << flush;
    ClearVectorRepository();
    cerr << "done." << endl;
  }

}

/* $Log$
 * Revision 1.28.1.5  1999/03/12 13:55:19  mkoeppe
 * Some clean-up with the vectors. BACKWARD_LEVEL option, but no impact
 * on performance.
 *
 * Revision 1.28.1.4  1999/03/12 11:57:09  mkoeppe
 * InnerNodes take less memory (no longer using STL vectors). n=19 takes
 * 57M to 61M, 7m43s user time.
 *
 * Revision 1.28.1.3  1999/03/11 18:59:54  mkoeppe
 * Uses a hash table for SimpleVectorSet (strong). Uses advance pointers
 * with digital trees (weak). Defines _NOTHREADS (strong). n=19 takes 8m36s.
 *
 * Revision 1.28.1.2  1999/03/11 16:27:28  mkoeppe
 * Avoids duplicates by `erasing sources'. Only slightly faster because
 * the lexicographical set<Vector> operations take lotsa time.
 *
 * Revision 1.28.1.1  1999/03/11 11:27:55  mkoeppe
 * Case `j=k=\frac{n+1}2' taken out because proven obsolete. Taraaa code
 * taken out because proven obsolete. Counts duplicates.
 *
 * Revision 1.28  1999/03/10 17:44:15  mkoeppe
 * The vectors in inner node of digital trees now grow on
 * demand. Insertion takes a bit longer, but the whole thing is much
 * cooler to memory, and lookups are also faster. n=18 takes 5m33s, n=19
 * takes 11m8s user time at 61MB memory use.
 *
 * Revision 1.27  1999/03/09 20:54:59  mkoeppe
 * Clean up.
 *
 * Revision 1.26  1999/03/09 20:33:37  mkoeppe
 * Only maintains P0 range-searchable, which saves memory. n=18 now takes
 * user 6m20.4s.
 *
 * Revision 1.25  1999/03/09 20:10:44  mkoeppe
 * n=18 took user 9m34.410s, but was at very memory edge.
 *
 * Revision 1.24  1999/03/09 17:36:01  mkoeppe
 * Clean up.
 *
 * Revision 1.23  1999/03/09 16:46:07  mkoeppe
 * Too cool to be true. Up to n=17 everything is computed in `zero
 * time'. From n=18 the memory load is too high.
 *
 * Revision 1.22  1999/03/09 15:12:56  mkoeppe
 * *** empty log message ***
 *
 * Revision 1.21  1999/03/09 00:53:49  mkoeppe
 * Using a simple `digital tree' instead. n=11 takes 3sec.
 *
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
