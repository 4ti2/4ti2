// Bestimmung der Primitiven Partitionsidentit"aten (PPI) nach [Urbaniak];
// Verwaltung der Testvektoren mit `Digital trees'.

// $Id$

#define _NOTHREADS
// STL would use slow mutexes 

#include <stdio.h>
#include <set>
#include <vector>
#include <iostream.h>
#include <iomanip.h>
#include <sys/times.h>
#include <time.h>
#include <limits.h>
#include <functional>

using namespace std;

double user_time()
{
  struct tms t;
  times(&t);
  return (double)t.tms_utime / CLK_TCK;
}

#define TALKATIVE 0
#define OUR_OWN_HASH
//#define HASH
#undef WITH_STATS
#undef ERASE_SOURCES_OF_IRREDUCIBLES
#define ERASE_ONLY_GUARANTEED
#define WITH_ADVANCE
#undef BACKWARD_LEVEL
#define COMPACT_VECTORS

#if defined(COMPACT_VECTORS)
struct VectorAux {
  /* order important */
  unsigned long Attribute;
  signed char Stuff[4];
  void *operator new(size_t s, int length) {
    if (length < 3) length=3;
    return malloc(s - 3 +length);
  }
  VectorAux(int length) { Length() = length; }
  VectorAux(const VectorAux &aux) {
    memcpy(this, &aux, sizeof(VectorAux) -4 + 1 + aux.Length());
  }
  inline unsigned char &Length() { return (unsigned char &) Stuff[0]; }
  inline int Length() const { return (unsigned char) Stuff[0]; }
  inline signed char *Data() { return Stuff+1; }
};

class Vector {
public:
  typedef const signed char *const_pointer;
  typedef signed char *pointer;
  typedef const signed char *const_iterator;
  typedef signed char *iterator;
  friend class LeafNode;
private:
  VectorAux *aux;
public:
  Vector() { aux = 0; }
  Vector(int length) { aux = new(length) VectorAux(length); }
  Vector(const Vector &v) { aux = new(v.aux->Length()) VectorAux(*v.aux); }
  ~Vector() { delete aux; }
  Vector &operator=(const Vector &v) {
    if (this!=&v) { 
      delete aux; 
      aux = new(v.aux->Length()) VectorAux(*v.aux);
    }
    return *this;
  }
  signed char &operator [](size_t i) { return aux->Data()[i]; }
  signed char &operator ()(int i) { return aux->Data()[i-1]; }
  const signed char &operator [](size_t i) const { return aux->Data()[i]; }
  const signed char &operator ()(int i) const { return aux->Data()[i-1];
  }
  const_iterator begin() const { return aux->Data(); }
  const_iterator end() const { return aux->Data() + aux->Length(); }
  size_t size() const { return aux->Length(); }
public:
  Vector operator-() const;
  unsigned long &Attribute() { return aux->Attribute; }
  friend bool operator==(const Vector &a, const Vector &b) {
    if (a.size() != b.size()) return false;
    Vector::const_iterator ai, bi;
    for (ai = a.begin(), bi = b.begin(); 
	 ai!=a.end(); ++ai, ++bi) if (*ai != *bi) return false;
    return true;
  }
public:    /* for use in SimpleVectorSet */
  size_t hash;
  Vector *next;
};  
#else
#  include "vec.h"
#endif

#if defined(OUR_OWN_HASH)

class SimpleVectorSet;

/* a simple hash table with chaining */

size_t vector_hash_fun (const Vector &v)
{
  signed long h = 0;
  for (Vector::const_iterator i = v.begin(); i!=v.end(); ++i)
    h = 5*h + *i;
  return size_t(h);
}

class SimpleVectorSetIterator {
public:
  size_t index;
  SimpleVectorSet *the_set;
  Vector *node;
  SimpleVectorSetIterator(SimpleVectorSet *theset, size_t ndex, Vector *ode) :
    the_set(theset), index(ndex), node(ode) {}
  Vector &operator*();
  friend bool operator == (const SimpleVectorSetIterator &a, const SimpleVectorSetIterator &b);
  friend bool operator != (const SimpleVectorSetIterator &a, const SimpleVectorSetIterator &b);
  SimpleVectorSetIterator &operator ++();
};

inline bool operator == (const SimpleVectorSetIterator &a, const SimpleVectorSetIterator &b)
{
  return a.index == b.index && a.node == b.node;
}
  
inline bool operator != (const SimpleVectorSetIterator &a, const SimpleVectorSetIterator &b)
{
  return a.index != b.index || a.node != b.node;
}

class SimpleVectorSet {
  SimpleVectorSet (const SimpleVectorSet &);
  SimpleVectorSet &operator=(const SimpleVectorSet &);
public:
  typedef SimpleVectorSetIterator iterator;
  Vector **buckets;
  size_t num_buckets;
  size_t count;
  SimpleVectorSet(size_t size) {
    count = 0;
    num_buckets = size; 
    buckets = (Vector**) malloc(sizeof(Vector *) * num_buckets);
    size_t i;
    for (i = 0; i<num_buckets; i++)
      buckets[i] = NULL;
  }
  ~SimpleVectorSet() {
    free(buckets);
  }
  SimpleVectorSetIterator begin() {
    size_t bucket;
    for (bucket = 0; bucket < num_buckets; bucket++) {
      if (buckets[bucket] != NULL)
	return SimpleVectorSetIterator(this, bucket, buckets[bucket]);
    }
    return end();
  }
  SimpleVectorSetIterator end() { return SimpleVectorSetIterator(this, num_buckets, NULL); }
  pair<iterator, bool> insert(const Vector& obj)
  {
    size_t h = vector_hash_fun(obj);
    size_t hm = h % num_buckets;
    /* See if it's there */
    Vector *v;
    for (v = buckets[hm]; v!=NULL; v=v->next) {
      if (v->hash == h) {
	if (*v == obj) {
	  return pair<iterator, bool>(SimpleVectorSetIterator(this, hm, v), false);
	}
      }
    }
    /* Insert it at the front of the list */
    v = new Vector(obj);
    v->next = buckets[hm];
    v->hash = h;
    buckets[hm] = v; count++;
    return pair<iterator, bool>(SimpleVectorSetIterator(this, hm, v), true);
  }
  void insert(const SimpleVectorSetIterator &a, const SimpleVectorSetIterator &b)
  {
    iterator i = a;
    for (; i!=b; ++i) {
      insert(*i);
    }
  }
  iterator find(const Vector& obj)
  {
    size_t h = vector_hash_fun(obj);
    size_t hm = h % num_buckets;
    Vector *v;
    for (v = buckets[hm]; v!=NULL; v=v->next) {
      if (v->hash == h) {
	if (*v == obj) {
	  return SimpleVectorSetIterator(this, hm, v);
	}
      }
    }
    return end();
  }
};

inline Vector &SimpleVectorSetIterator::operator*() {
  return *node;
}

inline SimpleVectorSetIterator &SimpleVectorSetIterator::operator ++()
{
  node = node->next;
  if (node == NULL) {
    /* end of list, next bucket */
    do {
      index++;
    } while (index < the_set->num_buckets && the_set->buckets[index] == NULL);
    if (index < the_set->num_buckets)
      node = the_set->buckets[index];
  }    
  return *this;
}

////// FIXME:

class VerySimpleVectorSet;

#else
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
#endif

Vector Vector::operator-() const 
{
  Vector v(size());
  for (int i = 0; i<size(); i++) v[i] = -(*this)[i];
  return v;
}

// Implementation with digital trees (FIXME: name?)

class Node {
public:
  unsigned char advance;
  signed char isleaf;
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
  const signed char *vec;
  LeafNode() : Node(true) {}
  LeafNode(const Vector &v, int length) : Node(true) {
    vec = v.aux->Data();
  }
};  

class DigitalTree {
  InnerNode *root;
  InnerNode *AdvanceNodes;
  bool DoSearch(const Vector &min, const Vector &max, int level, 
		Node *node);
  void Destroy(Node *node);
  void DestructiveStore(Node *node, SimpleVectorSet &S);
  void DoFinish(InnerNode *node);
  Node *&Put(InnerNode *&inner, int where, Node* n);
  Node *&Get(InnerNode *inner, int where) const { 
    static Node *null = 0;
    int i = where + inner->Delta;
    if (i < 0 || i>=inner->Size) return null;
    return inner->Children[i];
  }
public:
  int Dimension;
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
  void DestructiveCopy(SimpleVectorSet &S)
    { DestructiveStore(root, S); root = 0; }
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

bool DigitalTree::DoSearch(const Vector & min, const Vector & max, 
			   int level, Node *node)
{
  vector<Node*>::pointer ci;
  int mi = (min[level] + ((InnerNode*)node)->Delta);
  if (mi < 0) mi = 0;
  int ma = (max[level] + ((InnerNode*)node)->Delta);
  if (ma > ((int)((InnerNode*)node)->Size - 1))
    ma = ((int)((InnerNode*)node)->Size - 1);
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
	       mini = &min[pos],
	       maxi = &max[pos];
	     count; vi++, mini++, maxi++, count--)
	  if ((*vi)<(*mini) || (*vi)>(*maxi)) goto next;
#else
	int pos = level - 1;
	int count = level;
	for (vi = &((LeafNode*)(*ci))->vec[pos],
	       mini = &min[pos],
	       maxi = &max[pos];
	     count; vi--, mini--, maxi--, count--)
	  if ((*vi)<(*mini) || (*vi)>(*maxi)) goto next;
#endif
	return false; /*report(((LeafNode*)node)->leaf);*/
      next:
	void(0);
      }
      else if ((*ci)->isleaf == false) {
#if defined(BACKWARD_LEVEL)
	if (!DoSearch(min, max, level+1, *ci)) 
	  return false;
#else
	if (!DoSearch(min, max, level-1, *ci)) 
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
  for (level = 0; level<Dimension; level++) {
#else
  for (level = Dimension-1; level>=0; level--) {
#endif
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
	LeafNode *l = new LeafNode(v, level);
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
      LeafNode *l = new LeafNode(v, level);
      Put(*n, pos, l);
      // we are done
      return true; 
    }
  }
}

class VectorSet {
  DigitalTree *tree;
  VectorSet(const VectorSet &);
  VectorSet &operator=(const VectorSet &);
  
public:
  VectorSet(int level) : tree(new DigitalTree(level+1)) {}
  ~VectorSet() { delete tree; }
  void Finish() { tree->Finish(); }
  bool insert(const Vector &v) { 
    return tree->Insert(v); 
  }
  bool OrthogonalRangeSearch(const Vector & min, const Vector & max)
    { return tree->OrthogonalRangeSearch(min, max); }
  void DestructiveCopy(SimpleVectorSet &S)
    { tree->DestructiveCopy(S); }
};

ostream &operator<<(ostream &s, const Vector &z)
{
  for (int i = 0; i<z.size(); i++)
    s << setw(4) << (int)(z[i]);
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
      int f = (z[i] / y[i]);
      if (f < maxfactor)
	maxfactor = f;
    }
    else if (y[i] < 0) {
      if (y[i] < z[i]) return 0;
      // here is z[i]<=y[i]<0.
      int f = (z[i] / y[i]);
      if (f < maxfactor)
	maxfactor = f;
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

static int ppicount;
static int dupcount, simplecount, goodcount, redcount, hitcount[2], failcount[2];

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
  v.Attribute() = attrib;
}

// Erase u as a source, i.e. find it in Pold and erase its i'th flag.
inline void DoErase(Vector &u, int i, int n, SimpleVectorSet &Pold,
		    bool irreducible, bool negative)
{
  SimpleVectorSet::iterator ui = Pold.find(u);
  if (ui != Pold.end()) {
#if defined(WITH_STATS)
    hitcount[irreducible]++;
#endif
    const_cast<Vector&>(*ui).Attribute() 
      &= ~((negative ? (1<<16) : 1)<<(i-1));
  }
#if defined(WITH_STATS)
  else {
    failcount[irreducible]++;
    //    cerr << "Failed: " << u << endl;
  }
#endif
}

// Increase u at i and its complement if this leads to a valid source,
// and erase this source.
inline void DoIncAndErase(Vector &u, int i, int n, SimpleVectorSet &Pold,
			  bool irreducible)
{
#if defined(ERASE_ONLY_GUARANTEED)
  if (u(i) < 0 || u(n+1-i) < 0) // i in supp(u-)
#endif
    { 
      u(i)++, u(n+1-i)++;
      DoErase(u, i, n, Pold, irreducible, false);
      u(i)--, u(n+1-i)--;
    } 
}

inline void DoDecAndErase(Vector &u, int i, int n, SimpleVectorSet &Pold,
			  bool irreducible)
{
#if defined(ERASE_ONLY_GUARANTEED)
  if (u(i) > 0 || u(n+1-i) > 0)
#endif
    {
      u(i)--, u(n+1-i)--;
      DoErase(u, i, n, Pold, irreducible, true);
      u(i)++, u(n+1-i)++;
    } 
}

inline void EraseSource(const Vector &w, int n, SimpleVectorSet &Pold,
			bool irreducible, bool WithNegative)
{
#if !defined(ERASE_SOURCES_OF_IRREDUCIBLES)
  if (!irreducible) return; // they often fail, which is expensive
#endif
  //  cerr << "EraseSource " << w << endl; 
  Vector u = w;
  u(n+1)--;
  int LastNonzeroPos;
  if (WithNegative) { 
    // find out last nonzero pos
    // asserted is w(n+1)==1, i.e. u(n+1)==0
    for (LastNonzeroPos = n; !u(LastNonzeroPos); LastNonzeroPos--);
    if (u(LastNonzeroPos)<0 && LastNonzeroPos>n/2) { 
      // This is a complicated case.
      int i;
      for (i = 1; i<n+1-LastNonzeroPos; i++) {
	// We will put a positive value to the right of
	// LastNonzeroPos, so perform positive search.
	DoIncAndErase(u, i, n, Pold, irreducible);
      }
      Vector mu;
      if (2*LastNonzeroPos == n+1) { // Very special case:
	// Need not erase sources because this case is left out
	// anyway. (The code below doesn't grok this case.)
	mu = -u;
      }
      else {
	// Handle LastNonzeroPos:
	u(i)++, u(LastNonzeroPos)++; // This might zero it, so:
	mu = -u;
	for (; LastNonzeroPos && !u(LastNonzeroPos); LastNonzeroPos--);
	if (LastNonzeroPos) {
	  if (u(LastNonzeroPos) > 0) {
	    // Perform positive search
	    DoErase(u, i, n, Pold, irreducible, false);
	  }
	  else {
	    // Perform negative search
	    DoErase(mu, i, n, Pold, irreducible, true);
	  }
	}
	mu(i)++, mu(n+1-i)++; 
      }
      // (We have mu == -original_u here.)
      // Now for the rest:
      for (i++; i<=n/2; i++) {
	// LastNonzeroPos is left invariant here, so perform negative
	// search. 
	DoDecAndErase(mu, i, n, Pold, irreducible);
      }
      return;
    }
  }
  // we only have to do positive search!
  for (int i = 1; i<=n/2; i++) 
    DoIncAndErase(u, i, n, Pold, irreducible);
}

inline void RaisePPI(const Vector &v, int j, int k, int n, 
		     VectorSet &P, SimpleVectorSet &Pold,
		     SimpleVectorSet &Pnew, bool WithNegative,
		     bool KnownIrreducible = false)
{
  // assert(j<k);
  //  cerr << "RaisePPI " << v << " with " << j << " and " << k << endl;
  Vector w = v;
  w(n+1)++, w(j)--, w(k)--;
  // reducibility check
  if (KnownIrreducible || (w(j) >= 0 && w(k) >= 0)) { // w is irreducible 
#if defined(WITH_STATS)
    simplecount++;
#endif
    SetupAttribute(w, n, false);
    if (Pnew.insert(w).second) {
      EraseSource(w, n, Pold, true, WithNegative);
      reportx(w); 
    }
    else {
#if defined(WITH_STATS)
      dupcount++;
#endif
    }
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
#if defined(WITH_STATS)
	  goodcount++;
#endif
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

/* FIXME: Pn kann einfach ein Array werden, da darin nie gesucht
   werden mu"s. */
SimpleVectorSet *ExtendPPI(SimpleVectorSet *Pn, int n)
{
  SimpleVectorSet *Pold;
  int expected_count = 2 * Pn->count;
  if (expected_count < 20000) expected_count = 20000;
  SimpleVectorSet *Pnew = new SimpleVectorSet(expected_count);
  SimpleVectorSet *Pbase = new SimpleVectorSet(Pn->count);

  // Implementation of Algorithms 4.3.8, 4.3.11 from [Urbaniak] 
  
  {
    VectorSet P(n);

    // (1) Create the range-searchable set P of (n+1)-vectors from the
    // set Pn of n-vectors.
    cerr << "# Vectors copied from n = " << n << ": " << endl;
    for (SimpleVectorSet::iterator i = Pn->begin(); i!=Pn->end(); ++i) {
      Vector v(n+1);
      for (int j = 1; j <= n; j++) v(j) = (*i)(j);
      v(n+1) = 0;
      SetupAttribute(v, n, true); 
      P.insert(*Pbase->insert(v).first); reportx(v);
    }
    // destroy Pn to save memory
    delete Pn;
    Pn = new SimpleVectorSet(expected_count);
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

    // first a pass for the simple cases
    for (SimpleVectorSet::iterator i = Pbase->begin(); i!=Pbase->end(); ++i) {
      Vector v = *i;
      if (v.Attribute()) {
#if (TALKATIVE>=2)
	cerr << "Raising " << v << endl;
#endif
	for (int j = 1; j<=n/2; j++) { // yes n/2 is enough
	  int k = (n+1) - j;
	  if (v(j) > 0 && v(k) > 0)
	    RaisePPI(v, j, k, n, P, *Pbase, *Pnew, true, true);
	  if (v(j) < 0 && v(k) < 0)
	    RaisePPI(-v, j, k, n, P, *Pbase, *Pnew, true, true);
	}
      }
    }

    // then a pass for the complicated cases
    for (SimpleVectorSet::iterator i = Pbase->begin(); i!=Pbase->end(); ++i) {
      Vector v = *i;
      if (v.Attribute()) {
#if (TALKATIVE>=2)
	cerr << "Raising " << v << endl;
#endif
	for (int j = 1; j<=n/2; j++) { // yes n/2 is enough
	  int k = (n+1) - j;
	  if (v.Attribute() & (1<<(j-1)))
	    RaisePPI(v, j, k, n, P, *Pbase, *Pnew, true);
	  if (v.Attribute() & ((1<<16)<<(j-1)))
	    RaisePPI(-v, j, k, n, P, *Pbase, *Pnew, true);
	}
      }
    }
    //
    Pold = Pnew;
    Pnew = new SimpleVectorSet(expected_count);

    // (4) Build all other primitive identities with exactly (t+1)
    // components of (n+1).

    for (int t = 1; t<n; t++) {
      cerr << "# Vectors of P" << t+1 << "(" << n+1 << "):" << endl;
      // simple cases
      for (SimpleVectorSet::iterator i = Pold->begin(); i!=Pold->end(); ++i) {
	Vector v = *i;
	if (v.Attribute()) {
#if (TALKATIVE>=2)
	  cerr << "Raising " << v << endl;
#endif
	  for (int j = 1; j<=n/2; j++) {
	    int k = (n+1) - j;
	    if (v(j) > 0 && v(k) > 0) { // source erasing
	      RaisePPI(v, j, k, n, P, *Pold, *Pnew, false, true);
	    }
	  }
	}
      }
      // complicated cases
      for (SimpleVectorSet::iterator i = Pold->begin(); i!=Pold->end(); ++i) {
	Vector v = *i;
	if (v.Attribute()) {
#if (TALKATIVE>=2)
	  cerr << "Raising " << v << endl;
#endif
	  for (int j = 1; j<=n/2; j++) {
	    if (v.Attribute() & (1<<(j-1))) { // source erasing
	      int k = (n+1) - j;
	      RaisePPI(v, j, k, n, P, *Pold, *Pnew, false);
	    }
	  }
	}
      }

      Pn->insert(Pold->begin(), Pold->end());
      delete Pold;
      Pold = Pnew;
      Pnew = new SimpleVectorSet(expected_count);
    }
  } // P is dead now
  Pn->insert(Pold->begin(), Pold->end());
  delete Pold;
  delete Pnew;

  // Finally, add the base and destroy it
  Pn->insert(Pbase->begin(), Pbase->end());
  delete Pbase;
  return Pn;
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
  SimpleVectorSet *V;
  V = new SimpleVectorSet(1);
  Vector v(2); v(1) = -2, v(2) = +1; V->insert(v);
  for (int i = 2; i<n; i++) {
#if defined(WITH_STATS)
    hitcount[0] = failcount[0] = hitcount[1] = failcount[1] 
      = simplecount = goodcount = redcount = dupcount = 0;
#endif
    ppicount = 0;
    cerr << "### Extending to n = " << i+1 << endl;
    V = ExtendPPI(V, i);
    cerr << "### This makes " << ppicount 
	 << " PPI up to sign" 
#if defined(WITH_STATS)
	 << ", with " << redcount << " reduce ops, leading to "
	 << goodcount << " new vectors, "
	 << simplecount << " simple raises, "
	 << hitcount[0] << " red-erase hits, "
	 << failcount[0] << " red-erase fails, "
	 << hitcount[1] << " irr-erase hits, "
	 << failcount[1] << " irr-erase fails, and "
	 << dupcount << " duplicates. " 
#endif
	 << endl;
    if (i==n-1) { 
      cerr << "### Writing data file..." << flush;
      char fname[20];
      sprintf(fname, "ppi%d.dat", n);
      FILE *f = fopen(fname, "wb");
      char cn = n;
      fwrite(&cn, 1, 1, f);
      SimpleVectorSet::iterator j = V->begin();
      for (; j!=V->end(); ++j) {
	fwrite(&(*j)[0], 1, n, f);
      }
      fclose(f);
      cerr << "done." << endl;
    }
    cerr << "Elapsed time: " << user_time() << endl;
  }
}

/*
 * $Log$
 * Revision 1.32  2002/08/06 16:25:23  mkoeppe
 * Better memory mgmt
 *
 * Revision 1.31  2002/08/06 14:31:45  mkoeppe
 * Seems to work, but is slower
 *
 * Revision 1.30  2002/08/06 14:00:16  mkoeppe
 * OUR_OWN_HASH implementation. rudimentary.
 *
 * Revision 1.28.1.5.1.1.1.5.1.1  1999/10/14 14:31:30  mkoeppe
 * Added output to data file and time report.
 *
 * Revision 1.28.1.5.1.1.1.5  1999/03/23 16:05:51  mkoeppe
 * Clean-up.
 *
 * Revision 1.28.1.5.1.1.1.4  1999/03/23 12:03:41  mkoeppe
 * Digital trees no longer store the vectors but point to the Vectors
 * stored in the hash table. Reduces memory use to 64MB for n=20. Maybe
 * one should combine VectorAux and LeafNode to reduce memory even
 * more...
 *
 * Revision 1.28.1.5.1.1.1.3  1999/03/22 22:46:54  mkoeppe
 * Clean-up.
 *
 * Revision 1.28.1.5.1.1.1.2  1999/03/22 21:07:19  mkoeppe
 * New `compact' vectors. n=20 takes less than 70MB.
 *
 * Revision 1.28.1.5.1.1.1.1  1999/03/22 15:13:17  mkoeppe
 * Leaves of digital trees store the vectors directly, and only the coordinates
 * not known from the structure. n=19 now takes 47MB and 6m45s user
 * (everything measured on my Linux box with 47 Bogomips). n=20 seems to
 * take 74MB.
 *
 * Revision 1.28.1.5.1.1  1999/03/18 20:28:13  mkoeppe
 * Source erasing now `correct' (but not faster).
 *
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
