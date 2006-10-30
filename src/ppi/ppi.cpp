// Computation of the primitive partition identities (PPI)
// Copyright (C) 1998, 1999, 2002, 2003 Matthias Koeppe <mkoeppe@mail.math.uni-magdeburg.de>
// Copyright (C) 2006 4ti2 team.
// $Id$

// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; see the file COPYING. If not, write to the Free
// Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA.

// For an introduction to PPIs and for the algorithm used in this
// program, see Haus, Koeppe, Weismantel: A Primal All-Integer
// Algorithm Based on Irreducible Solutions, Math. Programming, Series
// B, 96 (2003), no. 2, pp. 205-246
//
// See also http://www.math.uni-magdeburg.de/~mkoeppe/art/ppi/

#define _NOTHREADS
// STL would use slow mutexes 

#include <cstdio>
#include <cassert>
#include <set>
#include <vector>
#include <iostream>
#include <iomanip>
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

typedef unsigned int attr_type; /* 32 bits are enough */

#if defined(COMPACT_VECTORS)
struct VectorAux {
  /* order important */
  attr_type Attribute;
  signed char Stuff[4];
  void *operator new(size_t s, int length) {
    if (length < 3) length=3;
    return malloc(s - 3 +length);
  }
  void operator delete(void *v) {
    free(v);
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
      if (aux == NULL) 
	aux = new(v.aux->Length()) VectorAux(*v.aux);
      else if (aux->Length() == v.aux->Length())
	memcpy(aux, v.aux, sizeof(VectorAux) -4 + 1 + v.aux->Length());
      else {
	delete aux;
	aux = new(v.aux->Length()) VectorAux(*v.aux);
      }
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
  attr_type &Attribute() { return aux->Attribute; }
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
    index(ndex), the_set(theset), node(ode) {}
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
  SimpleVectorSet ();
  SimpleVectorSet (const SimpleVectorSet &);
  SimpleVectorSet &operator=(const SimpleVectorSet &);
public:
  typedef SimpleVectorSetIterator iterator;
  Vector **buckets;
  size_t num_buckets;
  size_t count;
  SimpleVectorSet(size_t size) {
    count = 0;
    if (size == 0) size =1;
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
  pair<iterator, bool> insert(Vector *obj)  /* owns obj afterwards */
  {
    size_t h = vector_hash_fun(*obj);
    size_t hm = h % num_buckets;
    /* See if it's there */
    Vector *v;
    for (v = buckets[hm]; v!=NULL; v=v->next) {
      if (v->hash == h) {
	if (*v == *obj) {
	  return pair<iterator, bool>(SimpleVectorSetIterator(this, hm, v), false);
	}
      }
    }
    /* Insert it at the front of the list */
    v = obj;
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

class VerySimpleVectorSetIterator {
  Vector *p;
public:
  VerySimpleVectorSetIterator(Vector *pp) : p(pp) {}
  Vector &operator*() { return *p; }
  friend bool operator == (const VerySimpleVectorSetIterator &a, const VerySimpleVectorSetIterator &b)
  { return a.p == b.p; }
  friend bool operator != (const VerySimpleVectorSetIterator &a, const VerySimpleVectorSetIterator &b)
  { return a.p != b.p; }
  VerySimpleVectorSetIterator &operator ++() { p = p->next; return *this; }
};

class VerySimpleVectorSet { /* a linked list */
  Vector *head;
public:
  size_t count;
  typedef VerySimpleVectorSetIterator iterator;
  VerySimpleVectorSet() : head(NULL), count(0) {}
  ~VerySimpleVectorSet() {
    Vector *v, *n;
    for (v = head; v!=NULL; ) {
      n = v->next;
      delete v;
      v = n;
    }    
  }
  iterator begin() const { return VerySimpleVectorSetIterator(head); }
  iterator end() const { return NULL; }
  void insert_destructively(SimpleVectorSet &s) {
    /* Move all entries from s.  After that, s is empty. */
    size_t i;
    for (i = 0; i<s.num_buckets; i++) {
      Vector *p = s.buckets[i];
      if (p != NULL) {
	while (p->next != NULL) p = p->next;
	p->next = head;
	head = s.buckets[i];
	s.buckets[i] = NULL;
      }
    }
    count += s.count;
    s.count = 0;
  }
  void insert(const Vector &s)
  {
    Vector *v = new Vector(s);
    v->next = head;
    head = v;
    count++;
  }
};

#else
#if defined(HASH) // Use a hash table for SimpleVectorSet

#include <hash_set>

struct hash<Vector> 
{
  size_t operator()(const Vector &v) const {
    // FIXME: Is this a good enough hash function? 
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
    advance(adv), isleaf(leaf) {}
  void *operator new(size_t s) {
    return malloc(s);
  }
  void operator delete(void *v) {
    free(v);
  }
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
    root(new InnerNode(0)),
    AdvanceNodes(0),
    Dimension(dimension) {}
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
  assert(0);
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
  attr_type attrib = 0;
  attr_type mask = 1;
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

Vector EraseSource_u;
Vector EraseSource_mu;

inline void EraseSource(const Vector &w, int n, SimpleVectorSet &Pold,
			bool irreducible, bool WithNegative)
{
#if !defined(ERASE_SOURCES_OF_IRREDUCIBLES)
  if (!irreducible) return; // they often fail, which is expensive
#endif
  //  cerr << "EraseSource " << w << endl; 
  EraseSource_u = w;
  Vector &u = EraseSource_u;
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
      Vector &mu = EraseSource_mu;
      if (2*LastNonzeroPos == n+1) { // Very special case:
	// Need not erase sources because this case is left out
	// anyway. (The code below doesn't grok this case.)
	//mu = -u;
	for (int j = 1; j<=n+1; j++)
	  mu(j) = - u(j);
      }
      else {
	// Handle LastNonzeroPos:
	u(i)++, u(LastNonzeroPos)++; // This might zero it, so:
	//mu = -u;
	for (int j = 1; j<=n+1; j++)
	  mu(j) = - u(j);
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

Vector *temp = NULL;

inline void RaisePPI(const Vector &v, int j, int k, int n, 
		     VectorSet &P, SimpleVectorSet &Pold,
		     SimpleVectorSet &Pnew, bool WithNegative,
		     bool KnownIrreducible = false)
{
  // assert(j<k);
  //  cerr << "RaisePPI " << v << " with " << j << " and " << k << endl;
  /* Vector *ww = new Vector(v); */
  if (temp == NULL)
    temp = new Vector(v);
  else 
    *temp = v;
  Vector &w = *temp;
  
  w(n+1)++, w(j)--, w(k)--;
  // reducibility check
  if (KnownIrreducible || (w(j) >= 0 && w(k) >= 0)) { // w is irreducible 
#if defined(WITH_STATS)
    simplecount++;
#endif
    SetupAttribute(w, n, false);
    if (Pnew.insert(temp).second) {
      EraseSource(w, n, Pold, true, WithNegative);
      reportx(w);
      temp = NULL;
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
	if (Pnew.insert(temp).second) {
#if defined(WITH_STATS)
	  goodcount++;
#endif
	  EraseSource(w, n, Pold, true, WithNegative);
	  reportx(w);
	  temp = NULL;
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

VerySimpleVectorSet *ExtendPPI(VerySimpleVectorSet *Pn, int n)
{
  SimpleVectorSet *Pold;
  int expected_count = Pn->count/20; // 2 * Pn->count;
  if (expected_count < 20000) expected_count = 20000;
  SimpleVectorSet *Pnew = new SimpleVectorSet(expected_count);
  SimpleVectorSet *Pbase = new SimpleVectorSet(Pn->count);

  // Implementation of Algorithms 4.3.8, 4.3.11 from [Urbaniak] 
  
  {
    VectorSet P(n);
    /* initialize special vectors */
    delete temp; temp = NULL;
    rangemin = Vector(n+1), rangemax = Vector(n+1);
    EraseSource_u = Vector(n+1);
    EraseSource_mu = Vector(n+1);
    Vector minus_v(n+1);

    // (1) Create the range-searchable set P of (n+1)-vectors from the
    // set Pn of n-vectors.
    cerr << "# Vectors copied from n = " << n << ": " << endl;
    for (VerySimpleVectorSet::iterator i = Pn->begin(); i!=Pn->end(); ++i) {
      Vector *v = new Vector(n+1);
      for (int j = 1; j <= n; j++) (*v)(j) = (*i)(j);
      (*v)(n+1) = 0;
      SetupAttribute(*v, n, true); 
      P.insert(*Pbase->insert(v).first); reportx(*v);
    }
    // destroy Pn to save memory
    delete Pn;
    Pn = new VerySimpleVectorSet();
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
    cerr << "# Vectors of P" << 1 << "(" << n+1 << "):" << flush;

    // first a pass for the simple cases
    for (SimpleVectorSet::iterator i = Pbase->begin(); i!=Pbase->end(); ++i) {
      Vector &v = *i;
      if (v.Attribute()) {
#if (TALKATIVE>=2)
	cerr << "Raising " << v << endl;
#endif
	//Vector minus_v = -v;
	for (int j=1; j<=n+1; j++)
	  minus_v(j) = - v(j);
	for (int j = 1; j<=n/2; j++) { // yes n/2 is enough
	  int k = (n+1) - j;
	  if (v(j) > 0 && v(k) > 0)
	    RaisePPI(v, j, k, n, P, *Pbase, *Pnew, true, true);
	  if (v(j) < 0 && v(k) < 0)
	    RaisePPI(minus_v, j, k, n, P, *Pbase, *Pnew, true, true);
	}
      }
    }

    // then a pass for the complicated cases
    for (SimpleVectorSet::iterator i = Pbase->begin(); i!=Pbase->end(); ++i) {
      Vector &v = *i;
      if (v.Attribute()) {
#if (TALKATIVE>=2)
	cerr << "Raising " << v << endl;
#endif
	for (int j=1; j<=n+1; j++)
	  minus_v(j) = - v(j);
	for (int j = 1; j<=n/2; j++) { // yes n/2 is enough
	  int k = (n+1) - j;
	  if (v.Attribute() & (1<<(j-1)))
	    RaisePPI(v, j, k, n, P, *Pbase, *Pnew, true);
	  if (v.Attribute() & ((1<<16)<<(j-1)))
	    RaisePPI(minus_v, j, k, n, P, *Pbase, *Pnew, true);
	}
      }
    }
    cerr << Pnew->count << endl;
    //
    Pold = Pnew;
    Pnew = new SimpleVectorSet(expected_count);

    // (4) Build all other primitive identities with exactly (t+1)
    // components of (n+1).

    for (int t = 1; t<n; t++) {
      cerr << "# Vectors of P" << t+1 << "(" << n+1 << "):" << flush;
      // simple cases
      for (SimpleVectorSet::iterator i = Pold->begin(); i!=Pold->end(); ++i) {
	Vector &v = *i;
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
	Vector &v = *i;
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
      cerr << Pnew->count << endl;
      Pn->insert_destructively(*Pold);
      SimpleVectorSet *tmp = Pold; /* an empty set */
      Pold = Pnew;
#ifdef RECYCLE_SETS
      Pnew = tmp;
#else
      delete tmp;
      Pnew = new SimpleVectorSet(1 + Pold->count / 4);
#endif
    }
  } // P is dead now
  Pn->insert_destructively(*Pold);
  delete Pold;
  delete Pnew;

  // Finally, add the base and destroy it
  Pn->insert_destructively(*Pbase);
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
  VerySimpleVectorSet *V;
  V = new VerySimpleVectorSet();
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
      VerySimpleVectorSet::iterator j = V->begin();
      for (; j!=V->end(); ++j) {
	fwrite(&(*j)[0], 1, n, f);
      }
      fclose(f);
      cerr << "done." << endl;
    }
    cerr << "Elapsed time: " << user_time() << endl;
  }
  delete V;
}

