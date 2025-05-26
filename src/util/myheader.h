/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Raymond Hemmecke.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 
*/

#include <errno.h>
#include <stdio.h>
#include <math.h>
//#include <malloc.h>
#include <string.h>
#include <stdlib.h>

#define startPrintingTime -0.1
#define numOfBits 32
#define standardInfoLevel 10
#define maxNumOfVars 0
#define largePrime 32003
#define sizeX 3
#define sizeY 3
#define sizeZ 5

typedef int *vector;

typedef struct listVector {
  int *first;
  int *posSupport;
  int *negSupport;
  int posNorm;
  int negNorm;
  int sign;
  struct listVector *rest;
}
listVector;

typedef struct sparseVector {
  int posSize;
  int *posIndex;
  int *posEntry;
  int negSize;
  int *negIndex;
  int *negEntry;
}
sparseVector;

typedef struct listSparseVector {
  sparseVector *first;
  int sign;
  struct listSparseVector *rest;
}
listSparseVector;

typedef struct rationalVector {
  int* enumerator;
  int* denominator;
} rationalVector;

typedef struct listRationalVector {
  rationalVector* first;
  struct listRationalVector *rest;
} listRationalVector;

typedef struct tree {
  int index;
  int numOfSubTrees;
  vector entries;
  struct tree **subTrees;
}
tree;

typedef struct listPair {
  listVector* first;
  listVector* second;
  struct listPair* rest;
}
listPair;

typedef struct orbit {
  vector representative;
  listVector* posOrbit;
  listVector* negOrbit;
  int sizeOfOrbit;
  int posNorm;
  int negNorm;
  int shortNorm;
  int numOfPosEntries;
  int numOfNegEntries;
  vector posEntriesOrdered;
  vector negEntriesOrdered;
}
orbit;

typedef struct listOrbit {
  orbit* first;
  struct listOrbit* rest;
}
listOrbit;

typedef struct listCriticalPair {
  orbit *first;
  orbit *second;
  int maxShortNorm;
  struct listCriticalPair *rest;
}
listCriticalPair;

typedef struct heapElement {
  vector first;
  int x;
  int y;
  int summand;
}
heapElement;

typedef listVector **hashTable;

typedef struct problem {
  listVector* equations;
  listVector* inequalities;
  listVector* congruences;
}
problem;


typedef struct listProblem {
  problem* first;
  struct listProblem* rest;
}
listProblem;

typedef struct signPattern {
  char* types;
}
signPattern;

typedef struct listSignPattern {
  signPattern* first;
  struct listSignPattern* rest;
}
listSignPattern;

typedef struct listSignPatternPair {
  signPattern* first;
  listSignPattern* second;
  struct listSignPatternPair* rest;
}
listSignPatternPair;


#define myxasprintf(fILEnANEpTR,fORMAT,...) { \
    errno=0; \
    if (asprintf(fILEnANEpTR,fORMAT,__VA_ARGS__) < 0) { \
      printf("Error allocating file name (%s).\n",strerror(errno)); \
      exit(0); \
    } \
  }
