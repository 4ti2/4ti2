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

/* ------------------------------------------------------------------------ */
#include "myheader.h"
#include "print.h"
#include "vector.h"
/* ------------------------------------------------------------------------ */
int lengthListOrbit(listOrbit* LIST) {
  int len=0;

  while (LIST) {len++; LIST = LIST->rest;}
  return (len);
}
/* ----------------------------------------------------------------- */
orbit* createOrbit(vector v, listVector *L1, listVector *L2) {
  orbit *orb;

  orb=(orbit *)malloc(sizeof(orbit));
  orb->representative = v;
  orb->posOrbit = L1;
  orb->negOrbit = L2;
  orb->sizeOfOrbit = lengthListVector(L1)+lengthListVector(L2);
  orb->shortNorm = 0;
  orb->posNorm = 0;
  orb->negNorm = 0;
  orb->numOfPosEntries = 0;
  orb->numOfNegEntries = 0;
  orb->posEntriesOrdered = 0;
  orb->negEntriesOrdered = 0;
  return (orb);
}
/* ----------------------------------------------------------------- */
listOrbit* createListOrbit(vector v, listVector *L1, listVector *L2) {
  listOrbit *LIST;

  LIST=(listOrbit *)malloc(sizeof(listOrbit));
  LIST->first = createOrbit(v,L1,L2);
  LIST->rest = 0;
  return (LIST);
}
/* ----------------------------------------------------------------- */
orbit* softCopyOrbit(orbit *orb) {
  orbit *orb2;

  orb2=(orbit *)malloc(sizeof(orbit));
  orb2->representative = orb->representative;
  orb2->posOrbit = orb->posOrbit;
  orb2->negOrbit = orb->negOrbit;
  orb2->sizeOfOrbit = orb->sizeOfOrbit;
  orb2->shortNorm = orb->shortNorm;
  orb2->posNorm = orb->posNorm;
  orb2->negNorm = orb->negNorm;
  return (orb2);
}
/* ----------------------------------------------------------------- */
listOrbit** createArrayListOrbit(int numOfLists) {
  listOrbit** w;

  w = (listOrbit**)calloc(sizeof(listOrbit*),(numOfLists+1));
  if (w==0) exit(0);
  return (w);
}
/* ----------------------------------------------------------------- */
void printOrbit(orbit* orb, int numOfVars) {
  if (orb==0) printf("[]\n");
  printVector(orb->representative,numOfVars);
  printf("shortNorm = %d\n",orb->shortNorm);
  printf("posNorm = %d\n",orb->posNorm);
  printf("negNorm = %d\n",orb->negNorm);
/*     printf("Size of posOrbit = %d\n",lengthListVector(orb->posOrbit)); */
/*     printf("Size of negOrbit = %d\n",lengthListVector(orb->negOrbit)); */
/*   printf("\n"); */
/*   printListVector(orb->orbit,numOfVars); */
  printf("\n");
  return ;
}
/* ----------------------------------------------------------------- */
void printListOrbit(listOrbit* basis, int numOfVars) {
  if (basis==0) printf("[]\n");
  while(basis) {
    printOrbit(basis->first,numOfVars);
    basis = basis->rest;
  }
  return ;
}
/* ----------------------------------------------------------------- */
int lengthListCP(listCriticalPair* LIST) {
  int len=0;

  while (LIST) {len++; LIST = LIST->rest;}
  return (len);
}
/* ----------------------------------------------------------------- */
listCriticalPair* createListCriticalPair(orbit *v, orbit * w) {
  listCriticalPair *LIST;

  LIST=(listCriticalPair *)malloc(sizeof(listCriticalPair));
  LIST->first = v;
  LIST->second = w;
  LIST->maxShortNorm =0;
  LIST->rest = 0;
  return (LIST);
}
/* ------------------------------------------------------------------------ */
int isOrbitEqualToFullGroup(orbit *orb, listVector *S, int numOfVars) {
  vector v,p;
  listVector *tmp;

  v=orb->representative;
  tmp=S;
  while (tmp) {
    p=tmp->first;
    if ((isVectorEqualToPermutedVector(v,p,numOfVars)==1) 
	&& (isIdentityPermutation(p,numOfVars)==0))
	return (0);
    if ((isVectorEqualToNegativePermutedVector(v,p,numOfVars)==1) 
	&& (isIdentityPermutation(p,numOfVars)==0))
	return (0);
    tmp=tmp->rest;
  }

  return (1);
}
/* ------------------------------------------------------------------------ */
orbit* computeOrbitPermutationsGraver(orbit *orb, listVector *S, 
				      int numOfVars) {
  int i, lenS, lengthOfHeap;
  vector v, w, lastVector;
  vector *B, *P;
  listVector *tmp, *L, *endL;

/* Computes the permutations s from S such that orbit(v)=\cup {s}. */

/* Check, whether orbit(v)={s(v):s in S}. If so, let 
   orb->orbit=S. */

  if (isOrbitEqualToFullGroup(orb,S,numOfVars)==1) {
    orb->posOrbit=S;
    orb->sizeOfOrbit=lengthListVector(S);
    return (orb);
  }

/* If we are here, the orbit is NOT all of S. */
/* B contains the orbit vectors to check for multiple occurences.
   P contains the coresponding permutations. */

  lenS=lengthListVector(S);
  B=(vector*)calloc(sizeof(vector),lenS+1);
  P=(vector*)calloc(sizeof(vector),lenS+1);
  lengthOfHeap=0;

  v=orb->representative;
  tmp=S;
  while (tmp) {
    w=createVector(numOfVars);
    for (i=0; i<numOfVars; i++) w[i]=v[(tmp->first)[i]];
    w=lexPositiveVector(w,numOfVars);
    lengthOfHeap=updateHeapGraver(w,B,tmp->first,P,lengthOfHeap,numOfVars);
    tmp=tmp->rest;
  }

/* Now B contains the orbit, where vectors may occur more than once. */

  lastVector=copyVector(B[1],numOfVars);

/* L contains the permutations that form the orbit. */
  L=createListVector(P[1]);
  endL=L;
  lengthOfHeap=removeFirstHeapElement(B,P,lengthOfHeap,numOfVars);

  for (i=2; i<lenS+1; i++) {
    if (isVectorEqualToVector(lastVector,B[1],numOfVars)==0) {
      free(lastVector);
      lastVector=copyVector(B[1],numOfVars);
      endL->rest=createListVector(P[1]);
      endL=endL->rest;
    }
    lengthOfHeap=removeFirstHeapElement(B,P,lengthOfHeap,numOfVars);
  }
  free(lastVector);

  orb->posOrbit=L;
  orb->sizeOfOrbit=lengthListVector(L);

  free(B);
  free(P);

  return (orb);
}
/* ------------------------------------------------------------------------ */
orbit* computeOrbitPermutationsGroebner(orbit *orb, listVector *S, 
					int numOfVars) {
  int i, j, lenS, lengthOfHeap;
  vector v, w, lastVector;
  vector *B, *P;
  listVector *tmp, *tmp2, *L1, *endL1, *L2, *endL2;

/* Computes the permutations s from S such that orbit(v)=\cup {s}
   and puts them into two subsets such that the permutations in 
   the first set keep the ordering intact, whereas the permutations
   in the second set reverse the ordering. This helps in the
   reduction later. */

/* Check, whether orbit(v)={s(v):s in S}. If so, let 
   orb->orbit=S. */

  if ((orb->posOrbit)!=0) {
    tmp=orb->posOrbit;
    orb->posOrbit=0;
    while (tmp) {
      tmp2=tmp->rest;
      free(tmp->first);
      free(tmp);
      tmp=tmp2;
    }
  }

  if (isOrbitEqualToFullGroup(orb,S,numOfVars)==1) {
    L1=createListVector(0);
    endL1=L1;
    L2=createListVector(0);
    endL2=L2;

    v=orb->representative;
    tmp=S;
    w=createVector(numOfVars);
    while (tmp) {
      for(i=0; i<numOfVars; i++) w[i]=v[(tmp->first)[i]];
      if (isVectorCorrectlyOriented(w,0,numOfVars)==1) {
	endL1->rest=createListVector(tmp->first);
	endL1=endL1->rest;
      } else {
	endL2->rest=createListVector(tmp->first);
	endL2=endL2->rest;
      }
      tmp=tmp->rest;
    }
    free(w);

    orb->posOrbit=L1->rest;
    free(L1);
    orb->negOrbit=L2->rest;
    free(L2); 
    orb->sizeOfOrbit=lengthListVector(S);
    return (orb);
  }

/* If we are here, the orbit is NOT all of S.
   L1, L2 contain the permutations that form the orbit.
   B contains the orbit vectors to check for multiple occurences.
   P contains the coresponding permutations. */

  lenS=lengthListVector(S);
  B=(vector*)calloc(sizeof(vector),lenS);
  P=(vector*)calloc(sizeof(vector),lenS);
  lengthOfHeap=0;

  v=orb->representative;
  tmp=S;
  while (tmp) {
    w=createVector(numOfVars);
    for(i=0; i<numOfVars; i++) w[i]=v[(tmp->first)[i]];
    w=lexPositiveVector(w,numOfVars);
    lengthOfHeap=updateHeapGraver(w,B,tmp->first,P,lengthOfHeap,numOfVars);
    tmp=tmp->rest;
  }

/* Now B contains the orbit, where vectors may occur more than once. */

  L1=createListVector(0);
  endL1=L1;
  L2=createListVector(0);
  endL2=L2;

  lastVector=copyVector(B[1],numOfVars);

/* L1, L2 contain the permutations that form the orbit. */

  w=createVector(numOfVars);
  for(i=0; i<numOfVars; i++) w[i]=v[(P[1])[i]];
  if (isVectorCorrectlyOriented(w,0,numOfVars)==1) {
    endL1->rest=createListVector(P[1]);
    endL1=endL1->rest;
  } else {
    endL2->rest=createListVector(P[1]);
    endL2=endL2->rest;
  }
  free(w);

  lengthOfHeap=removeFirstHeapElement(B,P,lengthOfHeap,numOfVars);

  for (i=2; i<lenS+1; i++) {
    if (isVectorEqualToVector(lastVector,B[1],numOfVars)==0) {
      free(lastVector);
      lastVector=copyVector(B[1],numOfVars);
      w=createVector(numOfVars);
      for(j=0; j<numOfVars; j++) w[j]=v[(P[1])[j]];
      if (isVectorCorrectlyOriented(w,0,numOfVars)==1) {
	endL1->rest=createListVector(P[1]);
	endL1=endL1->rest;
      } else {
	endL2->rest=createListVector(P[1]);
	endL2=endL2->rest;
      }
      free(w);
    }
    lengthOfHeap=removeFirstHeapElement(B,P,lengthOfHeap,numOfVars);
  }
  free(lastVector);

  orb->posOrbit=L1->rest;
  free(L1);
  orb->negOrbit=L2->rest;
  free(L2);
  orb->sizeOfOrbit=lengthListVector(orb->posOrbit)+
      lengthListVector(orb->negOrbit);

  free(B);
  free(P);

  return (orb);
}
/* ------------------------------------------------------------------------ */
listVector* computeOrbit(vector v, listVector* permutations, int numOfVars) {
  int i;
  vector w;
  listVector *tmp, *L, *endL;

/* Computes the full orbit from representative v and permutations. */
/*   printf("Computing orbit of: "); */
/*   printVector(v,numOfVars); */

  L=createListVector(0);
  endL=L;

  tmp=permutations;

  while (tmp) {
    w=copyVector(v,numOfVars);
    for(i=0; i<numOfVars; i++) w[i]=v[(tmp->first)[i]];
    endL->rest=createListVector(w);
    endL=endL->rest;
    tmp=tmp->rest;
  }

  tmp=L;
  L=L->rest;
  free(tmp);

  return (L);
}
/* ------------------------------------------------------------------------ */
listOrbit* eliminateMultipleRepresentatives(listOrbit *basis, int numOfVars) {
  int ok;
  vector g;
  listVector *tmpOrb, *L, *L2;
  orbit *orb;
  listOrbit *tmp, *tmp2, *endBasis;

  tmp=basis;
  basis=createListOrbit(0,0,0);
  endBasis=basis;
  while (tmp) {
    orb=tmp->first;
    tmpOrb=computeOrbit(orb->representative,orb->posOrbit,numOfVars);
    ok=1;
    tmp2=basis->rest;
    while ((tmp2) && (ok==1)) {
      g=(tmp2->first)->representative;
      if ((isVectorInListVector(g,tmpOrb,numOfVars)!=0) || 
	  (isNegativeVectorInListVector(g,tmpOrb,numOfVars)!=0)) {
	ok=0;
      } 
      tmp2=tmp2->rest;
    }
    L=tmpOrb;
    while (tmpOrb) {
      tmpOrb=tmpOrb->rest;
      free(L->first);
      free(L);
      L=tmpOrb;
    }

    if (ok==1) {
      if (orb->negOrbit) {
	tmpOrb=computeOrbit(orb->representative,orb->negOrbit,numOfVars);
	tmp2=basis->rest;
	while ((tmp2) && (ok==1)) {
	  g=(tmp2->first)->representative;
	  if ((isVectorInListVector(g,tmpOrb,numOfVars)!=0) || 
	      (isNegativeVectorInListVector(g,tmpOrb,numOfVars)!=0)) {
	    ok=0;
	  } 
	  tmp2=tmp2->rest;
	}
	L=tmpOrb;
	while (tmpOrb) {
	  tmpOrb=tmpOrb->rest;
	  free(L->first);
	  free(L);
	  L=tmpOrb;
	}
      }
    }
    
    if (ok==1) {
      endBasis->rest=tmp;
      tmp=tmp->rest;
      endBasis=endBasis->rest;
      endBasis->rest=0;
    } else {
      tmp2=tmp;
      tmp=tmp->rest;
      orb=tmp2->first;
      tmp2->rest=0;
      free(tmp2);
      free(orb->representative);
      L=orb->posOrbit;
      while (L) {
	L2 = L->rest;
	free(L);
	L=L2;
      }
      L=orb->negOrbit;
      while (L) {
	L2 = L->rest;
	free(L);
	L=L2;
      }
      free(orb);
    }
    
  }

  tmp=basis;
  basis=basis->rest;
  free(tmp);

  return (basis);
}
/* ------------------------------------------------------------------------ */
vector canonicalRepresentative(vector v, listVector* permutations, 
			       int numOfVars) {
  int i;
  vector w, canonicalRepresentative;
  listVector *tmp;

  tmp=permutations;
  canonicalRepresentative=copyVector(v,numOfVars);

  while (tmp) {
    w=copyVector(v,numOfVars);
    for(i=0; i<numOfVars; i++) w[i]=v[(tmp->first)[i]];
    if (compareVectorsByLex(canonicalRepresentative,w,numOfVars)==-1) {
      free(canonicalRepresentative);
      canonicalRepresentative=w;
    } else free(w);
    w=copyVector(v,numOfVars);
    for(i=0; i<numOfVars; i++) w[i]=-v[(tmp->first)[i]];
    if (compareVectorsByLex(canonicalRepresentative,w,numOfVars)==-1) {
      free(canonicalRepresentative);
      canonicalRepresentative=w;
    } else free(w);
    tmp=tmp->rest;
  }

  return (canonicalRepresentative);
}
/* ------------------------------------------------------------------------ */
vector canonicalRepresentativeAndShortNorm(vector v, listVector* permutations,
					   vector normIndices, int *shortNorm, 
					   int numOfVars, int dimension, 
					   int minNorm) {
  int i,norm,s;
  vector w, canonicalRepresentative;
  listVector *tmp;

  tmp=permutations;
  canonicalRepresentative=copyVector(v,numOfVars);

  norm=normOfVector(v,numOfVars);

  while (tmp) {
    w=copyVector(v,numOfVars);
    for(i=0; i<numOfVars; i++) w[i]=v[(tmp->first)[i]];

/* Check short norm. */
    s=0;
    for (i=0; i<dimension; i++) s=s+abs(w[normIndices[i]]);
    if (s<norm) {
      norm=s;
      if (norm<minNorm) {
	free(w);
	free(canonicalRepresentative);
	(*shortNorm)=0;
	return (0);
      }
    }
/* Check lex ordering to get canonical representative. */
    if (compareVectorsByLex(canonicalRepresentative,w,numOfVars)==-1) {
      free(canonicalRepresentative);
      canonicalRepresentative=w;
    } else free(w);
    tmp=tmp->rest;
  }

  (*shortNorm)=norm;

  return (canonicalRepresentative);
}
/* ----------------------------------------------------------------- */
listVector* extractSymmetryRepresentatives(listVector *basis,
					   listVector *permutations, 
					   int numOfVars) {
  int i,numOfVectorsConsidered,numOfRepresentatives,lenBasis,maxNorm,norm;
  vector v;
  listVector *tmp, *tmp2, *orbitV, *representatives, *endRepresentatives;
  listVector **arrayBasis;

  lenBasis=lengthListVector(basis);
  printf("basis = %d elements, symmGroup = %d permutations\n",lenBasis,
	 lengthListVector(permutations));

  maxNorm=maximalNormInListVector(basis,numOfVars);
  printf("Maximum appearing norm: %d\n\n",maxNorm);

  arrayBasis=createArrayListVector(maxNorm+1);
  for (i=0;i<maxNorm+1;i++) arrayBasis[i]=0;

  tmp=basis;
  while (tmp) {
    v=copyVector(tmp->first,numOfVars);
    norm=normOfVector(v,numOfVars);
    tmp2=createListVector(v);
    tmp2->rest=arrayBasis[norm];
    arrayBasis[norm]=tmp2;
    tmp=tmp->rest;
  }

  representatives=createListVector(0);
  endRepresentatives=representatives;

  numOfVectorsConsidered=0;
  numOfRepresentatives=0;


  for (i=0;i<maxNorm+1;i++) {
    printf("Considering norm: %d,   Number of vectors: %d\n",i,
	   lengthListVector(arrayBasis[i]));
    tmp=arrayBasis[i];
    while (tmp) {
      if (numOfVectorsConsidered==100*(numOfVectorsConsidered/100)) {
	printf("%d / %d considered.   %d representatives found so far.\n",
	       numOfVectorsConsidered,lenBasis,numOfRepresentatives);
      }
      v=tmp->first;
      if (v!=0) {
	numOfRepresentatives++;
	endRepresentatives->rest=createListVector(v);
	endRepresentatives=endRepresentatives->rest;
	orbitV=computeOrbit(v,permutations,numOfVars);
	tmp2=tmp->rest;
	while (tmp2) {
	  if (tmp2->first!=0)
	    if (isVectorInListVector(tmp2->first,orbitV,numOfVars)==1) {
	      free(tmp2->first);
	      tmp2->first=0;
	    }
	  tmp2=tmp2->rest;
	}
      }
      numOfVectorsConsidered++;
      tmp=tmp->rest;
    }
  }
  return (representatives->rest);
}
/* ----------------------------------------------------------------- */
listVector* expandRepresentativeIntoFullOrbits(listVector *representatives,
					       listVector *permutations, 
					       int numOfVars, int infoLevel) {
  int i,numOfVectorsConsidered,lenRepresentatives;
  vector v, p, rep;
  listVector *tmp, *basis, *endBasis, *newBasis, *endNewBasis;

  lenRepresentatives=lengthListVector(representatives);
if (infoLevel>0) {
  printf("Number of representatives = %d\n",lenRepresentatives);
}

  v=0;
  basis=createListVector(0);
  endBasis=basis;

  numOfVectorsConsidered=0;
  while (representatives) {
    rep=representatives->first;
    if (isVectorInListVector(rep,basis->rest,numOfVars)==0) {
      tmp=permutations;
      newBasis=createListVector(0);
      endNewBasis=newBasis;
      while (tmp) {
	p=tmp->first;
	v=createVector(numOfVars);
	for(i=0;i<numOfVars;i++) v[i]=rep[p[i]];
	v=lexPositiveVector(v,numOfVars);
	if (isVectorInListVector(v,newBasis->rest,numOfVars)==0) {
	  endNewBasis->rest=createListVector(v);
	  endNewBasis=endNewBasis->rest;
	} else {
	  free(v); 
	}
	tmp=tmp->rest;
      }
if (infoLevel>0) {
      printf("%d new basis vectors found.\n",lengthListVector(newBasis->rest));
}
      endBasis->rest=newBasis->rest;
      while (endBasis->rest) endBasis=endBasis->rest;
    }
    representatives=representatives->rest;
    numOfVectorsConsidered++;
if (infoLevel>0) {
    if (numOfVectorsConsidered==1*(numOfVectorsConsidered/1)) {
      printf("%d / %d considered.   %d basis vectors found so far.\n",
	     numOfVectorsConsidered,lenRepresentatives,
	     lengthListVector(basis->rest));
    }
}
  }

if (infoLevel>0) {
  printf("Done.   %d basis vectors found.\n",
	 lengthListVector(basis->rest));
}

  return (basis->rest);
}
/* ----------------------------------------------------------------- */
