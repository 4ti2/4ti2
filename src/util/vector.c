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

/* ----------------------------------------------------------------- */
/*                                                                   */
/* Functions on numbers, vectors and listVectors.                    */
/*                                                                   */
/* Author   : Raymond Hemmecke   (implementation)                    */
/* Co-Author: Ralf Hemmecke      (data structure, code optimization) */
/*                                                                   */
/*------------------------------------------------------------------ */
#include "myheader.h"
#include "print.h"
/* ----------------------------------------------------------------- */
int gcd(int a, int b) {
  int c;

  if ((a<0) || (b<0)) {
    printf("Error in GCD, a = %d, b = %d\n",a,b);
    exit(0);
  }
  if (b==0) return (a);

  if (a<b) return (gcd(b,a));

  c=a/b;
  return (gcd(b,a-c*b));
}
/* ----------------------------------------------------------------- */
int lcm(int a, int b) {
  int c;

  if ((a==0) || (b==0)) return (0);

  /*  printf("a = %d, b = %d\n",a,b); */

  c=gcd(a,b);
  a=a/c;

  if (a<0) {
    printf("LCM, a = %d, b = %d\n",a,b);
    exit(0);
  }

  if (a*b<0) {
    printf("LCM, a*b, a = %d, b = %d\n",a,b);
    exit(0);
  }

  return (a*b);
}
/* ----------------------------------------------------------------- */
vector createVector(int numOfVars) {
  vector w;

  w = (vector)malloc(sizeof(int)*(numOfVars+1));
  if (w==0) {
    printf("Could not allocate memory in function createVector.\n");
    printf("Please report this problem.\n");
    exit(0);
  }
  return (w);
}
/* ----------------------------------------------------------------- */
int dotProduct(vector v, vector w, int numOfVars) {
  int i,s;

  s=0;
  for (i=0; i<numOfVars; i++) s=s+v[i]*w[i];
  
  return(s);
}
/* ----------------------------------------------------------------- */
vector posVector(vector v, int numOfVars) {
  int i;
  vector w;

  w=createVector(numOfVars);
  for (i=0; i<numOfVars; i++) {
    if (v[i]>0) w[i]=v[i]; else w[i]=0;
  }

  return (w);
}
/* ----------------------------------------------------------------- */
vector negVector(vector v, int numOfVars) {
  int i;
  vector w;

  w=createVector(numOfVars);
  for (i=0; i<numOfVars; i++) {
    if (v[i]<0) w[i]=-v[i]; else w[i]=0;
  }

  return (w);
}
/* ----------------------------------------------------------------- */
vector removeGCDfromVector(vector w, int numOfVars) {
  int i,g;

  g=0;

  for (i=0; i<numOfVars; i++) g=gcd(g,abs(w[i]));

/*   if (g>1) { */
/*     printf("g = %d\n",g); */
/*     printVector(w,numOfVars); */
/*   } */

  for (i=0; i<numOfVars; i++) w[i]=w[i]/g;

/*   if (g>1) { */
/*     printVector(w,numOfVars); */
/*   } */

  return(w); 
}
/* ----------------------------------------------------------------- */
vector supportOfVector(vector v, int numOfVars, int numOfBytes) {
  int i,j,pos,s;
  vector supp;

  supp=createVector(numOfBytes);

  if (numOfVars==(numOfBytes*numOfBits)) {
    for (j=0; j<numOfBytes; j++) {
      pos=0;
      for (i=0; i<numOfBits; i++) {
        pos = pos << 1;
	if (v[j*numOfBits+i]!=0) pos++;
      }
      supp[j]=pos;
    }
  } else {
    for (j=0;j<numOfBytes-1;j++) {
      pos=0;
      for (i=0; i<numOfBits; i++) {
        pos = pos << 1;
	if (v[j*numOfBits+i]!=0) pos++;
     }
      supp[j]=pos;
    }
    pos=0;
    s=(numOfBytes-1)*numOfBits;
    for (i=0; i<numOfVars-s; i++) {
      pos = pos << 1;
      if (v[s+i]!=0) pos++;
    }
    supp[numOfBytes-1]=pos;
  }

  return (supp);
}
/* ----------------------------------------------------------------- */
vector negateSupportVector(vector v, int numOfBytes) {
  int i;

  for (i=0; i<numOfBytes; i++) {
      v[i]=(~(v[i]));
  }

  return (v);
}
/* ----------------------------------------------------------------- */
vector positiveSupportOfVector(vector v, int numOfVars, int numOfBytes) {
  int i,j,pos,s;
  vector supp;

  supp=createVector(numOfBytes);

  if (numOfVars==(numOfBytes*numOfBits)) {
    for (j=0; j<numOfBytes; j++) {
      pos=0;
      for (i=0; i<numOfBits; i++) {
        pos = pos << 1;
	if (v[j*numOfBits+i]>0) pos++;
      }
      supp[j]=pos;
    }
  } else {
    for (j=0;j<numOfBytes-1;j++) {
      pos=0;
      for (i=0; i<numOfBits; i++) {
        pos = pos << 1;
	if (v[j*numOfBits+i]>0) pos++;
     }
      supp[j]=pos;
    }
    pos=0;
    s=(numOfBytes-1)*numOfBits;
    for (i=0; i<numOfVars-s; i++) {
      pos = pos << 1;
      if (v[s+i]>0) pos++;
    }
    supp[numOfBytes-1]=pos;
  }

  return (supp);
}
/* ----------------------------------------------------------------- */
vector negativeSupportOfVector(vector v, int numOfVars, int numOfBytes) {
  int i,j,neg,s;
  vector supp;

  supp=createVector(numOfBytes);

  if (numOfVars==(numOfBytes*numOfBits)) {
    for (j=0; j<numOfBytes; j++) {
      neg=0;
      for (i=0; i<numOfBits; i++) {
        neg = neg << 1;
	if (v[j*numOfBits+i]<0) neg++;
      }
      supp[j]=neg;
    }
  } else {
    for (j=0; j<numOfBytes-1; j++) {
      neg=0;
      for (i=0; i<numOfBits; i++) {
        neg = neg << 1;
	if (v[j*numOfBits+i]<0) neg++;
      }
      supp[j]=neg;
    }
    neg=0;
    s=(numOfBytes-1)*numOfBits;
    for (i=0; i<numOfVars-s; i++) {
      neg = neg << 1;
      if (v[s+i]<0) neg++;
    }
    supp[numOfBytes-1]=neg;
  }

  return (supp);
}
/* ----------------------------------------------------------------- */
int normOfBinaryVector(vector v, int numOfVars) {
  int i,j,k,m,norm=0;

  for (i=0; i<numOfVars; i++) {
    k=v[i];
    for (j=0; j<numOfBits; j++) {
      m=k-2*(k/2);
      k=(k-m)/2;
      norm=norm+m;
    }
  }

  return (norm);
}
/* ----------------------------------------------------------------- */
int positiveNormOfVector(vector v, int numOfVars) {
  int i, norm=0;

  for (i=0; i<numOfVars; i++) if (v[i]>0) norm+=v[i];
  return (norm);
}
/* ----------------------------------------------------------------- */
int negativeNormOfVector(vector v, int numOfVars) {
  int i, norm=0;

  for (i=0; i<numOfVars; i++) if (v[i]<0) norm-=v[i];
  return (norm);
}
/* ----------------------------------------------------------------- */
int hasCommonFactor(vector x, vector y, int numOfVars) {
  int i;

  for (i=0; i<numOfVars; i++)
    if ((x[i]>0) && (y[i]>0)) return (1);

  return (0);
}
/* ----------------------------------------------------------------- */
listVector* appendVectorToListVector(vector v, listVector *REST) {
  listVector *LIST;

  LIST=(listVector *)malloc(sizeof(listVector));
  if (LIST==0) {
    printf("Could not allocate memory in appendVectorToListVector.\n");
    printf("Please report this problem.\n");
    exit(0);
  }
  LIST->first = v;
  LIST->sign = 0;
  LIST->rest = REST;
  return (LIST);
}
/* ----------------------------------------------------------------- */
void freeVector(vector v) {
  free(v);
  return;
}
/* ----------------------------------------------------------------- */
void freeListVector(listVector *v) {
  if ((v!=0) && (v->first!=0)) {
    freeVector(v->first);
    free(v);
  }
  return;
}
/* ----------------------------------------------------------------- */
void freeAllOfListVector(listVector *v) {
  listVector *tmp;

  while (v) {
    if ((v->first)!=0) {
      free(v->first);
    }
    tmp=v;
    v=v->rest;
    free(tmp);
  }

  return;
}
/* ----------------------------------------------------------------- */
vector* createArrayVector(int numOfVectors) {
  vector* w;

  w = (vector*)malloc(sizeof(vector)*(numOfVectors+1));
  if (w==0) {
    printf("Could not allocate memory in createArrayVector.\n");
    printf("Please report this problem.\n");
    exit(0);
  }
  return (w);
}
/* ----------------------------------------------------------------- */
listVector** createArrayListVector(int numOfLists) {
  listVector** w;

  w = (listVector**)calloc(sizeof(listVector*),(numOfLists+1));
  if (w==0) {
    printf("Could not allocate memory in createArrayListVector.\n");
    printf("Please report this problem.\n");
    exit(0);
  }
  return (w);
}
/* ----------------------------------------------------------------- */
int isZeroVector(vector v, int numOfVars) {
  int i;

  for (i=0; i<numOfVars; i++) if (!(v[i]==0)) return (0);
  return (1);
}
/* ----------------------------------------------------------------- */
int isAllOneVector(vector v, int numOfVars) {
  int i;

  for (i=0; i<numOfVars; i++) if (!(v[i]==1)) return (0);
  return (1);
}
/* ----------------------------------------------------------------- */
int isVectorEqualToVector(vector v, vector w, int numOfVars) {
  int i;

  if ((v==0) || (w==0)) return (0);
  for (i=0; i<numOfVars; i++) if (!(v[i]==w[i])) return (0);
  return (1);
}
/* ----------------------------------------------------------------- */
int isVectorEqualToNegativeVector(vector v, vector w, int numOfVars) {
  int i;
  if ((v==0) || (w==0)) return (0);
  for (i=0; i<numOfVars; i++) if (!(v[i]==-w[i])) return (0);
  return (1);
}
/* ----------------------------------------------------------------- */
vector copyVector(vector v, int numOfVars) {
  int i;
  vector w;

  w = createVector(numOfVars);
  for (i=0; i<numOfVars; i++) w[i] = v[i];
  return (w);
}
/* ----------------------------------------------------------------- */
vector negativeVector(vector v, int numOfVars) {
  int i;

  for (i=0; i<numOfVars; i++) v[i]=-v[i];
  return (v);
}
/* ----------------------------------------------------------------- */
vector lexPositiveVector(vector v, int numOfVars) {
  int i=0;

  while (i<numOfVars && v[i]==0) i++;
  if (v[i]<0) return(negativeVector(v, numOfVars));

  return (v);
}
/* ----------------------------------------------------------------- */
int isVectorLexPositive(vector v, int numOfVars) {
  int i=0;

  while (i<numOfVars && v[i]==0) i++;
  if (v[i]<0) return(0);

  return (1);
}
/* ----------------------------------------------------------------- */
vector addVector(vector v, vector w, int numOfVars) {
  int i;

  for (i=0; i<numOfVars; i++) v[i]=v[i]+w[i];
  return (v);
}
/* ----------------------------------------------------------------- */
vector addZeroOneVector(vector v, vector w, int numOfVars) {
  int i;

  for (i=0; i<numOfVars-1; i++) {
    if ((v[i] & w[i])!=0) {
      free(v);
      return(0);
    } else {
      v[i]=(v[i] | w[i]);
    }
  }
  v[numOfVars-1]=v[numOfVars-1]+w[numOfVars-1];

  return (v);
}
/* ----------------------------------------------------------------- */
vector subVector(vector v, vector w, int numOfVars) {
  int i;
  for (i=0; i<numOfVars; i++) v[i]=v[i]-w[i];
  return (v);
}
/* ----------------------------------------------------------------- */
listVector* createListVector(vector v) {
  return (appendVectorToListVector(v,0));
}
/* ----------------------------------------------------------------- */
int lengthListVector(listVector* LIST) {
  int len=0;

  while (LIST) {len++; LIST = LIST->rest;}
  return (len);
}
/* ----------------------------------------------------------------- */
listVector* copyListVector(listVector* LIST, int numOfVars) {
  vector v;
  listVector *tmp, *endtmp;

/* Copies LIST and vectors therein! */

  if (LIST==0) return(0);

  tmp=createListVector(copyVector(LIST->first,numOfVars));
  tmp->sign=LIST->sign;
  LIST=LIST->rest;
  endtmp=tmp;

  while (LIST) {
    v=copyVector(LIST->first,numOfVars);
    endtmp->rest=createListVector(v);
    endtmp=endtmp->rest;
    endtmp->sign=LIST->sign;
    LIST = LIST->rest;
  }
  return (tmp);
}
/* ----------------------------------------------------------------- */
listVector* copyListVectorWithoutVectors(listVector* LIST) {
  listVector *tmp, *endtmp;

/* Copies LIST and vectors therein! */

  if (LIST==0) return(0);

  tmp=createListVector(LIST->first);
  tmp->sign=LIST->sign;
  LIST=LIST->rest;
  endtmp=tmp;

  while (LIST) {
    endtmp->rest=createListVector(LIST->first);
    endtmp=endtmp->rest;
    endtmp->sign=LIST->sign;
    LIST = LIST->rest;
  }

  printf("len tmp = %d\n", lengthListVector(tmp));

  return (tmp);
}
/* ----------------------------------------------------------------- */
vector permuteVector(vector w, vector p, int numOfVars) {
  int i;
  vector v;

  if (p==0) return (w);
  if (w==0) return (w);
  v=createVector(numOfVars);
  for (i=0; i<numOfVars; i++) v[p[i]]=w[i];
  free(w);
  return (v);
}
/* ----------------------------------------------------------------- */
listVector* permuteListVector(listVector* LIST, vector p, int numOfVars) {
  vector v,tmp;
  listVector *tmp2;

  tmp2=LIST;
  v=createVector(numOfVars);
  while (LIST) {
    tmp=LIST->first;
    v=permuteVector(LIST->first,p,numOfVars);
    LIST->first=v;
    v=tmp;
    LIST = LIST->rest;
  }
  return (tmp2);
}
/* ----------------------------------------------------------------- */
vector rePermuteVector(vector v, vector w, vector p, int numOfVars) {
  int i;

  if (v==0) return (v);
  for (i=0; i<numOfVars; i++) v[i]=w[p[i]];
  return (v);
}
/* ----------------------------------------------------------------- */
listVector* rePermuteListVector(listVector* LIST, vector p, int numOfVars) {
  vector v,tmp;
  listVector *tmp2;

  tmp2=LIST;
  v=createVector(numOfVars);
  while (LIST) {
    tmp=LIST->first;
    v=rePermuteVector(v,LIST->first,p,numOfVars);
    LIST->first=v;
    v=tmp;
    LIST = LIST->rest;
  }
  return (tmp2);
}
/* ----------------------------------------------------------------- */
int isVectorInListVector(vector v, listVector* LIST, int numOfVars) {
  vector w;

  while (LIST) {
    w = LIST->first;
    LIST = LIST->rest;
    if (isVectorEqualToVector(v,w,numOfVars)==1) return (1);
  }
  return (0);
}
/* ----------------------------------------------------------------- */
int isNegativeVectorInListVector(vector v, listVector* LIST, int numOfVars) {
  vector w;

  if (LIST==0) return (0);
  while (LIST) {
    w = LIST->first;
    LIST = LIST->rest;
    if (isVectorEqualToNegativeVector(v,w,numOfVars)==1) return (1);
  }
  return (0);
}
/* ----------------------------------------------------------------- */
listVector* updateBasis(listVector *v, listVector *endBasis) {
  endBasis->rest = v;
  v->rest=0;
  endBasis = endBasis->rest;
  return (endBasis);
}
/* ----------------------------------------------------------------- */
int isBelowUpperBounds(vector v, vector upperBounds, int numOfVars) {
  int i;

  if (upperBounds==0) return (1);
  for (i=0; i<numOfVars; i++)
    if ((upperBounds[i]!=0) && (abs(v[i])>upperBounds[i])) return (0);
  return (1);
}
/* ----------------------------------------------------------------- */
listVector* appendListVectorToListVector(listVector *u, listVector *v) {
  listVector *w;
  /* Append u to v by forming a new listVector which points to the
     vectors in u and the end of listVector points to v. */
  while (u) {
    w=createListVector(u->first);
    w->rest=v;
    v=w;
    u=u->rest;
  }
  return (v);
}
/* ----------------------------------------------------------------- */
vector addMultipleVector(vector v, int a, vector w, int numOfVars) {
  int i;

  for (i=0; i<numOfVars; i++) v[i]=v[i]+a*w[i];
  return (v);
}
/* ----------------------------------------------------------------- */
vector subMultipleVector(vector v, int a, vector w, int numOfVars) {
  int i;

  for (i=0; i<numOfVars; i++) v[i]=v[i]-a*w[i];
  return (v);
}
/* ----------------------------------------------------------------- */
vector permuteMatrix(vector M, vector p, int numOfRows, int numOfVars) {
  int i,j;
  vector tmp;

  if (p==0) return(M);

  tmp=createVector(numOfRows*numOfVars);
  for (i=0; i<numOfRows; i++)
    for (j=0; j<numOfVars; j++)
      tmp[i*numOfVars+p[j]]=M[i*numOfVars+j];
  return (tmp);
}
/* ----------------------------------------------------------------- */
vector permuteTransposedMatrix(vector M, vector p, int numOfRows,
                   int numOfVars) {
  int i,j;
  vector tmp;

  if (p==0) return(M);
  if (M==0) return(M);

  tmp=createVector(numOfRows*numOfVars);
  for (i=0; i<numOfRows; i++)
    for (j=0; j<numOfVars; j++)
      tmp[(p[j])*numOfRows+i]=M[j*numOfRows+i];
  return (tmp);
}
/* ----------------------------------------------------------------- */
vector transpose(vector mat,int numOfVars,int numOfRows){
  int i,j,k,lenOfMatrix;
  vector transposedMat;

/*   printf("Transposing matrix\n"); */

  lenOfMatrix = numOfVars*numOfRows;
  transposedMat=createVector(lenOfMatrix);
  k=0;
  for (j=0; j<numOfVars; j++)
    for (i=0; i<numOfRows; i++)
      transposedMat[k++] = mat[i*numOfVars+j];
  return(transposedMat);
}
/* ----------------------------------------------------------------- */
vector matrixTimesVector(listVector *Mat, vector v, int numOfRows,
                         int numOfVars) {
  int i;
  vector b;

  b=createVector(numOfRows);
  for (i=0;i<numOfRows;i++) {
    b[i]=dotProduct(Mat->first,v,numOfVars);
    Mat=Mat->rest;
  }

  return (b);
}
/* ----------------------------------------------------------------- */
int hasSmallerSupport(vector u, vector v, int numOfVars) {
  int i;
  for (i=0; i<numOfVars; i++) if ((v[i]==0) && (u[i]!=0)) return (0);
  return(1);
}
/* ----------------------------------------------------------------- */
int isCircuit(listVector *basis, vector v, int numOfVars) {
  while (basis) {
    if ((hasSmallerSupport(basis->first,v,numOfVars)==1) && 
	(isVectorEqualToVector(basis->first,v,numOfVars)==0)) {
      return(0);
    }
    basis=basis->rest;
  }
  return(1);
}
/* ----------------------------------------------------------------- */
listVector* extractCircuits(listVector *basis, int numOfVars) {
  listVector* circuits=0;
  listVector* tmp;

  tmp=basis;
  while (tmp) {
    if (isCircuit(basis,tmp->first,numOfVars)==1) {
      if (circuits==0)
        circuits=createListVector(tmp->first);
      else
        circuits=appendVectorToListVector(tmp->first,circuits);
    }
    tmp=tmp->rest;
  }
  return(circuits);
}
/* ------------------------------------------------------------------------ */
int compareVectorsByLex(vector v, vector w, int numOfVars) {
  int i;
  /* Returns -1 if v is lex smaller, 0 if v=w, and 1 if w is smaller. */
  for (i=0; i<numOfVars; i++) 
    if (v[i]!=w[i]) {
      if (v[i]<w[i]) {
        return (-1);
      } else {
        return (1);
      }
    }
  return(0);
}
/* ----------------------------------------------------------------- */
listVector* combineOrderedListVectors(listVector *L1, listVector *L2, 
				      int numOfVars) {
  int r;
  listVector *tmp, *tmp2, *endtmp;
  
  if (L1==0) return (L2);
  if (L2==0) return (L1);

  r=compareVectorsByLex(L1->first,L2->first,numOfVars);

  tmp=0;
  endtmp=0;

  if (r==0) {
    tmp2=L2;
    L2=L2->rest;
    free(tmp2->first);
    free(tmp2);
    tmp=L1;
    endtmp=L1;
    L1=L1->rest;
    endtmp->rest=0;
  }
  if (r==-1) {
    tmp=L1;
    endtmp=L1;
    L1=L1->rest;
    endtmp->rest=0;
  }
  if (r==1) {
    tmp=L2;
    endtmp=L2;
    L2=L2->rest;
    endtmp->rest=0;
  }

  while ((L1!=0) && (L2!=0)) {
    r=compareVectorsByLex(L1->first,L2->first,numOfVars);
    if (r==0) {
      tmp2=L2;
      L2=L2->rest;
      free(tmp2->first);
      free(tmp2);
      endtmp->rest=L1;
      L1=L1->rest;
      endtmp=endtmp->rest;
      endtmp->rest=0;
    }
    if (r==-1) {
      endtmp->rest=L1;
      L1=L1->rest;
      endtmp=endtmp->rest;
      endtmp->rest=0;
    }
    if (r==1) {
      endtmp->rest=L2;
      L2=L2->rest;
      endtmp=endtmp->rest;
      endtmp->rest=0;
    }
  }

  if (L1!=0) endtmp->rest=L1;
  if (L2!=0) endtmp->rest=L2;

  return (tmp);
}
/* ----------------------------------------------------------------- */
listVector* swapColumnsInListVector(listVector* LIST, int a, int b,
				    int numOfVars) {
  int s;
  vector v;
  listVector *tmp;

/* printf("Swapping columns %d and %d.\n",a,b); */
  tmp=LIST;

  while (tmp) {
    v=tmp->first;
    s=v[a];
    v[a]=v[b];
    v[b]=s;
    tmp=tmp->rest;
  }
  return (LIST);
}
/* ----------------------------------------------------------------- */
listVector* extractVectorsBelowUpperBounds(listVector *basis, 
					   vector upperBounds, 
					   int numOfVars) {
  vector v;
  listVector *tmp, *tmp2, *endBasis;

  tmp=basis;
  basis=0;
  endBasis=0;

  while (tmp) {
    tmp2=tmp;
    v=tmp->first;
    if (isBelowUpperBounds(v,upperBounds,numOfVars)==1) {
      if (basis==0) {
	basis=createListVector(v);
	endBasis=basis;
      } else {
	endBasis->rest=createListVector(v);
	endBasis=endBasis->rest;
      }
    } else free(v);
    tmp=tmp->rest;
    free(tmp2);
  }

  return(basis);
}
/* ----------------------------------------------------------------- */
int normOfVector(vector v, int numOfVars) {
  int i, norm=0;

  for (i=0; i<numOfVars; i++) norm+=abs(v[i]);
  return (norm);
}
/* ----------------------------------------------------------------- */
int signedNormOfVector(vector v, int numOfVars) {
  int i, norm=0;

  for (i=0; i<numOfVars; i++) norm+=v[i];
  return (norm);
}
/* ----------------------------------------------------------------- */
int weightedNormOfVector(vector v, vector c, int numOfVars) {
  int i, norm=0;

  for (i=0; i<numOfVars; i++) norm+=abs(c[i]*v[i]);
  return (norm);
}
/* ----------------------------------------------------------------- */
int minimalNormInListVector(listVector* L, int numOfVars) {
  int s,norm;

  norm=-1;
  while (L) {
    s=normOfVector(L->first,numOfVars);   
    if (norm==-1) {
      norm=s;
    } else {
      if (s<norm) norm=s;
    }
    L=L->rest;
  }

  return (norm);
}
/* ----------------------------------------------------------------- */
int maximalNormInListVector(listVector* L, int numOfVars) {
  int s,norm;

  norm=-1;
  while (L) {
    s=normOfVector(L->first,numOfVars);   
    if (norm==-1) {
      norm=s;
    } else {
      if (s>norm) norm=s;
    }
    L=L->rest;
  }

  return (norm);
}
/* ----------------------------------------------------------------- */
listVector* projectListVectorDown(listVector* basis, int numOfVars, 
				  int dimension) {
  int i;
  vector v;
  listVector *tmp, *endtmp;

/*   printf("Projecting vectors down to dimension %d.\n",dimension); */

  if (basis==0) return (0);

  v=createVector(dimension);
  for (i=0; i<dimension; i++) v[i]=(basis->first)[i];
  tmp=createListVector(v);
  endtmp=tmp;
  basis=basis->rest;

  while (basis) {
    v=createVector(dimension);
    for (i=0; i<dimension; i++) v[i]=(basis->first)[i];
    endtmp=updateBasis(createListVector(v),endtmp);
    basis=basis->rest;
  }
  return (tmp);
}
/* ----------------------------------------------------------------- */
listVector* reduceLastComponent(listVector* basis, vector v, 
				int numOfVars) {
  int a;
  listVector* tmp;

  tmp=basis;
  while (tmp) {
    a=((tmp->first)[numOfVars-1])/v[numOfVars-1];
    if ((tmp->first)[numOfVars-1]<0) {
      if (a*v[numOfVars-1]==(tmp->first)[numOfVars-1]) a++;
      a--;
    }
    tmp->first=subMultipleVector(tmp->first,a,v,numOfVars);
    tmp=tmp->rest;
  }
  return(basis);
}
/* ----------------------------------------------------------------- */
listVector* appendNewComponentToListVector(listVector* basis, 
					   listVector* originalInput, 
					   int currentVars,int dimension,
					   int infoLevel) {
  int i,j,a;
  vector v,w;
  listVector *tmp, *endtmp, *tmp2;

  setbuf(stdout,0);

  /* We heavily rely on a triangularized originalInput!!! */

if (infoLevel>0) {
  printf("Appending component %d.\n",currentVars);
}
  if (basis==0) return (0);

  v=createVector(currentVars);
  for (i=0; i<currentVars-1; i++) v[i]=(basis->first)[i];
  v[currentVars-1]=0;

  tmp2=originalInput;
  w=copyVector(v,currentVars);
  j=0;
  while (tmp2) {
    if (j<currentVars-1) {
      a=w[j]/((tmp2->first)[j]);
if (w[j]!=a*((tmp2->first)[j])) {
  printf("Foul play in appendNewComponentToListVector!\n");
  printVector(w,currentVars);
  printVector(tmp2->first,currentVars);
  printf("a = %d\n",a);
  exit(0);
}
      v[currentVars-1]=v[currentVars-1]+(tmp2->first)[currentVars-1]*a;
      w=subMultipleVector(w,a,tmp2->first,currentVars);
    }
    tmp2=tmp2->rest;
    j++;
  }
  freeVector(w);
  tmp=createListVector(v);
  endtmp=tmp;

  tmp2=basis;
  basis=basis->rest;
  freeListVector(tmp2);

  while (basis) {
    v=createVector(currentVars);
    for (i=0; i<currentVars-1; i++) v[i]=(basis->first)[i];
    v[currentVars-1]=0;
    tmp2=originalInput;
    w=copyVector(v,currentVars);
    j=0;
    while (tmp2) {
      if (j<currentVars-1) {
	a=w[j]/((tmp2->first)[j]);
if (w[j]!=a*((tmp2->first)[j])) {
  printf("Foul play in appendNewComponentToListVector!\n");
  printVector(w,currentVars);
  printVector(tmp2->first,currentVars);
  printf("a = %d\n",a);
  exit(0);
}
        v[currentVars-1]+=(tmp2->first)[currentVars-1]*a;
	w=subMultipleVector(w,a,tmp2->first,currentVars);
      }
      tmp2=tmp2->rest;
      j++;
    }
    freeVector(w);
    endtmp=updateBasis(createListVector(v),endtmp);
    tmp2=basis;
    basis=basis->rest;
    freeListVector(tmp2);
  }
  return (tmp);
}
/* ----------------------------------------------------------------- */
listVector* appendRemainingComponentsToListVector(listVector* basis, 
						  listVector* originalInput, 
						  int numOfVars,
						  int numOfAllVars) {
  int i,j,a;
  vector v,w;
  listVector *tmp, *endtmp, *tmp2;

  /* We heavily rely on a triangularized originalInput!!! */

  if (basis==0) return (0);

  v=createVector(numOfAllVars);
  for (i=0; i<numOfVars; i++) v[i]=(basis->first)[i];
  for (i=numOfVars; i<numOfAllVars; i++) v[i]=0;

  tmp2=originalInput;
  w=copyVector(v,numOfVars);
  j=0;

  while (tmp2) {
    a=w[j]/((tmp2->first)[j]);
    for (i=numOfVars; i<numOfAllVars; i++)
      v[i]=v[i]+(tmp2->first)[i]*a;
    w=subMultipleVector(w,a,tmp2->first,numOfVars);
    tmp2=tmp2->rest;
    j++;
  }
  freeVector(w);
  tmp=createListVector(v);
  endtmp=tmp;

  tmp2=basis;
  basis=basis->rest;
  freeListVector(tmp2);

  while (basis) {
    v=createVector(numOfAllVars);
    for (i=0; i<numOfVars; i++) v[i]=(basis->first)[i];
    for (i=numOfVars; i<numOfAllVars; i++) v[i]=0;
    tmp2=originalInput;
    w=copyVector(v,numOfVars);
    j=0;

    while (tmp2) {
      a=w[j]/((tmp2->first)[j]);
      for (i=numOfVars; i<numOfAllVars; i++)
	v[i]=v[i]+(tmp2->first)[i]*a;
      w=subMultipleVector(w,a,tmp2->first,numOfVars);
      tmp2=tmp2->rest;
      j++;
    }
    freeVector(w);
    endtmp=updateBasis(createListVector(v),endtmp);
    tmp2=basis;
    basis=basis->rest;
    freeListVector(tmp2);
  }
  return (tmp);
}
/* ------------------------------------------------------------------------ */
int isZeroOneVector(vector v, int numOfVars) {
  int i;
  for (i=0; i<numOfVars; i++) if (abs(v[i])>1) return (0);
  return (1);
}
/* ----------------------------------------------------------------- */
listVector* extractZeroOneVectors(listVector *basis, int numOfVars) {
  int s;
  vector v;
  listVector *tmp, *endBasis;

  tmp=basis;
  basis=0;
  endBasis=0;

  while (tmp) {
    v=tmp->first;
    s=isZeroOneVector(v,numOfVars);
    if (s==1) {
      if (basis==0) {
	basis=createListVector(v);
	endBasis=basis;
      } else {
	endBasis->rest=createListVector(v);
	endBasis=endBasis->rest;
      }
    } else free(v);
    tmp=tmp->rest;
  }

  return(basis);
}
/* ----------------------------------------------------------------- */
listVector* extractZeroOneVectorsInLastComponent(listVector *basis, 
						 int numOfVars) {
  vector v;
  listVector *tmp, *endBasis;

  tmp=basis;
  basis=0;
  endBasis=0;

  while (tmp) {
    v=tmp->first;
    if (abs(v[numOfVars-1])>1) {
      if (basis==0) {
	basis=createListVector(v);
	endBasis=basis;
      } else {
	endBasis->rest=createListVector(v);
	endBasis=endBasis->rest;
      }
    } else free(v);
    tmp=tmp->rest;
  }

  return(basis);
}
/* ------------------------------------------------------------------------ */
int isNonnegativeVector(vector v, int numOfVars) {
  int i;
  for (i=0; i<numOfVars; i++) if (v[i]<0) return (0);
  return (1);
}
/* ----------------------------------------------------------------- */
listVector* extractNonnegativeVectors(listVector *basis, int numOfVars) {
  int s;
  vector v;
  listVector *tmp, *endBasis;

  tmp=basis;
  basis=0;
  endBasis=0;

  while (tmp) {
    v=tmp->first;
    s=isNonnegativeVector(v,numOfVars);
    if (s==1) {
      if (basis==0) {
	basis=createListVector(v);
	endBasis=basis;
      } else {
	endBasis->rest=createListVector(v);
	endBasis=endBasis->rest;
      }
    } else free(v);
    tmp=tmp->rest;
  }

  return(basis);
}
/* ----------------------------------------------------------------- */
void swapGraver(vector *B, int i, int j) {
  /* swap pointers in B[i] and in B[j], using B[0] as a temporary storage */
  B[0] = B[i];
  B[i] = B[j];
  B[j] = B[0];
  return;
}
/* ------------------------------------------------------------------------ */
void liftGraver(vector *B, vector *P, int m, int numOfVars) {
  /* lift B[m] up, at most up to B[1] */
  int j;
  j=m;
  while (1<j) {
    if (compareVectorsByLex(B[j/2],B[j],numOfVars)==1) {
      swapGraver(B,j/2,j);
      swapGraver(P,j/2,j);
      j=j/2;
    } else return;
  }
  return;
}
/* ------------------------------------------------------------------------ */
int removeFirstHeapElement(vector *B, vector *P, int m, int numOfVars) {
  int i,j;

  free(B[1]);
  B[1]=B[m];
  P[1]=P[m];
  m--;
  i=1;

  while (2*i <= m) {
    j = 2*i;
    if (j < m) {   /* take the smaller one of B[j], B[j+1] */
      if (compareVectorsByLex(B[j],B[j+1],numOfVars)==1)
        j = j+1;
    }
    if (compareVectorsByLex(B[i],B[j],numOfVars)==1) {
      swapGraver(B,i,j);
      swapGraver(P,i,j);
      i = j;
    } else return (m);
  }

  return (m);
}
/* ----------------------------------------------------------------- */
int updateHeapGraver(vector v, vector *B, vector w, vector *P,
		     int lengthOfHeap, int numOfVars) {
  lengthOfHeap++;
  B[lengthOfHeap]=v;
  P[lengthOfHeap]=w;

  if (lengthOfHeap>1) liftGraver(B,P,lengthOfHeap,numOfVars);

  return (lengthOfHeap);
}
/* ----------------------------------------------------------------- */
int isInSameOrthant(vector v, vector w, int numOfVars) {
  int i;

  if ((v==0) || (w==0)) return (1);
  for (i=0; i<numOfVars; i++) {
    if ((v[i]>0) && (w[i]<0)) return (0);
    if ((v[i]<0) && (w[i]>0)) return (0);
  }
  return (1);
}
/* ----------------------------------------------------------------- */
int isInOppositeOrthant(vector v, vector w, int numOfVars) {
  int i;

  if ((v==0) || (w==0)) return (1);
  for (i=0; i<numOfVars; i++) {
    if ((v[i]>0) && (w[i]>0)) return (0);
    if ((v[i]<0) && (w[i]<0)) return (0);
  }
  return (1);
}
/* ----------------------------------------------------------------- */
int isIdentityPermutation(vector p, int numOfVars) {
  int i;

  if (p==0) return (0);
  for (i=0; i<numOfVars; i++) if ((p[i])!=i) return (0);
  return (1);
}
/* ----------------------------------------------------------------- */
int isVectorEqualToPermutedVector(vector v, vector p, int numOfVars) {
  int i;

  if ((v==0) || (p==0)) return (0);
  for (i=0; i<numOfVars; i++) if (!(v[i]==v[p[i]])) return (0);
  return (1);
}
/* ----------------------------------------------------------------- */
int isVectorEqualToNegativePermutedVector(vector v, vector p, int numOfVars) {
  int i;
  if ((v==0) || (p==0)) return (0);
  for (i=0; i<numOfVars; i++) if (!(v[i]==-v[p[i]])) return (0);
  return (1);
}
/* ----------------------------------------------------------------- */
vector subMultiplePermutedVector(vector v, int a, vector w, vector p,
				 int numOfVars) {
  int i;

  for (i=0; i<numOfVars; i++) v[i]=v[i]-a*w[p[i]];
  return (v);
}
/* ----------------------------------------------------------------- */
vector multiplyPermutation(vector v, vector w, int numOfVars) {
  int i;

/* Computes permutation w(v(x)). */

  for(i=0; i<numOfVars; i++) v[i]=w[v[i]];

  return (v);
}
/* ------------------------------------------------------------------------ */
listVector* vTimesS(vector v, listVector *S, int signV, int numOfVars) {
  vector w;
  listVector *tmp, *L, *endL;

/* Multiplies permutation v to all permutations in S. */

  L=createListVector(0);
  endL=L;

  tmp=S;
  while (tmp) {
    w=copyVector(tmp->first,numOfVars);
    w=multiplyPermutation(w,v,numOfVars);
    w=lexPositiveVector(w,numOfVars);
    if (isVectorInListVector(w,L->rest,numOfVars)==0) {
      endL->rest=createListVector(w);
      endL=endL->rest;
      endL->sign=signV*(tmp->sign);
    } else free(w);
    tmp=tmp->rest;
  }

  return (L->rest);
}
/* ------------------------------------------------------------------------ */
listVector* generateSymmetryGroup(listVector *S, int numOfVars) {
  listVector *G, *endG, *tmp, *L;

printf("Generating symmetry group.\n");
  G=copyListVector(S,numOfVars);
  endG=G;
  tmp=G;
  while (tmp) {
    endG=tmp;
    tmp=tmp->rest;
  }

  tmp=G;
  while (tmp) {
    L=vTimesS(tmp->first,S,tmp->sign,numOfVars);
    while (L) {
      if (isVectorInListVector(L->first,G,numOfVars)==0) {
	endG->rest=createListVector(L->first);
	endG=endG->rest;
	endG->sign=L->sign;
      }
      L=L->rest;
    }
    tmp=tmp->rest;
  }

  printf("Symmetry group has been computed.\n");

  return (G);
}
/* ------------------------------------------------------------------------ */
vector orientVector(vector v, int variable, int numOfVars) {
  int i;

/* Orients the vector v according to degrevlex with "variable" as
   smallest variable. */

/*   if (signedNormOfVector(v,numOfVars)>0) return (v); */
/*   if (signedNormOfVector(v,numOfVars)<0) return (negativeVector(v,numOfVars)); */

  if (v[variable]<0) {
    return (v);
  }
  if (v[variable]>0) {
    v=negativeVector(v,numOfVars);
    return (v);
  }

  for (i=0; i<numOfVars; i++) {
    if (i!=variable) {
      if (v[i]<0) return (v);
      if (v[i]>0) return (negativeVector(v,numOfVars));
    }
  }

  printVector(v,numOfVars);
  printf("I should never have ended up here! (orientVector)\n");
  exit(0);

  return (v);
}
/* ----------------------------------------------------------------- */
vector orientVectorAccordingToOrdering(vector v, vector costVector, 
				       int numOfVars) {
  int s;

/* Orients the vector v according to degrevlex with "variable" as
   smallest variable. */

/*   if (signedNormOfVector(v,numOfVars)>0) return (v); */
/*   if (signedNormOfVector(v,numOfVars)<0) return (negativeVector(v,numOfVars)); */

  s=dotProduct(v,costVector,numOfVars);

  if (s>0) return (v);
  if (s<0) return (negativeVector(v,numOfVars));

  return (orientVector(v,0,numOfVars));
}
/* ----------------------------------------------------------------- */
int isVectorCorrectlyOriented(vector v, int variable, int numOfVars) {
  int i;

/* Tests whether the vector v is correctly oriented according to 
   degrevlex with "variable" as smallest variable. */

/*   if (signedNormOfVector(v,numOfVars)>0) return (v); */
/*   if (signedNormOfVector(v,numOfVars)<0) return (negativeVector(v,numOfVars)); */

  if (v[variable]<0) return (1);
  if (v[variable]>0) return (0);

  for (i=0; i<numOfVars; i++) {
    if (i!=variable) {
      if (v[i]<0) return (1);
      if (v[i]>0) return (0);
    }
  }

  printf("I should never have ended up here! (isVectorCorrectlyOriented)\n");
  exit(0);

  return (1);
}
/* ------------------------------------------------------------------------ */
int isZeroVectorInListVector(listVector *L, int numOfVars) {
  while (L) {
    if (L->first!=0) {
      if (isZeroVector(L->first,numOfVars)==1) return (1);
    }
    L=L->rest;
  }

  return(0);
}
/* ----------------------------------------------------------------- */
