/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.

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
/* Print functions                                                   */
/*                                                                   */
/* Author   : Raymond Hemmecke   (implementation)                    */
/* Co-Author: Ralf Hemmecke      (data structure, code optimization) */
/*                                                                   */
/*------------------------------------------------------------------ */
#include "myheader.h"
#include "orbit.h"
#include "vector.h"
#include <stdlib.h>
/* ----------------------------------------------------------------- */
void printVersionInfo() {

  printf("-------------------------------------------------\n");
  printf("4ti2 version 1.3.1, Copyright (C) 2006 4ti2 team.\n");
  printf("4ti2 comes with ABSOLUTELY NO WARRANTY.\n");
  printf("This is free software, and you are welcome\n");
  printf("to redistribute it under certain conditions.\n");
  printf("For details, see the file COPYING.\n");
  printf("-------------------------------------------------\n");
  return ;
}
/* ----------------------------------------------------------------- */
void printVector(vector v, int numOfVars) {
  int i;

  /*  lexPositiveVector(v,numOfVars);*/

  if (v==0) {
    printf("[]\n");
    return ;
  }
  printf("%s", "[");
  for (i=0; i<(numOfVars-1); i++) {
    printf("%d ", v[i]);
  }
  printf("%d]\n", v[i]);
  return ;
}
/* ----------------------------------------------------------------- */
void printListVector(listVector* basis, int numOfVars) {
  if (basis==0) printf("[]\n");
  while(basis) {
    printVector(basis->first,numOfVars);
    basis = basis->rest;
  }
  printf("\n");
  return ;
}
/* ----------------------------------------------------------------- */
void printVectorToFile(FILE *out, vector v, int numOfVars) {
  int i,multiplier=1;

  if (v==0) return;
  if (isVectorLexPositive(v,numOfVars)==0) multiplier=-1;

  multiplier=1;
  for (i=0; i<(numOfVars); i++) fprintf(out,"%d ", multiplier*v[i]);
  fprintf(out,"\n");
  return ;
}
/* ----------------------------------------------------------------- */
void printListVectorToFile(char* fileName, listVector* basis, int numOfVars) {
  int len;
  FILE* out;

  if (!(out = fopen(fileName,"w"))) {
    printf("Error opening output file!");
    exit (0);
  }
  if (basis==0) {
    fprintf(out,"0 %d\n",numOfVars);
    fclose(out);
    return;
  }

  len=lengthListVector(basis);
  fprintf(out,"%d %d\n",len,numOfVars);
  while(basis) {
    printVectorToFile(out,basis->first,numOfVars);
    basis = basis->rest;
  }
  fprintf(out,"\n");
  fclose(out);
  return ;
}
/* ----------------------------------------------------------------- */
void printTransposedListVectorToFile(char* fileName, listVector* basis, 
				     int numOfVars) {
  int i,len;
  listVector *tmp;
  FILE* out;

  if (!(out = fopen(fileName,"w"))) {
    printf("Error opening output file!");
    exit (0);
  }
  if (basis==0) {
    fclose(out);
    return;
  }

  len=lengthListVector(basis);
  fprintf(out,"%d %d\n",numOfVars,len);

  for (i=0; i<numOfVars; i++) {
    tmp=basis;
    while(tmp) {
      fprintf(out,"%d ",(tmp->first)[i]);
      tmp=tmp->rest;
    }
    fprintf(out,"\n");
  }
  fprintf(out,"\n");
  fclose(out);

  return ;
}
/* ----------------------------------------------------------------- */
void printVectorToFileWithBrackets(FILE *out, vector v, int numOfVars) {
  int i,multiplier=1;

  if (isVectorLexPositive(v,numOfVars)==0) multiplier=-1;

  multiplier=1;
  fprintf(out,"[");
  for (i=0; i<(numOfVars-1); i++) fprintf(out,"%d ", multiplier*v[i]);
  fprintf(out,"%d]\n", multiplier*v[numOfVars-1]);
  return ;
}
/* ----------------------------------------------------------------- */
void printBinomialToFile(FILE *out, vector v, int numOfVars, char** labels) {
  int i,tmp,posNorm=0,negNorm=0;

  for (i=0; i<numOfVars; i++) {
    if (v[i]>0) posNorm = posNorm + v[i];
    else if (v[i]<0) negNorm = negNorm - v[i];
  }

  if (posNorm==0) fprintf(out,"1-");
  else {
    tmp=posNorm;
    i=0;
    while (i<numOfVars && tmp>0) {
      if (v[i]>0) {
	if (v[i]>1) 
	  if (labels!=0) { 
	    fprintf(out,"%s^%d",labels[i],v[i]);
	  } else {
	    fprintf(out,"x[%d]^%d",i+1,v[i]);
	  }
	else
	  if (labels!=0) { 
	    fprintf(out,"%s",labels[i]);
	  } else {
	    fprintf(out,"x[%d]",i+1);
	  }
	tmp=tmp-v[i];
	if (tmp>0) fprintf(out,"*");
      }
      i++;
    }
    fprintf(out,"-");
  }
  if (negNorm==0) fprintf(out,"1");
  else {
    tmp=negNorm;
    i=0;
    while (i<numOfVars && tmp>0) {
      if (v[i]<0) {
	if (v[i]<-1) 
	  if (labels!=0) { 
	    fprintf(out,"%s^%d",labels[i],-v[i]);
	  } else {
	    fprintf(out,"x[%d]^%d",i+1,-v[i]);
	  }
	else
	  if (labels!=0) { 
	    fprintf(out,"%s",labels[i]);
	  } else {
	    fprintf(out,"x[%d]",i+1);
	  }
        tmp=tmp+v[i];
        if (tmp>0) fprintf(out,"*");
      }
      i++;
    }
  }
  return ;
}
/* ----------------------------------------------------------------- */
void printListBinomialsToFile(char* fileName, listVector* basis, 
			      int numOfVars, char** labels) {
  FILE* out;

  if (!(out = fopen(fileName,"w"))) {
    printf("Error opening binomial file!");
    exit (0);
  }
  if (basis==0) {
    fprintf(out,"[]\n");
    fclose(out);
    return;
  }

  printf("Writing binomials to file\n\n");
  fprintf(out,"[\n");
  while(basis->rest) {
    printBinomialToFile(out,basis->first,numOfVars,labels);
    fprintf(out,",\n");
    basis = basis->rest;
  }
  printBinomialToFile(out,basis->first,numOfVars,labels);
  fprintf(out,"\n]\n");
  fclose(out);
  return ;
}
/* ----------------------------------------------------------------- */
void printMonomialToFile(FILE *out, vector v, int numOfVars, char** labels) {
  int i,tmp,posNorm=0;

  for (i=0; i<numOfVars; i++) posNorm = posNorm + v[i];

  if (posNorm==0) fprintf(out,"1");
  else {
    tmp=posNorm;
    i=0;
    while (i<numOfVars && tmp>0) {
      if (v[i]>0) {
	if (v[i]>1) 
	  if (labels!=0) { 
	    fprintf(out,"%s^%d",labels[i],v[i]);
	  } else {
	    fprintf(out,"x[%d]^%d",i+1,v[i]);
	  }
	else
	  if (labels!=0) { 
	    fprintf(out,"%s",labels[i]);
	  } else {
	    fprintf(out,"x[%d]",i+1);
	  }
	tmp=tmp-v[i];
	if (tmp>0) fprintf(out,"*");
      }
      i++;
    }
  }
  return ;
}
/* ----------------------------------------------------------------- */
void printListMonomialsAndBinomialsToFile(char* fileName, listVector* basis, 
					  int numOfVars, char** labels) {
  FILE* out;

  if (!(out = fopen(fileName,"w"))) {
    printf("Error opening binomial file!");
    exit (0);
  }
  if (basis==0) {
    fprintf(out,"[]\n");
    fclose(out);
    return;
  }

  printf("Writing monomials and binomials to file\n\n");
  fprintf(out,"[\n");
  while(basis->rest) {
    if (isNonnegativeVector(basis->first,numOfVars)==1) {
      printMonomialToFile(out,basis->first,numOfVars,labels);
    } else {
      printBinomialToFile(out,basis->first,numOfVars,labels);
    }
    fprintf(out,",\n");
    basis = basis->rest;
  }
  printBinomialToFile(out,basis->first,numOfVars,labels);
  fprintf(out,"\n]\n");
  fclose(out);
  return ;
}
/* ----------------------------------------------------------------- */
void printMatrix(vector v, int numOfRows, int numOfVars) {
  int i,j;

  for (j=0; j<(numOfRows); j++) {
    printf("%s", "[");
    for (i=0; i<(numOfVars-1); i++) {
      printf("%d ", v[i]);
    }
    printf("%d]\n", v[i]);
    v+=numOfVars;
  }
  printf("\n");
  return ;
}
/* ----------------------------------------------------------------- */
void printVectorToFileMaple(FILE *out, vector v, int numOfVars) {
  int i;

  if (v==0) {
    fprintf(out,"[]\n");
    return ;
  }

  fprintf(out,"%s", "[");
  for (i=0; i<(numOfVars-1); i++) {
    fprintf(out,"%d,", v[i]);
  }
  fprintf(out,"%d]", v[i]);
  return ;
}
/* ----------------------------------------------------------------- */
void printVectorToFileMacaulay2(FILE *out, vector v, int numOfVars) {
  int i;

  if (v==0) {
    fprintf(out,"{}\n");
    return ;
  }

  fprintf(out,"%s", "{");
  for (i=0; i<(numOfVars-1); i++) {
    fprintf(out,"%d,", v[i]);
  }
  fprintf(out,"%d}", v[i]);
  return ;
}
/* ----------------------------------------------------------------- */
void printListVectorMaple(char* fileName, listVector* basis, int numOfVars) {
  FILE* out;

  if (!(out = fopen(fileName,"w"))) {
    printf("Error opening output file");
    exit (0);
  }
  if (basis==0) {
    fprintf(out,"[]\n");
    fclose(out);
    return;
  }

  fprintf(out,"[");
  while(basis->rest) {
    printVectorToFileMaple(out,basis->first,numOfVars);
    basis = basis->rest;
    fprintf(out,",");
  }
  printVectorToFileMaple(out,basis->first,numOfVars);
  fprintf(out,"]\n\n");
  fclose(out);
  return ;
}
/* ----------------------------------------------------------------- */
void printListVectorMacaulay2(char* fileName, listVector* basis, 
			      int numOfVars) {
  FILE* out;

  if (!(out = fopen(fileName,"w"))) {
    printf("Error opening output file");
    exit (0);
  }
  if (basis==0) {
    fprintf(out,"{}\n");
    fclose(out);
    return;
  }

  fprintf(out,"{");
  while (basis->rest) {
    printVectorToFileMacaulay2(out,basis->first,numOfVars);
    basis = basis->rest;
    fprintf(out,",");
  }
  printVectorToFileMacaulay2(out,basis->first,numOfVars);
  fprintf(out,"}\n\n");
  fclose(out);
  return ;
}
/* ----------------------------------------------------------------- */
void print3wayTables(char* fileName, listVector* basis, int x, int y, 
		     int z, int numOfVars) {
  int i,j,k;
  vector v;
  FILE* out;

  if (!(out = fopen(fileName,"w"))) {
    printf("Error opening output file");
    exit (0);
  }

  fprintf(out,"%d %d %d %d\n",x,y,z,lengthListVector(basis));

  if (basis==0) {
    fprintf(out,"[]\n");
    fclose(out);
    return;
  }

  while (basis) {
    fprintf(out,"===============\n");
    v=basis->first;
    printVectorToFile(out,v,numOfVars);
    fprintf(out,"===============\n");
    for (i=0; i<z; i++) {
      for (j=0; j<y; j++) {
	for (k=0; k<x; k++) {
	  if (v[i*x*y+j*x+k]<0) {
	    fprintf(out,"%d ",v[i*x*y+j*x+k]);
	  } else {
	    fprintf(out," %d ",v[i*x*y+j*x+k]);
	  }
	}
	fprintf(out,"\n");
      }
      if (i<z-1) fprintf(out,"\n");
    } 
    basis = basis->rest;
  }
  fprintf(out,"===============\n");
  fclose(out);
  return ;
}
/* ----------------------------------------------------------------- */
void printL1NormOfListVector(listVector *basis, int numOfVars) {
  int i,s,maxs;
  vector sum;
  listVector* tmp;

  sum=createVector(100000);
  maxs=0;

  for (i=0; i<100000; i++) sum[i]=0;  

  tmp=basis;
  while (tmp) {
    s=0;
    for (i=0; i<numOfVars; i++) s=s+abs((tmp->first)[i]);  
    sum[s]=sum[s]+1;
    if (s>maxs) maxs=s;
    tmp=tmp->rest;
  }

/*   tmp=basis; */
/*   while (tmp) { */
/*     s=0; */
/*     for (i=0; i<numOfVars; i++) s=s+abs((tmp->first)[i]);   */
/*     if (s==maxs) printVector(tmp->first,numOfVars); */
/*     tmp=tmp->rest; */
/*   } */

  for (i=0; i<100000; i++) 
    if (sum[i]>0) printf("Norm = %d,   number of elements = %d\n",i,sum[i]); 

  return;
}
/* ----------------------------------------------------------------- */
void printListVectorWithGivenDegreeToFile(char *outFileName, 
					  listVector *basis, int numOfVars,
					  int degree){
  int i,s,len;
  listVector *tmp;
  FILE* out;

  if (!(out = fopen(outFileName,"w"))) {
    printf("Error opening output file!");
    exit (0);
  }
  if (basis==0) {
    fprintf(out,"0 %d\n",numOfVars);
    fclose(out);
    return;
  }

  len=0;
  tmp=basis;
  while (tmp) {
    s=0;
    for (i=0; i<numOfVars; i++) s=s+abs((tmp->first)[i]);
    if (s==degree) len++;
    tmp=tmp->rest;
  }
  fprintf(out,"%d %d\n",len,numOfVars);

  tmp=basis;
  while (tmp) {
    s=0;
    for (i=0; i<numOfVars; i++) s=s+abs((tmp->first)[i]);
    if (s==degree) printVectorToFile(out,tmp->first,numOfVars);
    tmp=tmp->rest;
  }
  fclose(out);

  return;
}
/* ----------------------------------------------------------------- */
void printListVectorWithGivenNonzeroEntryToFile(char *outFileName, 
						listVector *basis, 
						int numOfVars,
						int coord){
  int len;
  listVector *tmp;
  FILE* out;

  if (!(out = fopen(outFileName,"w"))) {
    printf("Error opening output file!");
    exit (0);
  }
  if (basis==0) {
    fprintf(out,"0 %d\n",numOfVars);
    fclose(out);
    return;
  }

  len=0;
  tmp=basis;
  while (tmp) {
    if ((tmp->first)[coord-1]!=0) len++;
    tmp=tmp->rest;
  }
  fprintf(out,"%d %d\n",len,numOfVars);

  tmp=basis;
  while (tmp) {
    if ((tmp->first)[coord-1]!=0) printVectorToFile(out,tmp->first,numOfVars);
    tmp=tmp->rest;
  }
  fclose(out);

  return;
}
/* ----------------------------------------------------------------- */
void writeResult(listVector *basis, int numOfVars, char *fileName, 
		 char *basisType, int infoLevel) {
  char outFileName[127];

  /* Write result to screen and files. */

if (infoLevel>0) {
  printf("Writing result to files: ");
}
  /* Write Graver basis file. */
  if (basisType[0]=='g') {
if (infoLevel>0) {
    printf("Graver basis elements: %d\n\n",lengthListVector(basis));
}
    strcpy(outFileName,fileName);
    printListVectorToFile(outFileName,basis,numOfVars);
  } 

  /* Write Hilbert basis file. */
  if (basisType[0]=='h') {
if (infoLevel>0) {
    printf("Hilbert basis elements: %d\n\n",lengthListVector(basis));
}

    strcpy(outFileName,fileName);
    printListVectorToFile(outFileName,basis,numOfVars);
  }

  /* Write Hilbert basis file of dual cone. */
  if (basisType[0]=='d') {
if (infoLevel>0) {
    printf("Hilbert basis elements: %d\n\n",lengthListVector(basis));
}

    strcpy(outFileName,fileName);
    strcat(outFileName,".dual.hil");
    printListVectorToFile(outFileName,basis,numOfVars);
  }

  /* Write extreme rays file. */
  if (basisType[0]=='r') {
if (infoLevel>0) {
    printf("Extreme rays: %d\n\n",lengthListVector(basis));
}
    strcpy(outFileName,fileName);
    strcat(outFileName,".ray");
    printListVectorToFile(outFileName,basis,numOfVars);
  }

  return;
}
/* ----------------------------------------------------------------- */
void printListRepresentativesToFile(char* fileName, listOrbit* basis, 
				    int numOfVars) {
  int len;
  FILE* out;

  if (!(out = fopen(fileName,"w"))) {
    printf("Error opening output file!");
    exit (0);
  }
  if (basis==0) {
    fclose(out);
    return;
  }

  len=lengthListOrbit(basis);
  fprintf(out,"%d %d\n",len,numOfVars);
  while(basis) {
    printVectorToFile(out,(basis->first)->representative,numOfVars);
    basis = basis->rest;
  }
  fprintf(out,"\n");
  fclose(out);
  return ;
}
/* ----------------------------------------------------------------- */
