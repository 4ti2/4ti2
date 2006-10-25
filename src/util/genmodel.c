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
/* Generation of toric models in statistics                          */
/*                                                                   */
/* Author   : Raymond Hemmecke                                       */
/*                                                                   */
/* Created    : 01-JUL-03                                            */
/*                                                                   */
/* ----------------------------------------------------------------- */
#include <stdio.h>
#include <string.h>
#include "myheader.h"
#include "vector.h"
#include "print.h"
/* ----------------------------------------------------------------- */
listVector* readSimplicialComplex(char* fileName, int* numOfNodes) {
  int i,j,numOfFaces,sizeOfFace;
  listVector *basis, *endBasis;
  vector b;
  FILE *in;

  setbuf(stdout,0);
  if (!(in = fopen(fileName,"r"))) {
    printf("Error opening file %s containing the simplicial complex.\n",
	   fileName);
    exit(0);
  }

  fscanf(in,"%d",numOfNodes);

  b=createVector(*numOfNodes);
  for (j=0; j<(*numOfNodes); j++) fscanf(in,"%d",&b[j]);
  basis = createListVector(b);
  endBasis = basis;

  fscanf(in,"%d",&numOfFaces);

  for (i=0; i<numOfFaces; i++) {
    fscanf(in,"%d",&sizeOfFace);
    b=createVector(sizeOfFace);
    b[0]=sizeOfFace;
    for (j=0; j<sizeOfFace; j++) fscanf(in,"%d",&b[j+1]);
    endBasis = updateBasis(createListVector(b), endBasis);
  }
  fclose(in);

  return (basis);
}
/* ----------------------------------------------------------------- */
vector decomposeIntegerIntoLevelIndices(int z, int numOfVars, vector face, 
					vector levels) {
  int i;
  vector v;

  v=createVector(numOfVars);
  for(i=0;i<numOfVars;i++) {
    v[i]=z-(levels[face[i+1]-1])*(z/(levels[face[i+1]-1]));
    z=(z-v[i])/(levels[face[i+1]-1]);
  }

  return(v);
}
/* ----------------------------------------------------------------- */
int isSubString(vector x, vector y, vector positions) {
  int i;

  for(i=0;i<positions[0];i++) {
    if (x[i]!=y[positions[i+1]-1]) {
      return(0);
    }
  }
  return(1);
}
/* ----------------------------------------------------------------- */
int genmodel_main(int argc, char *argv[]) {
  int i,j,maxIndexFace,numOfNodes,numOfRows,numOfColumns,infoLevel;
  vector face, faceValues, column, levels, nodes;
  listVector *faces, *tmp;
  char fileName[127],outFileName[127];
  FILE *out;

  infoLevel=standardInfoLevel;
  for (i=1; i<argc-1; i++) {
    if (strncmp(argv[i],"--",2)==0) {
      if (strncmp(argv[i],"--qui",5)==0) {
	infoLevel=-1;
      }
    }
  }

if (infoLevel>-1) {
  printVersionInfo();
}

  strcpy(fileName,argv[argc-1]);
  strcat(fileName,".mod");
  strcpy(outFileName,argv[argc-1]);

if (infoLevel>-1) {
  printf("Creating file %s.\n",outFileName);
}

  numOfNodes=0;
  faces=readSimplicialComplex(fileName,&numOfNodes);
  levels=faces->first;
  faces=faces->rest;

/*   strcpy(outFileName,fileName); */
/*   strcat(outFileName,".mod"); */

  if (!(out = fopen(outFileName,"w"))) {
    printf("Error opening file for output.");
    exit (0);
  }

  numOfColumns=1;
  for(i=0;i<numOfNodes;i++) numOfColumns*=levels[i];

  numOfRows=0;
  tmp=faces;
  while (tmp) {
    face=tmp->first;
    maxIndexFace=1;
    for(i=0;i<face[0];i++) maxIndexFace*=levels[face[i+1]-1];
    numOfRows+=maxIndexFace;
    tmp=tmp->rest;
  }

  fprintf(out,"%d %d\n",numOfRows,numOfColumns);

  nodes=createVector(numOfNodes+1);
  for (i=0;i<numOfNodes+1;i++) nodes[i]=i;
  nodes[0]=numOfNodes;

  while (faces) {
    face=faces->first;
    maxIndexFace=1;
    for(i=0;i<face[0];i++) maxIndexFace*=levels[face[i+1]-1];
    for(i=0;i<maxIndexFace;i++) {
      faceValues=decomposeIntegerIntoLevelIndices(i,face[0],face,levels);
      for(j=0;j<numOfColumns;j++) {
        column=decomposeIntegerIntoLevelIndices(j,numOfNodes,nodes,levels);
	fprintf(out,"%d ",isSubString(faceValues,column,face));
      }
    fprintf(out,"\n");
    }
    faces=faces->rest;
  }

  fclose(out);

  return(0);
}
/* ----------------------------------------------------------------- */
