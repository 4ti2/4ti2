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
/* Generation of toric models in statistics                          */
/*                                                                   */
/* Author   : Raymond Hemmecke                                       */
/*                                                                   */
/* Created    : 01-JUL-03                                            */
/*                                                                   */
/* ----------------------------------------------------------------- */
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include "myheader.h"
#include "vector.h"
#include "print.h"

#include <getopt.h>
#include "banner.h"

/* ----------------------------------------------------------------- */
listVector* readSimplicialComplex(char* inpFileName, int* numOfNodes) {
  int i,j,numOfFaces,sizeOfFace;
  listVector *basis, *endBasis;
  vector b;
  FILE *in;

  setbuf(stdout,0);
  if (!(in = fopen(inpFileName,"r"))) {
    printf("Error opening file %s containing the simplicial complex.\n",
	   inpFileName);
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
static const struct option longopts[] = {
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'v'},
  {"quiet", no_argument, NULL, 'q'},
  {NULL, 0, NULL, 0}
};

static void print_version()
{
  printf("%s", FORTY_TWO_BANNER);
}
 
static void print_usage()
{
  printf("usage: genmodel [--options] FILENAME\n"
	 "\n"
	 "Computes the problem matrix corresponding to graphical statistical models\n"
	 "given by a simplicial complex and levels on the nodes.\n"
	 "\n"
	 "Options:\n"
	 " -q, --quiet	No output is written to the screen\n"
	 "\n"
	 "Input file:\n"
	 "FILENAME.mod    Simplicial complex and levels on the nodes\n"
	 "\n"
	 "Output file:\n"
	 "FILENAME.mat    Matrix file\n"
	 "\n"
	 "Example: Consider the problem of 3x3x3 tables with 2-marginals. These\n"
	 "are given by K_3 as the simplicial complex on 3 nodes and with levels\n"
	 "of 3 on each node.  In '333.mod' write:\n"
	 "3\n"
	 "3 3 3\n"
	 "3\n"
	 "2 1 2\n"
	 "2 2 3\n"
	 "2 3 1\n"
	 "Calling 'genmodel 333' produces the following file '333.mat':\n"
	 "27 27\n"
	 "1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0\n"
	 "0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0\n"
	 "[...]\n"
	 /* "0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0\n" */
	 /* "0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0\n" */
	 /* "0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0\n" */
	 /* "0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0\n" */
	 /* "0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0\n" */
	 /* "0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0\n" */
	 /* "0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1\n" */
	 "1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n"
	 "0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n"
	 "0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n"
	 "[...]\n"
	 /* "0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0\n" */
	 /* "0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0\n" */
	 /* "0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0\n" */
	 /* "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0\n" */
	 /* "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0\n" */
	 /* "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1\n" */
	 "1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n"
	 "0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n"
	 "[...]\n"
	 /* "0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n" */
	 /* "0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n" */
	 /* "0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0\n" */
	 /* "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0\n" */
	 /* "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0\n" */
	 /* "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0\n" */
	 /* "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1\n" */
	 "\n");
}

/* ----------------------------------------------------------------- */
int genmodel_main(int argc, char *argv[]) {
  int i,j,maxIndexFace,numOfNodes,numOfRows,numOfColumns,infoLevel;
  vector face, faceValues, column, levels, nodes;
  listVector *faces, *tmp;
  const char *fileName=NULL;
  char *inpFileName=NULL, *outFileName=NULL;
  FILE *out;

  int optc;

  infoLevel=standardInfoLevel;


  while ((optc = getopt_long (argc, argv, "hvq", longopts, NULL)) != -1)
    switch (optc) {
    case 'v':
      print_version();
      exit(0);
      break;
    case 'h':
      print_usage();
      exit(0);
      break;
    case 'q':
      infoLevel=-1;
      break;
    default:
      print_usage();
      exit(1);
      break;
    }

  if (optind != argc - 1) {
    print_usage();
    exit(1);
  }

  if (infoLevel>-1) {
    printVersionInfo();
  }

  fileName=argv[argc-1];
  myxasprintf(&inpFileName,"%s.mod",fileName);
  myxasprintf(&outFileName,"%s.mat",fileName);

  if (infoLevel>-1) {
    printf("Creating file %s.\n",outFileName);
  }

  numOfNodes=0;
  faces=readSimplicialComplex(inpFileName,&numOfNodes);
  levels=faces->first;
  faces=faces->rest;

  errno=0;
  if (!(out = fopen(outFileName,"w"))) {
    printf("Error opening file for output (%s).",strerror(errno));
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

  free(outFileName);
  free(inpFileName);

  return(0);
}
/* ----------------------------------------------------------------- */
