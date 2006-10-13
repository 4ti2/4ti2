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
/* Transformation of the standard output to something else.          */
/*                                                                   */
/* Author   : Raymond Hemmecke   (implementation)                    */
/* Co-Author: Ralf Hemmecke      (data structure, code optimization) */
/*                                                                   */
/* $Id$             */
/*                                                                   */
/*------------------------------------------------------------------ */
#include "myheader.h"
#include "orbit.h"
#include "print.h"
#include "vector.h"
#include <stdlib.h>
/* ----------------------------------------------------------------- */
listVector* readListVector(int *numOfVars, char *fileName) {
  int i,j,numOfVectors;
  listVector *basis, *endBasis;
  vector b;
  FILE *in;

  setbuf(stdout,0);
  if (!(in = fopen(fileName,"r"))) {
    printf("Error opening file \"%s\" containing list of vectors!",fileName);
    (*numOfVars)=0;
    return(0);
  }

  fscanf(in,"%d",&numOfVectors);
  fscanf(in,"%d",numOfVars);

  if (numOfVectors==0) return (0);

  b=createVector(*numOfVars);
  for (j=0; j<(*numOfVars); j++) fscanf(in,"%d",&b[j]);
  basis = createListVector(b);
  endBasis = basis;

  for (i=1; i<numOfVectors; i++) {
    b=createVector(*numOfVars);
    for (j=0; j<(*numOfVars); j++) fscanf(in,"%d",&b[j]);
    endBasis = updateBasis(createListVector(b), endBasis);
  }
  fclose(in);
  return(basis);
}
/* ----------------------------------------------------------------- */
listVector* extractPositivePartsOfVectors(listVector *basis, 
					  int numOfVars) {
  int i;
  vector v;
  listVector *tmp;

  tmp=basis;

  while (tmp) {
    v=tmp->first;
    for (i=0;i<numOfVars;i++) if (v[i]<0) v[i]=0;
    tmp->first=v;
    tmp=tmp->rest;
  }

  return(basis);
}
/* ----------------------------------------------------------------- */
listVector* extractVectorsWithFirstEntryEqualToOne(listVector *basis, 
						   int numOfVars) {
  vector v;
  listVector *tmp, *endBasis;

  tmp=basis;
  basis=0;
  endBasis=0;

  while (tmp) {
    v=tmp->first;
    if (abs(v[0])==1) {
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
int output_main(int argc, char *argv[]) {
  int i,x,y,z,numOfVars,numOfLabels,infoLevel,degree,coord;
  char *s;
  char fileName[127],outFileName[127],symFileName[127],varFileName[127];
  char **labels;
  vector v;
  listVector *basis, *circuits, *rays, *tmpV, *symmGroup;
  FILE *in;

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
  basis=readListVector(&numOfVars,fileName);

  for (i=1; i<argc; i++) {
    if (strncmp(argv[i],"--pos",5)==0) {
      basis=extractPositivePartsOfVectors(basis,numOfVars);
      strcpy(outFileName,fileName);
      strcat(outFileName,".pos");
      printListVectorToFile(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--sol",5)==0) {
      basis=extractVectorsWithFirstEntryEqualToOne(basis,numOfVars);
      strcpy(outFileName,fileName);
      strcat(outFileName,".sol");
      printListVectorToFile(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--rep",5)==0) {
      strcpy(symFileName,fileName);
      strcat(symFileName,".sym.full");
      symmGroup=readListVector(&numOfVars,symFileName);

      if (symmGroup==0) { 
	strcpy(symFileName,fileName);
	strcat(symFileName,".sym");
	symmGroup=readListVector(&numOfVars,symFileName);

	tmpV=symmGroup;
	while (tmpV) {
	  v=tmpV->first;
	  for (i=0; i<numOfVars; i++) v[i]=v[i]-1;
	  tmpV->first=v;
	  tmpV=tmpV->rest;
	}
	symmGroup=generateSymmetryGroup(symmGroup,numOfVars);
      }

      basis=extractSymmetryRepresentatives(basis,symmGroup,numOfVars);
      strcpy(outFileName,fileName);
      strcat(outFileName,".rep");
      printListVectorToFile(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--exp",5)==0) {
      strcpy(symFileName,fileName);
      strcat(symFileName,".sym.full");
      symmGroup=readListVector(&numOfVars,symFileName);

      if (symmGroup==0) { 
	strcpy(symFileName,fileName);
	strcat(symFileName,".sym");
	symmGroup=readListVector(&numOfVars,symFileName);
	tmpV=symmGroup;
	while (tmpV) {
	  v=tmpV->first;
	  for (i=0; i<numOfVars; i++) v[i]=v[i]-1;
	  tmpV->first=v;
	  tmpV=tmpV->rest;
	}
	symmGroup=generateSymmetryGroup(symmGroup,numOfVars);
      }
      basis=expandRepresentativeIntoFullOrbits(basis,symmGroup,numOfVars,10);
      strcpy(outFileName,fileName);
      strcat(outFileName,".exp");
      printListVectorToFile(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--deg",5)==0) {
      if (argc==3) {
	printL1NormOfListVector(basis,numOfVars);
      } else {
	degree=atoi(argv[i+1]);
	printf("degree = %d\n",degree);
	strcpy(outFileName,fileName);
	strcat(outFileName,".deg.");
	strcat(outFileName,argv[i+1]);
	printListVectorWithGivenDegreeToFile(outFileName,basis,numOfVars,
					     degree);
      }
      return(0);
    }
    if (strncmp(argv[i],"--non",5)==0) {
      if (argc==3) {
	printf("You need to specify a coordinate!\n");
	return(0);
      } else {
	coord=atoi(argv[i+1]);
	strcpy(outFileName,fileName);
	strcat(outFileName,".nonzero.");
	strcat(outFileName,argv[i+1]);
	printListVectorWithGivenNonzeroEntryToFile(outFileName,basis,
						   numOfVars,coord);
      }
      return(0);
    }
    if (strncmp(argv[i],"--0-1",5)==0) {
      basis=extractZeroOneVectors(basis,numOfVars);
      strcpy(outFileName,fileName);
      strcat(outFileName,".0-1");
      printListVectorToFile(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--3wa",5)==0) {
      x=atoi(argv[i+1]);
      y=atoi(argv[i+2]);
      z=atoi(argv[i+3]);
      strcpy(outFileName,fileName);
      strcat(outFileName,".3way");
      print3wayTables(outFileName,basis,x,y,z,numOfVars);
    }
    if (strncmp(argv[i],"--tra",5)==0) {
      strcpy(outFileName,fileName);
      strcat(outFileName,".tra");
      printTransposedListVectorToFile(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--map",5)==0) {
      strcpy(outFileName,fileName);
      strcat(outFileName,".maple");
      printListVectorMaple(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--mac",5)==0) {
      strcpy(outFileName,fileName);
      strcat(outFileName,".macaulay2");
      printListVectorMacaulay2(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--mat",5)==0) {
      strcpy(outFileName,fileName);
      strcat(outFileName,".mathematica");
      printListVectorMaple(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--coc",5)==0) {
      strcpy(outFileName,fileName);
      strcat(outFileName,".cocoa");
      printListVectorMaple(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--bin",5)==0) {
      strcpy(outFileName,fileName);
      strcat(outFileName,".bin");

      labels=0;
      strcpy(varFileName,fileName);
      strcat(varFileName,".vars");
      if ((in = fopen(varFileName,"r"))) {
	printf("File \"%s\" found. 4ti2 will use it.\n\n",varFileName);
	fscanf(in,"%d",&numOfLabels);
	if (numOfLabels<numOfVars) {
	  printf("There are not enough names in \"%s\" to rename variables.\n",
		 varFileName);
	  printf("4ti2 will use the standard names x[1], ..., x[%d].\n",
		 numOfVars);
	} else {
	  if (numOfLabels>numOfVars) {
	    printf("There too many names in \"%s\" to rename variables.\n",
		   varFileName);
	    printf("4ti2 will use the first %d names for renaming.\n",
		   numOfVars);
	  }
	  labels = (char **)malloc(sizeof(char*)*(numOfVars));
	  for (i=0; i<numOfVars; i++) {
	    s=(char *)malloc(sizeof(char)*127);
	    fscanf(in,"%s",s);
	    labels[i]=s;
	  }
	}
	fclose(in);
      }
      printListBinomialsToFile(outFileName,basis,numOfVars,labels);
    }
    if (strncmp(argv[i],"--cir",5)==0) {
      strcpy(outFileName,fileName);
      strcat(outFileName,".cir");
      circuits=extractCircuits(basis,numOfVars);
      printListVectorToFile(outFileName,circuits,numOfVars);
    }
    if (strncmp(argv[i],"--ray",5)==0) {
      strcpy(outFileName,fileName);
      strcat(outFileName,".ray");
      rays=extractCircuits(basis,numOfVars);
      printListVectorToFile(outFileName,rays,numOfVars);
    }
  }

  return (0);
}
/* ----------------------------------------------------------------- */
