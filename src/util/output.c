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
/*------------------------------------------------------------------ */
#include "myheader.h"
#include "orbit.h"
#include "print.h"
#include "vector.h"
#include <stdlib.h>
/* ----------------------------------------------------------------- */
listVector* readListVector(int *numOfVars, char *fileName) {
  int numOfVectors;
  listVector *basis, *endBasis;
  vector b;
  FILE *in;

  setbuf(stdout,0);
  if (!(in = fopen(fileName,"r"))) {
    printf("File \"%s\" not found for reading!\n",fileName);
    return(0);
  }

  fscanf(in,"%d",&numOfVectors);
  fscanf(in,"%d",numOfVars);

  if (numOfVectors==0) return (0);

  b=createVector(*numOfVars);
  for (int j=0; j<(*numOfVars); j++) fscanf(in,"%d",&b[j]);
  basis = createListVector(b);
  endBasis = basis;

  for (int i=1; i<numOfVectors; i++) {
    b=createVector(*numOfVars);
    for (int j=0; j<(*numOfVars); j++) fscanf(in,"%d",&b[j]);
    endBasis = updateBasis(createListVector(b), endBasis);
  }
  fclose(in);
  return(basis);
}
/* ----------------------------------------------------------------- */
listVector* extractPositivePartsOfVectors(listVector *basis, 
					  int numOfVars) {
  vector v;
  listVector *tmp;

  tmp=basis;

  while (tmp) {
    v=tmp->first;
    for (int i=0;i<numOfVars;i++) if (v[i]<0) v[i]=0;
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
listVector* extractInitialForms(listVector *basis, vector w, int numOfVars) {
  int s;
  vector v;
  listVector *tmp;

  tmp=basis;

  while (tmp) {
    v=tmp->first;
    s=dotProduct(v,w,numOfVars);
    if (s>0)
      for (int i=0;i<numOfVars;i++) { if (v[i]<0) v[i]=0; }
    tmp->first=v;
    tmp=tmp->rest;
  }

  return(basis);
}
/* ----------------------------------------------------------------- */
int output_main(int argc, char *argv[]) {
  int x,y,z,numOfVars,numOfRows, numOfLabels,infoLevel,degree,coord;
  char *s;
  char fileName[127],outFileName[127],symFileName[127],varFileName[127],
    groFileName[127],costFileName[127];
  char **labels;
  vector v,w;
  listVector *basis, *tmpV, *symmGroup, *weights;
  FILE *in;

  infoLevel=standardInfoLevel;
  for (int i=1; i<argc-1; i++) {
    if (strncmp(argv[i],"--",2)==0) {
      if (strncmp(argv[i],"--qui",5)==0) {
	infoLevel=-1;
      }
    }
  }

if (infoLevel>-1) {
  printVersionInfo();
}

  for (int i=1; i<argc; i++) {
    if (strncmp(argv[i],"--pos",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      basis=extractPositivePartsOfVectors(basis,numOfVars);
      strcpy(outFileName,fileName);
      strcat(outFileName,".pos");
      printListVectorToFile(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--rep",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(symFileName,fileName);
      strcat(symFileName,".sym.full");
      symmGroup=readListVector(&numOfVars,symFileName);

      if (symmGroup==0) { 
        strcpy(fileName,argv[argc-1]);
        basis=readListVector(&numOfVars,fileName);
	strcpy(symFileName,fileName);
	strcat(symFileName,".sym");
	symmGroup=readListVector(&numOfVars,symFileName);

	tmpV=symmGroup;
	while (tmpV) {
	  v=tmpV->first;
	  for (int i=0; i<numOfVars; i++) v[i]=v[i]-1;
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
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
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
	  for (int i=0; i<numOfVars; i++) v[i]=v[i]-1;
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
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
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
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
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
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      basis=extractZeroOneVectors(basis,numOfVars);
      strcpy(outFileName,fileName);
      strcat(outFileName,".0-1");
      printListVectorToFile(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--3wa",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      x=atoi(argv[i+1]);
      y=atoi(argv[i+2]);
      z=atoi(argv[i+3]);
      strcpy(outFileName,fileName);
      strcat(outFileName,".3way");
      print3wayTables(outFileName,basis,x,y,z,numOfVars);
    }
    if (strncmp(argv[i],"--tra",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(outFileName,fileName);
      strcat(outFileName,".tra");
      printTransposedListVectorToFile(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--map",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(outFileName,fileName);
      strcat(outFileName,".maple");
      printListVectorMaple(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--mac",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(outFileName,fileName);
      strcat(outFileName,".macaulay2");
      printListVectorMacaulay2(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--mat",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(outFileName,fileName);
      strcat(outFileName,".mathematica");
      printListVectorMaple(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--coc",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(outFileName,fileName);
      strcat(outFileName,".cocoa");
      printListVectorMaple(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--bin",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(outFileName,fileName);
      strcat(outFileName,".bin");

      labels=0;
      strcpy(varFileName,fileName);
      strcat(varFileName,".vars");
      if ((in = fopen(varFileName,"r"))) {
	printf("File \"%s\" found. 4ti2 will use it.\n\n",varFileName);
	if (fscanf(in,"%d %d",&numOfRows, &numOfLabels)!=2 || numOfRows!=1) {
          printf("ERROR: Unrecognised file format for \"%s\".\n", varFileName);
          exit(1);
        }
	if (numOfLabels != numOfVars) {
	  printf("ERROR: Incorrect number of variable names in \"%s\".\n",
                          varFileName);
          exit(1);
	}
	labels = (char **)malloc(sizeof(char*)*(numOfVars));
	for (int i=0; i<numOfVars; i++) {
	  s=(char *)malloc(sizeof(char)*127);
	  if (fscanf(in,"%s",s) != 1) {
            printf("ERROR: Unrecognised file format for \"%s\".\n",
                            varFileName);
            exit(1);
          }
	  labels[i]=s;
	}
	fclose(in);
      }
      printListBinomialsToFile(outFileName,basis,numOfVars,labels);
    }
    if (strncmp(argv[i],"--ini",5)==0) {
      strcpy(fileName,argv[argc-1]);
      strcpy(groFileName,fileName);
      strcat(groFileName,".gro");
      basis=readListVector(&numOfVars,groFileName);
      strcpy(costFileName,argv[argc-1]);
      strcat(costFileName,".cost");
      weights=readListVector(&numOfVars,costFileName);
      if (weights!=0) {
        w=weights->first;
      } else {
	w=createVector(numOfVars);
	for (int i=0;i<numOfVars;i++) w[i]=1;
      }
      basis=extractInitialForms(basis,w,numOfVars);

      strcpy(outFileName,fileName);
      strcat(outFileName,".ini");
      printListVectorToFile(outFileName,basis,numOfVars);

      labels=0;
      strcpy(varFileName,fileName);
      strcat(varFileName,".vars");
      if ((in = fopen(varFileName,"r"))) {
	printf("File \"%s\" found. 4ti2 will use it.\n\n",varFileName);
	fscanf(in,"%d %d",&numOfRows, &numOfLabels);
	labels = (char **)malloc(sizeof(char*)*(numOfVars));
	if (fscanf(in,"%d %d",&numOfRows, &numOfLabels)!=2 || numOfRows!=1) {
          printf("ERROR: Unrecognised file format for \"%s\".\n", varFileName);
          exit(1);
        }
	if (numOfLabels != numOfVars) {
	  printf("ERROR: Incorrect number of variable names in \"%s\".\n",
                          varFileName);
          exit(1);
	}
	for (int i=0; i<numOfVars; i++) {
	  s=(char *)malloc(sizeof(char)*127);
	  if (fscanf(in,"%s",s) != 1) {
            printf("ERROR: Unrecognised file format for \"%s\".\n",
                            varFileName);
            exit(1);
          }
	  labels[i]=s;
	}
	fclose(in);
      }
      strcpy(outFileName,fileName);
      strcat(outFileName,".ini.bin");
      printListMonomialsAndBinomialsToFile(outFileName,basis,numOfVars,labels);
    }
  }

  return (0);
}
/* ----------------------------------------------------------------- */
