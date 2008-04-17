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
listVector* extractStabilizer(listVector *basis, listVector *S, int numOfVars) {
  int check;
  listVector *tmp, *stab, *endStab;

  stab=createListVector(0);
  endStab=stab;

  while (S) {
    check=1;
    tmp=basis;
    while (tmp) {
      if (isVectorEqualToPermutedVector(tmp->first,S->first,numOfVars)==0) {
	tmp=0;
	check=0;
      }
      if (tmp) tmp=tmp->rest;
    }
    if (check==1) {
      endStab->rest=createListVector(S->first);
      endStab=endStab->rest;
    }
    S=S->rest;
  }

  return(stab->rest);
}
/* ----------------------------------------------------------------- */
int isVectorFixed(vector v, vector fix, int numOfFixedPoints) {
  int i;

  for (i=0;i<numOfFixedPoints;i++)
    if (v[fix[i]]!=fix[i]) return (0);

  return (1);
}
/* ----------------------------------------------------------------- */
int isVectorRelaxedFixed(vector v, vector fix, int numOfFixedPoints) {
  int i,j,ok;

  for (i=0;i<numOfFixedPoints;i++) {
    ok=0;
    for (j=0;j<numOfFixedPoints;j++) if (v[fix[i]]==fix[j]) ok=1;
    if (ok==0) return (0);
  }
  return (1);
}
/* ----------------------------------------------------------------- */
listVector* extractFixedVectors(listVector *basis, vector fix, 
				int numOfFixedPoints) {
  listVector *F, *endF;

  F=createListVector(0);
  endF=F;

  while (basis) {
    if (isVectorFixed(basis->first,fix,numOfFixedPoints)==1) {
      endF->rest=createListVector(basis->first);
      endF=endF->rest;
    }
    basis=basis->rest;
  }

  return(F->rest);
}
/* ----------------------------------------------------------------- */
listVector* extractRelaxedFixedVectors(listVector *basis, vector fix, 
				       int numOfFixedPoints) {
  listVector *F, *endF;

  F=createListVector(0);
  endF=F;

  while (basis) {
    if (isVectorRelaxedFixed(basis->first,fix,numOfFixedPoints)==1) {
      endF->rest=createListVector(basis->first);
      endF=endF->rest;
    }
    basis=basis->rest;
  }

  return(F->rest);
}
/* ----------------------------------------------------------------- */
int isVectorDominatedByVector(vector v, vector w, int numOfVars) {
  int i;

  for (i=0; i<numOfVars; i++) {
    if (v[i]>w[i]) return (0);     
  }
  return (1);
}
/* ----------------------------------------------------------------- */
int isVectorDominatedByListVector(vector v, listVector *dom, int numOfVars) {


  while (dom) {
    if (isVectorDominatedByVector(v,dom->first,numOfVars)==1) return (1);
    dom=dom->rest;
  }

  return (0);
}
/* ----------------------------------------------------------------- */
listVector* extractNonDominatedVectors(listVector *basis, listVector *dom, 
				       int numOfVars) {
  int lenBasis, count;
  listVector *F, *endF, *tmp;

  F=createListVector(0);
  endF=F;

  lenBasis=lengthListVector(basis);
  count=0;

  while (basis) {
    count ++;
    if (count==(100000*(count/100000))) 
      printf("Considering vector %d/%d\n",count,lenBasis);
    if (isVectorDominatedByListVector(basis->first,dom,numOfVars)==0) {
      endF->rest=createListVector(basis->first);
      endF=endF->rest;
    } else {
      free(basis->first);
    }
    tmp=basis;
    basis=basis->rest;
    free(tmp);
  }

  return(F->rest);
}
/* ----------------------------------------------------------------- */
listVector* extractMaximalNonDominatedVectors(listVector *basis, 
					      listVector *symmGroup,
					      int numOfVars) {
  int count, maxNorm;
  vector maxVector;
  listVector *tmp, *M, *endM, *orbit;

  M=createListVector(0);
  endM=M;
  maxVector=0;

  count=0;
  printf("%d nondominated vectors found, %d vectors left to consider\n",
	 count,lengthListVector(basis));

  while (basis) {
    maxNorm=maximalNormInListVector(basis,numOfVars);
    tmp=basis;
    while (tmp) {
      if (normOfVector(tmp->first,numOfVars)==maxNorm) {
	count++;
	maxVector=tmp->first;
	endM->rest=createListVector(maxVector);
	endM=endM->rest;
	tmp=0;
      } else {
	tmp=tmp->rest;
      }
    }

    orbit=expandRepresentativeIntoFullOrbits(createListVector(maxVector),
					     symmGroup,numOfVars,10); 

    basis=extractNonDominatedVectors(basis,orbit,numOfVars);
    printf("%d nondominated vectors found, %d vectors left to consider\n",
	   count,lengthListVector(basis));
  }

  return(M->rest);
}
/* ----------------------------------------------------------------- */
int output_main(int argc, char *argv[]) {
  int i,j,x,y,z,numOfVars,numOfRows,numOfFixPoints,numOfLabels,infoLevel,
    degree,lowdegree,highdegree,coord,sizeOfLayer,val;
  char *s;
  char fileName[127],outFileName[127],domFileName[127],symFileName[127],
    varFileName[127],groFileName[127],costFileName[127];
  char **labels;
  vector v,w,fixpoints;
  listVector *A, *B, *C, *basis, *domBasis, *tmp, *tmpV, *symmGroup, *weights;
  FILE *in;

  infoLevel=standardInfoLevel;
  for (i=1;i<argc-1;i++) {
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
	strcat(symFileName,".full");
	printListVectorToFile(symFileName, symmGroup, numOfVars);
      }

      basis=extractSymmetryRepresentatives(basis,symmGroup,numOfVars);
      strcpy(outFileName,fileName);
      strcat(outFileName,".rep");
      printListVectorToFile(outFileName,basis,numOfVars);
      printf("%d representatives found.\n",lengthListVector(basis));
    }
    if (strncmp(argv[i],"--dom",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(domFileName,argv[argc-2]);
      domBasis=readListVector(&numOfVars,domFileName);

      basis=extractNonDominatedVectors(basis,domBasis,numOfVars); 
      strcpy(outFileName,fileName);
      strcat(outFileName,".nondom");
      printListVectorToFile(outFileName,basis,numOfVars);
      printf("%d non-dominated vectors found.\n",
	     lengthListVector(basis));
    }
    if (strncmp(argv[i],"--max",5)==0) {
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
	strcat(symFileName,".full");
	printListVectorToFile(symFileName, symmGroup, numOfVars);
      }
      basis=extractMaximalNonDominatedVectors(basis,symmGroup,numOfVars); 
      strcpy(outFileName,fileName);
      strcat(outFileName,".maxnondom");
      printListVectorToFile(outFileName,basis,numOfVars);
      printf("%d maximal non-dominated vectors found.\n",
	     lengthListVector(basis));
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
	if (argc==4) {
	  degree=atoi(argv[i+1]);
	  printf("degree = %d\n",degree);
	  strcpy(outFileName,fileName);
	  strcat(outFileName,".deg.");
	  strcat(outFileName,argv[i+1]);
	  printListVectorWithGivenDegreesToFile(outFileName,basis,numOfVars,
						degree,degree);
	} else {
	  lowdegree=atoi(argv[i+1]);
	  highdegree=atoi(argv[i+2]);
	  printf("degrees %d..%d\n",lowdegree,highdegree);
	  strcpy(outFileName,fileName);
	  strcat(outFileName,".deg.");
	  strcat(outFileName,argv[i+1]);
	  strcat(outFileName,"-");
	  strcat(outFileName,argv[i+2]);
	  printListVectorWithGivenDegreesToFile(outFileName,basis,numOfVars,
						lowdegree,highdegree);
	}
      }
      return(0);
    }
    if (strncmp(argv[i],"--sup",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      if (argc==3) {
	printSupportsOfListVector(basis,numOfVars);
      } else {
	if (argc==4) {
	  degree=atoi(argv[i+1]);
	  printf("size of support = %d\n",degree);
	  strcpy(outFileName,fileName);
	  strcat(outFileName,".supp.");
	  strcat(outFileName,argv[i+1]);
	  printListVectorWithGivenSupportsToFile(outFileName,basis,numOfVars,
						 degree,degree);
	} else {
	  lowdegree=atoi(argv[i+1]);
	  highdegree=atoi(argv[i+2]);
	  printf("sizes of support %d..%d\n",lowdegree,highdegree);
	  strcpy(outFileName,fileName);
	  strcat(outFileName,".supp.");
	  strcat(outFileName,argv[i+1]);
	  strcat(outFileName,"-");
	  strcat(outFileName,argv[i+2]);
	  printListVectorWithGivenSupportsToFile(outFileName,basis,numOfVars,
						 lowdegree,highdegree);
	}
      }
      return(0);
    }
    if (strncmp(argv[i],"--typ",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      sizeOfLayer=atoi(argv[i+1]);
      printTypesOfListVector(basis,sizeOfLayer,numOfVars);
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
    if (strncmp(argv[i],"--AxB",5)==0) {
      strcpy(fileName,argv[argc-3]);
      A=readListVector(&numOfVars,fileName);
      strcpy(fileName,argv[argc-2]);
      B=readListVector(&numOfVars,fileName);
      v=matrixTimesVector(A,B->first,lengthListVector(A),numOfVars);
      C=createListVector(v);
      strcpy(fileName,argv[argc-1]);
      printListVectorToFile(fileName,C,lengthListVector(A));
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
    if (strncmp(argv[i],"--sum",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      v=createVector(numOfVars);
      for (i=0; i<numOfVars; i++) v[i]=0;
      while (basis) {
	for (i=0; i<numOfVars; i++) v[i]=v[i]+(basis->first)[i];
	basis=basis->rest;
      }
      printVector(v,numOfVars);
    }
    if (strncmp(argv[i],"--sub",5)==0) {
      strcpy(fileName,argv[argc-2]);
      basis=readListVector(&numOfVars,fileName);
      v=basis->first;
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(outFileName,fileName);
      strcat(outFileName,".submat");
      printSubsetOfListVectorToFile(outFileName,basis,v,numOfVars);
    }
    if (strncmp(argv[i],"--rem",5)==0) {
      if (strncmp(argv[i],"--remcol",8)==0) {
	/* Remove column. */
        strcpy(fileName,argv[argc-1]);
        basis=readListVector(&numOfVars,fileName);
	coord=atoi(argv[argc-2]);
        strcpy(outFileName,fileName);
        strcat(outFileName,".remcol");
        printListVectorWithoutColumnToFile(outFileName,basis,coord,numOfVars);
      }
    }
    if (strncmp(argv[i],"--sta",5)==0) {
      /* Extracting those symmetries from a list of given symmetries 
	 that keep a given list of vectors fixed. */
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(symFileName,argv[argc-2]);
      symmGroup=readListVector(&numOfVars,symFileName);
      symmGroup=extractStabilizer(basis,symmGroup,numOfVars);
      strcpy(outFileName,symFileName);
      strcat(outFileName,".stab");
      printListVectorToFile(outFileName,symmGroup,numOfVars);
    }
    if (strncmp(argv[i],"--fil",5)==0) {
      /* Fill specified column with specified value. */

      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      coord=atoi(argv[argc-3]);
      val=atoi(argv[argc-2]);
      tmp=basis;
      while (tmp) {
	(tmp->first)[coord-1]=val;
	tmp=tmp->rest;
      }
      strcpy(outFileName,fileName);
      strcat(outFileName,".fil");
      printListVectorToFile(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--add",5)==0) {
      /* Add column with specified value. */
        strcpy(fileName,argv[argc-1]);
        basis=readListVector(&numOfVars,fileName);
	coord=atoi(argv[argc-3]);
	val=atoi(argv[argc-2]);
        strcpy(outFileName,fileName);
        strcat(outFileName,".addcol");
        printListVectorWithAdditionalColumnToFile(outFileName,basis,
						 coord,val,numOfVars);
    }
    if (strncmp(argv[i],"--fix",5)==0) {
      /* Extracting those vectors that have given coordinates x[i]=i. */

      numOfFixPoints=argc-3;
      fixpoints=createVector(numOfFixPoints);
      for (j=2;j<argc-1;j++) fixpoints[j-2]=atoi(argv[j]);
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      basis=extractFixedVectors(basis,fixpoints,numOfFixPoints);
      strcpy(outFileName,fileName);
      strcat(outFileName,".fix");
      printListVectorToFile(outFileName,basis,numOfVars);
    }
    if (strncmp(argv[i],"--fox",5)==0) {
      /* Extracting those vectors that have given coordinates x[i]=i. */

      numOfFixPoints=argc-3;
      fixpoints=createVector(numOfFixPoints);
      for (j=2;j<argc-1;j++) fixpoints[j-2]=atoi(argv[j]);
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      basis=extractRelaxedFixedVectors(basis,fixpoints,numOfFixPoints);
      strcpy(outFileName,fileName);
      strcat(outFileName,".fox");
      printListVectorToFile(outFileName,basis,numOfVars);
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
