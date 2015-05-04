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
#include <limits.h>

#include "banner.h"

/* ----------------------------------------------------------------- */
listVector* readListVector(int *numOfVars, char *fileName) {
  int numOfVectors;
  listVector *basis, *endBasis;
  vector b;
  FILE *in;
  int i, j;

  setbuf(stdout,0);
  if (!(in = fopen(fileName,"r"))) {
    printf("File \"%s\" not found for reading!\n",fileName);
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
  vector v;
  listVector *tmp;

  tmp=basis;

  while (tmp) {
    int i;
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
listVector* extractInitialForms(listVector *basis, vector w, int numOfVars) {
  int s;
  vector v;
  listVector *tmp;

  tmp=basis;

  while (tmp) {
    int i;
    v=tmp->first;
    s=dotProduct(v,w,numOfVars);
    if (s>0)
      for (i=0;i<numOfVars;i++) { if (v[i]<0) v[i]=0; }
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
static void print_version()
{
  printf("%s", FORTY_TWO_BANNER);
}
 
static void print_usage()
{
  printf("usage: output [--options] FILENAME.EXT\n"
         "\n"
         "Transforms a 4ti2 matrix file to something else.\n"
         "\n"
         "General options:       \n"
         " --quiet        No output is written to the screen.\n"
         "\n"
         "Options that control what to output:                 their output files:\n"
         " --binomials    Write vectors as binomials.          FILENAME.EXT.bin\n"
         "                Use an optional input file\n"
         "                'FILENAME.EXT.vars'\n"
         "                to specify variable names.\n"
         "\n"
         " --maple        Write vectors as Maple list.         FILENAME.EXT.maple\n"
         "                This format is suitable also for\n"
         "                CoCoA, Mathematica, Macaulay2.\n"
         "\n"
         " --0-1          Extract vectors with 0-1\n"
         "                components only.                     FILENAME.EXT.0-1\n"
         "\n"
         " --transpose    Transpose matrix and write it        FILENAME.EXT.tra\n"
         "                in 4ti2 format.\n"
         "\n"
         " --degree       Print 1-norms of all vectors.\n"
         "\n"
         " --degree N     Extract all vectors of 1-norm        FILENAME.EXT.deg.N\n"
         "                equal to N.\n"
         "\n"
         " --degree N1 N2 Extract all vectors of 1-norm        FILENAME.EXT.deg.N1-N2\n"
         "                between N1 and N2 (inclusive).\n"
         "\n"
         " --support      Print supports of all vectors.\n"
         "\n"
         " --support S    Extract all vectors of support       FILENAME.EXT.supp.S\n"
	 "                size equal to S.\n"
         "\n"
         " --support S1 S2    Extract all vectors of support   FILENAME.EXT.supp.S1-S2\n"
	 "                    size between S1 and S2 (incl.)\n"
         "\n"
         " --positive     Extract positive parts of vectors.   FILENAME.EXT.pos\n"
         "                Corresponds to leading terms of\n"
         "                binomials.\n"
         "\n"
         " --3way A B C   Write vectors as 3-way tables        FILENAME.EXT.3way\n"
         "                of size A x B x C.\n"
         "\n"
         " --nonzero-at K Extract all vectors that have        FILENAME.EXT.nonzero.K\n"
	 "                nonzero K-th coordinate.\n"
         "\n"
         "Undocumented or obscure options for experts:\n"
         " --representatives\n"
         " --dominated    Extract all non-dominated vectors    FILENAME.EXT.nondom\n"
         " --maximal-non-dominated                             FILENAME.EXT.maxnondom\n"
         " --expand-representatives-to-full-orbits\n"
         " --type T\n"
         " --AxB          Computes a matrix-vector product.\n"
         " --macaulay2\n"
         " --mathematica\n"
         " --cocoa\n"
         " --sum          Print the sum of the columns.\n"
         " --submatrix LISTFILENAME                            FILENAME.EXT.submat\n"
         " --remove-column I                                   FILENAME.EXT.remcol\n"
         " --remcol I                                          FILENAME.EXT.remcol\n"
         " --stabilizer SYMMFILENAME                           FILENAME.EXT.stab\n"
         " --fill-column                                       FILENAME.EXT.fil\n"
         " --add-column                                        FILENAME.EXT.addcol\n"
         " --fix I1 ... IK    Extract fixed vectors,           FILENAME.EXT.fix\n"
	 "                    that is, those vectors that\n"
	 "                    have x[i]=i for the given i.\n"
         " --fox I1 ... IK    Extract relaxed fixed vectors.   FILENAME.EXT.fox\n"
         " --initial-forms    Extract initial forms.           FILENAME.ini\n"
	 "                    (Call with FILENAME rather       FILENAME.ini.bin\n"
         "                    than FILENAME.EXT.  Reads\n"
         "                    FILENAME.gro and\n"
         "                    optionally FILENAME.cost and\n"
	 "                    FILENAME.vars.\n"
         "\n"
         "Examples:\n"
         " 'output --binomials file.gra' writes the Graver basis elements as\n"
         " binomials in 'file.gra.bin'.\n"
         "\n"
         " 'output --0-1 foo.gra' extracts the 0-1 elements from the Graver basis\n"
	 " elements and writes them into 'foo.gra.0-1'.\n"
         "\n");
}
  
/* ----------------------------------------------------------------- */
int output_main(int argc, char *argv[]) {
  int i,j,x,y,z,numOfVars,numOfRows,numOfFixPoints,numOfLabels,infoLevel,
    degree,lowdegree,highdegree,coord,sizeOfLayer,val;
  char *s;
  char fileName[PATH_MAX],outFileName[PATH_MAX],domFileName[PATH_MAX],symFileName[PATH_MAX],
    varFileName[PATH_MAX],groFileName[PATH_MAX],costFileName[PATH_MAX];
  char **labels;
  vector v,w,fixpoints;
  listVector *A, *B, *C, *basis, *domBasis, *tmp, *tmpV, *symmGroup, *weights;
  FILE *in;
  int did_something = 0;
  
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

  for (i=1; i<argc; i++) {
    if (strncmp(argv[i],"--version",9)==0) {
      print_version();
      exit(0);
    }
    else if (strncmp(argv[i],"--help", 6)==0) {
      print_usage();
      exit(0);
    }
    else if (strncmp(argv[i],"--pos",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      basis=extractPositivePartsOfVectors(basis,numOfVars);
      strcpy(outFileName,fileName);
      strcat(outFileName,".pos");
      printListVectorToFile(outFileName,basis,numOfVars);
    }
    else if (strncmp(argv[i],"--rep",5)==0) {
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
	  int i;
	  v=tmpV->first;
	  for (i=0; i<numOfVars; i++) v[i]=v[i]-1;
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
    else if (strncmp(argv[i],"--dom",5)==0) {
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
    else if (strncmp(argv[i],"--max",5)==0) {
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
	  int i;
	  v=tmpV->first;
	  for (i=0; i<numOfVars; i++) v[i]=v[i]-1;
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
    else if (strncmp(argv[i],"--exp",5)==0) {
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
	  int i;
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
    else if (strncmp(argv[i],"--deg",5)==0) {
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
    else if (strncmp(argv[i],"--sup",5)==0) {
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
    else if (strncmp(argv[i],"--typ",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      sizeOfLayer=atoi(argv[i+1]);
      printTypesOfListVector(basis,sizeOfLayer,numOfVars);
      return(0);
    }
    else if (strncmp(argv[i],"--non",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      if (argc==3) {
	printf("ERROR: You need to specify a coordinate!\n");
	exit(1);
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
    else if (strncmp(argv[i],"--AxB",5)==0) {
      strcpy(fileName,argv[argc-3]);
      A=readListVector(&numOfVars,fileName);
      strcpy(fileName,argv[argc-2]);
      B=readListVector(&numOfVars,fileName);
      v=matrixTimesVector(A,B->first,lengthListVector(A),numOfVars);
      C=createListVector(v);
      strcpy(fileName,argv[argc-1]);
      printListVectorToFile(fileName,C,lengthListVector(A));
    }
    else if (strncmp(argv[i],"--0-1",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      basis=extractZeroOneVectors(basis,numOfVars);
      strcpy(outFileName,fileName);
      strcat(outFileName,".0-1");
      printListVectorToFile(outFileName,basis,numOfVars);
    }
    else if (strncmp(argv[i],"--3wa",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      x=atoi(argv[i+1]);
      y=atoi(argv[i+2]);
      z=atoi(argv[i+3]);
      strcpy(outFileName,fileName);
      strcat(outFileName,".3way");
      print3wayTables(outFileName,basis,x,y,z,numOfVars);
    }
    else if (strncmp(argv[i],"--tra",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(outFileName,fileName);
      strcat(outFileName,".tra");
      printTransposedListVectorToFile(outFileName,basis,numOfVars);
    }
    else if (strncmp(argv[i],"--map",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(outFileName,fileName);
      strcat(outFileName,".maple");
      printListVectorMaple(outFileName,basis,numOfVars);
    }
    else if (strncmp(argv[i],"--mac",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(outFileName,fileName);
      strcat(outFileName,".macaulay2");
      printListVectorMacaulay2(outFileName,basis,numOfVars);
    }
    else if (strncmp(argv[i],"--mat",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(outFileName,fileName);
      strcat(outFileName,".mathematica");
      printListVectorMaple(outFileName,basis,numOfVars);
    }
    else if (strncmp(argv[i],"--coc",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(outFileName,fileName);
      strcat(outFileName,".cocoa");
      printListVectorMaple(outFileName,basis,numOfVars);
    }
    else if (strncmp(argv[i],"--bin",5)==0) {
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(outFileName,fileName);
      strcat(outFileName,".bin");

      labels=0;
      strcpy(varFileName,fileName);
      strcat(varFileName,".vars");
      if ((in = fopen(varFileName,"r"))) {
	int i;
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
	for (i=0; i<numOfVars; i++) {
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
    else if (strncmp(argv[i],"--sum",5)==0) {
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
    else if (strncmp(argv[i],"--sub",5)==0) {
      strcpy(fileName,argv[argc-2]);
      basis=readListVector(&numOfVars,fileName);
      v=basis->first;
      strcpy(fileName,argv[argc-1]);
      basis=readListVector(&numOfVars,fileName);
      strcpy(outFileName,fileName);
      strcat(outFileName,".submat");
      printSubsetOfListVectorToFile(outFileName,basis,v,numOfVars);
    }
    else if (strncmp(argv[i],"--rem",5)==0) {
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
    else if (strncmp(argv[i],"--sta",5)==0) {
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
    else if (strncmp(argv[i],"--fil",5)==0) {
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
    else if (strncmp(argv[i],"--add",5)==0) {
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
    else if (strncmp(argv[i],"--fix",5)==0) {
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
    else if (strncmp(argv[i],"--fox",5)==0) {
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
    else if (strncmp(argv[i],"--ini",5)==0) {
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
	int i;
	w=createVector(numOfVars);
	for (i=0;i<numOfVars;i++) w[i]=1;
      }
      basis=extractInitialForms(basis,w,numOfVars);

      strcpy(outFileName,fileName);
      strcat(outFileName,".ini");
      printListVectorToFile(outFileName,basis,numOfVars);

      labels=0;
      strcpy(varFileName,fileName);
      strcat(varFileName,".vars");
      if ((in = fopen(varFileName,"r"))) {
	int i;
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
	for (i=0; i<numOfVars; i++) {
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
    else if (strncmp(argv[i],"--",2)==0) {
      printf("ERROR: Unknown option: %s\n\n", argv[i]);
      print_usage();
      exit(1);
    }
    else {
      /* All standard options take the FILENAME.EXT argument that
	 appears last.  Others that take several arguments exit
	 before we get here. So signal an error if there are too many
	 arguments left. */
      if (i == argc-1) {
	if (!did_something) {
	  printf("ERROR: At least one option that controls what to output needs to be given.\n\n");
	  print_usage();
	  exit(1);
	}
      }
      else {
	printf("ERROR: Unhandled command-line argument: %s\n\n", argv[i]);
	print_usage();
	exit(1);
      }
    }
    did_something = 1;
  }

  return (0);
}
/* ----------------------------------------------------------------- */
