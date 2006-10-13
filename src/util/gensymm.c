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
/* Generation of symmetry groups of 2- and 3-way tables              */
/*                                                                   */
/* Author   : Raymond Hemmecke                                       */
/*                                                                   */
/* Created    : 26-AUG-02                                            */
/* $Id$             */
/*                                                                   */
/* ----------------------------------------------------------------- */
#include <stdio.h>
#include <string.h>
#include "myheader.h"
#include "print.h"
#include "vector.h"
/* ----------------------------------------------------------------- */
void printPermutationToFile(FILE *out, vector v, int numOfVars) {
  int i;

  if (v==0) return;

  for (i=0; i<numOfVars; i++) {
    fprintf(out,"%d ", v[i]);
  }
  fprintf(out,"\n");
  return ;
}
/* ----------------------------------------------------------------- */
int gensymm_main(int argc, char *argv[]) {
  int a,i,j,k,l,x,y,z,w,numOfVars,numOfGenerators,infoLevel;
  vector v;
  char fileName[127],outFileName[127];
  FILE *out;

  setbuf(stdout,0);

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

/* We require that all 1's come last in the tuple (x,y,z,w). */

  strcpy(fileName,argv[argc-1]);
  x=atoi(argv[argc-5]);
  y=atoi(argv[argc-4]);
  z=atoi(argv[argc-3]);
  w=atoi(argv[argc-2]);

/* Sigh. I did not code the stuff below correctly. Damn. To fix it,
   I need to bring the numbers >1 in (x,y,z,w) into opposite
   order. */

  if (w>1) {
/* Switch order of x,y,z,w. */
    a=w; w=x; x=a;
    a=z; z=y; y=a; 
  } else {
    if (z>1) {
/* Switch order of x,y,z. */
      a=z; z=x; x=a;
    } else {
/* Switch order of x,y. */
      a=y; y=x; x=a;
    }
  }

/* Now we are back in business. */

  strcpy(outFileName,fileName);
  strcat(outFileName,".sym");

  if (!(out = fopen(outFileName,"w"))) {
    printf("Error opening generator file for output.");
    exit (0);
  }

/* Write dimensions. */

  numOfGenerators=4;
  if (z>1) numOfGenerators=numOfGenerators+2;
  if (w>1) numOfGenerators=numOfGenerators+2;
  if (x==y) numOfGenerators++;
  if (x==z) numOfGenerators++;
  if (x==w) numOfGenerators++;
  if (y==z) numOfGenerators++;
  if (y==w) numOfGenerators++;
  if ((z>1) && (z==w)) numOfGenerators++;

  numOfVars = x*y*z*w;
  fprintf(out,"%d %d\n",numOfGenerators,numOfVars);

  v = (vector)malloc(sizeof(int)*(numOfVars));

/* Permuting first component. */
  for (i=0; i<x; i++) {
    for (j=0; j<y; j++) {
      for (k=0; k<z; k++) {
	for (l=0; l<w; l++) {
	  if (i==x-1) {
	    v[i*(y*z*w)+j*(z*w)+k*w+l]=0*(y*z*w)+j*(z*w)+k*w+l+1;
	  } else{
	    v[i*(y*z*w)+j*(z*w)+k*w+l]=(i+1)*(y*z*w)+j*(z*w)+k*w+l+1;
	  }
	}
      }
    }
  }
  printPermutationToFile(out, v, numOfVars);

  for (i=0; i<x; i++) {
    for (j=0; j<y; j++) {
      for (k=0; k<z; k++) {
	for (l=0; l<w; l++) {
	  if (i==0) v[i*(y*z*w)+j*(z*w)+k*w+l]=1*(y*z*w)+j*(z*w)+k*w+l+1;
	  if (i==1) v[i*(y*z*w)+j*(z*w)+k*w+l]=0*(y*z*w)+j*(z*w)+k*w+l+1;
	  if (i>1)  v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+j*(z*w)+k*w+l+1;
	}
      }
    }
  }
  printPermutationToFile(out, v, numOfVars);

/* Permuting second component. */
  for (i=0; i<x; i++) {
    for (j=0; j<y; j++) {
      for (k=0; k<z; k++) {
	for (l=0; l<w; l++) {
	  if (j==y-1) {
	    v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+0*(z*w)+k*w+l+1;
	  } else{
	    v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+(j+1)*(z*w)+k*w+l+1;
	  }
	}
      }
    }
  }
  printPermutationToFile(out, v, numOfVars);

  for (i=0; i<x; i++) {
    for (j=0; j<y; j++) {
      for (k=0; k<z; k++) {
	for (l=0; l<w; l++) {
	  if (j==0) v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+1*(z*w)+k*w+l+1;
	  if (j==1) v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+0*(z*w)+k*w+l+1;
	  if (j>1)  v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+j*(z*w)+k*w+l+1;
	}
      }
    }
  }
  printPermutationToFile(out, v, numOfVars);

/* Permuting third component. */
  if (z>1) {
    for (i=0; i<x; i++) {
      for (j=0; j<y; j++) {
	for (k=0; k<z; k++) {
	  for (l=0; l<w; l++) {
	    if (k==z-1) {
	      v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+j*(z*w)+0*w+l+1;
	    } else{
	      v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+j*(z*w)+(k+1)*w+l+1;
	    }
	  }
	}
      }
    }
    printPermutationToFile(out, v, numOfVars);

    for (i=0; i<x; i++) {
      for (j=0; j<y; j++) {
	for (k=0; k<z; k++) {
	  for (l=0; l<w; l++) {
	    if (k==0) v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+j*(z*w)+1*w+l+1;
	    if (k==1) v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+j*(z*w)+0*w+l+1;
	    if (k>1)  v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+j*(z*w)+k*w+l+1;
	  }
	}
      }
    }
    printPermutationToFile(out, v, numOfVars);
  }

/* Permuting fourth component. */
  if (w>1) {
    for (i=0; i<x; i++) {
      for (j=0; j<y; j++) {
	for (k=0; k<z; k++) {
	  for (l=0; l<w; l++) {
	    if (l==w-1) {
	      v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+j*(z*w)+k*w+0+1;
	    } else{
	      v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+j*(z*w)+k*w+(l+1)+1;
	    }
	  }
	}
      }
    }
    printPermutationToFile(out, v, numOfVars);

    for (i=0; i<x; i++) {
      for (j=0; j<y; j++) {
	for (k=0; k<z; k++) {
	  for (l=0; l<w; l++) {
	    if (l==0) v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+j*(z*w)+k*w+1+1;
	    if (l==1) v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+j*(z*w)+k*w+0+1;
	    if (l>1)  v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+j*(z*w)+k*w+l+1;
	  }
	}
      }
    }
    printPermutationToFile(out, v, numOfVars);
  }


/* We are left with defining mirror symmetries. */

  if (x==y) {
    for (i=0; i<x; i++) {
      for (j=0; j<y; j++) {
	for (k=0; k<z; k++) {
	  for (l=0; l<w; l++) {
	    v[i*(y*z*w)+j*(z*w)+k*w+l]=j*(y*z*w)+i*(z*w)+k*w+l+1;
	  }
	}
      }
    }
    printPermutationToFile(out, v, numOfVars);
  }

  if (x==z) {
    for (i=0; i<x; i++) {
      for (j=0; j<y; j++) {
	for (k=0; k<z; k++) {
	  for (l=0; l<w; l++) {
	    v[i*(y*z*w)+j*(z*w)+k*w+l]=k*(y*z*w)+j*(z*w)+i*w+l+1;
	  }
	}
      }
    }
    printPermutationToFile(out, v, numOfVars);
  }

  if (x==w) {
    for (i=0; i<x; i++) {
      for (j=0; j<y; j++) {
	for (k=0; k<z; k++) {
	  for (l=0; l<w; l++) {
	    v[i*(y*z*w)+j*(z*w)+k*w+l]=l*(y*z*w)+j*(z*w)+k*w+i+1;
	  }
	}
      }
    }
    printPermutationToFile(out, v, numOfVars);
  }

  if (y==z) {
    for (i=0; i<x; i++) {
      for (j=0; j<y; j++) {
	for (k=0; k<z; k++) {
	  for (l=0; l<w; l++) {
	    v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+k*(z*w)+j*w+l+1;
	  }
	}
      }
    }
    printPermutationToFile(out, v, numOfVars);
  }

  if (y==w) {
    for (i=0; i<x; i++) {
      for (j=0; j<y; j++) {
	for (k=0; k<z; k++) {
	  for (l=0; l<w; l++) {
	    v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+l*(z*w)+k*w+j+1;
	  }
	}
      }
    }
    printPermutationToFile(out, v, numOfVars);
  }

  if ((z>1) && (z==w)) {
    for (i=0; i<x; i++) {
      for (j=0; j<y; j++) {
	for (k=0; k<z; k++) {
	  for (l=0; l<w; l++) {
	    v[i*(y*z*w)+j*(z*w)+k*w+l]=i*(y*z*w)+j*(z*w)+l*w+k+1;
	  }
	}
      }
    }
    printPermutationToFile(out, v, numOfVars);
  }


  return(0);
}
/* ----------------------------------------------------------------- */
