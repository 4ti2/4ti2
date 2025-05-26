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
/* Generation of symmetry groups of 2- and 3-way tables              */
/*                                                                   */
/* Author   : Raymond Hemmecke                                       */
/*                                                                   */
/* Created    : 26-AUG-02                                            */
/*                                                                   */
/* ----------------------------------------------------------------- */
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include "myheader.h"
#include "print.h"
#include "vector.h"

#include <getopt.h>
#include "banner.h"

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
  printf("usage: gensymm [--options] A B C D FILENAME\n"
	 "\n"
	 "Computes the generators for the symmetry group acting on 4-way tables\n"
	 "with 3-marginals. By putting one side length to 1, this includes\n"
	 "3-way tables with 2-marginals.\n"
	 "\n"
	 "Options:\n"
	 " -q, --quiet       No output is written to the screen\n"
	 "\n"
	 "Output file:\n"
	 " FILENAME.sym      generators for the symmetry group\n"
	 "\n"
	 "Example:  Consider the problem of 3x3x3 tables with 2-marginals. Calling\n"
	 "  gensymm 3 3 3 1 333\n"
	 "produces the file '333.sym' containing the following lines.\n"
	 "\n"
	 "9 27\n"
	 "10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 1 2 3 4 5 6 7 8 9 \n"
	 "10 11 12 13 14 15 16 17 18 1 2 3 4 5 6 7 8 9 19 20 21 22 23 24 25 26 27 \n"
	 "4 5 6 7 8 9 1 2 3 13 14 15 16 17 18 10 11 12 22 23 24 25 26 27 19 20 21 \n"
	 "4 5 6 1 2 3 7 8 9 13 14 15 10 11 12 16 17 18 22 23 24 19 20 21 25 26 27 \n"
	 "2 3 1 5 6 4 8 9 7 11 12 10 14 15 13 17 18 16 20 21 19 23 24 22 26 27 25 \n"
	 "2 1 3 5 4 6 8 7 9 11 10 12 14 13 15 17 16 18 20 19 21 23 22 24 26 25 27 \n"
	 "1 2 3 10 11 12 19 20 21 4 5 6 13 14 15 22 23 24 7 8 9 16 17 18 25 26 27 \n"
	 "1 10 19 4 13 22 7 16 25 2 11 20 5 14 23 8 17 26 3 12 21 6 15 24 9 18 27 \n"
	 "1 4 7 2 5 8 3 6 9 10 13 16 11 14 17 12 15 18 19 22 25 20 23 26 21 24 27 \n"
	 "\n");
}

/* ----------------------------------------------------------------- */
int gensymm_main(int argc, char *argv[]) {
  int a,i,j,k,l,x,y,z,w,numOfVars,numOfGenerators,infoLevel;
  vector v;
  const char *fileName=NULL;
	char *outFileName=NULL;
  FILE *out;

  int optc;

  setbuf(stdout,0);

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

  if (optind != argc - 5) {
    print_usage();
    exit(1);
  }

  if (infoLevel>-1) {
    printVersionInfo();
  }

/* We require that all 1's come last in the tuple (x,y,z,w). */

  fileName=argv[argc-1];
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

  myxasprintf(&outFileName,"%s.sym",fileName);

  errno=0;
  if (!(out = fopen(outFileName,"w"))) {
    printf("Error opening generator file for output (%s).\n",strerror(errno));
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

  fclose(out);

  free(outFileName);

  return(0);
}
/* ----------------------------------------------------------------- */
