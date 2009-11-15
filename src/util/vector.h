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

int gcd(int, int);
int lcm(int, int);
vector posVector(vector, int);
vector negVector(vector, int);
vector removeGCDfromVector(vector, int);
vector supportOfVector(vector, int, int);
vector negateSupportVector(vector, int);
vector positiveSupportOfVector(vector v, int, int);
vector negativeSupportOfVector(vector v, int, int);
int normOfBinaryVector(vector v, int);
int positiveNormOfVector(vector v, int);
int negativeNormOfVector(vector v, int);
int hasCommonFactor(vector, vector, int);
listVector* appendVectorToListVector(vector, listVector*);
void freeVector(vector);
void freeListVector(listVector*);
void freeAllOfListVector(listVector*);
vector createVector(int);
vector* createArrayVector(int);
listVector** createArrayListVector(int);
int isZeroVector(vector, int);
int isAllOneVector(vector, int);
int isVectorEqualToVector(vector, vector, int);
int isVectorEqualToNegativeVector(vector, vector, int);
vector copyVector(vector, int);
vector negativeVector(vector, int);
vector lexPositiveVector(vector, int);
int isVectorLexPositive(vector, int);
vector addVector(vector, vector, int);
vector addZeroOneVector(vector, vector, int);
vector subVector(vector, vector, int);
listVector* createListVector(vector);
int lengthListVector(listVector*);
listVector* copyListVector(listVector*, int);
listVector* copyListVectorWithoutVectors(listVector*);
vector permuteVector(vector, vector, int);
listVector* permuteListVector(listVector*, vector, int);
vector rePermuteVector(vector, vector, vector, int);
listVector* rePermuteListVector(listVector*, vector, int);
int isVectorInListVector(vector, listVector*, int);
int isNegativeVectorInListVector(vector, listVector*, int);
listVector* updateBasis(listVector*, listVector*);
int isBelowUpperBounds(vector, vector, int);
listVector* appendListVectorToListVector(listVector*, listVector*);
vector addMultipleVector(vector, int, vector, int);
vector subMultipleVector(vector, int, vector, int);
vector permuteMatrix(vector, vector, int, int);
vector permuteTransposedMatrix(vector, vector, int, int);
vector transpose(vector,int,int);
vector matrixTimesVector(listVector*, vector, int, int);
int hasSmallerSupport(vector, vector, int);
int isCircuit(listVector*, vector, int);
listVector* extractCircuits(listVector*, int);
int compareVectorsByLex(vector, vector, int);
listVector* combineOrderedListVectors(listVector*, listVector*, int);
listVector* swapColumnsInListVector(listVector*, int, int, int);
listVector* extractVectorsBelowUpperBounds(listVector*, vector, int);
int normOfVector(vector, int);
int signedNormOfVector(vector, int);
int weightedNormOfVector(vector, vector, int);
int minimalNormInListVector(listVector*, int);
int maximalNormInListVector(listVector*, int);
listVector* projectListVectorDown(listVector*, int, int);
listVector* reduceLastComponent(listVector*, vector, int);
listVector* appendNewComponentToListVector(listVector*, listVector*, 
					   int, int, int);
listVector* appendRemainingComponentsToListVector(listVector*, listVector*, 
						  int, int);
int isZeroOneVector(vector, int);
listVector* extractZeroOneVectors(listVector*, int);
listVector* extractZeroOneVectorsInLastComponent(listVector*, int);
int isNonnegativeVector(vector, int);
listVector* extractNonnegativeVectors(listVector*, int);
void swapGraver(vector*, int, int);
void liftGraver(vector*, vector*, int, int);
int removeFirstHeapElement(vector*, vector*, int, int);
int updateHeapGraver(vector, vector*, vector, vector*, int, int);
int isInSameOrthant(vector, vector, int);
int isInOppositeOrthant(vector, vector, int);
int isIdentityPermutation(vector, int);
int isVectorEqualToPermutedVector(vector, vector, int);
int isVectorEqualToNegativePermutedVector(vector, vector, int);
vector subMultiplePermutedVector(vector, int, vector, vector, int);
vector multiplyPermutation(vector, vector, int);
listVector* vTimesS(vector, listVector*, int, int);
listVector* generateSymmetryGroup(listVector*, int);
vector orientVector(vector, int, int);
vector orientVectorAccordingToOrdering(vector, vector, int);
int isVectorCorrectlyOriented(vector, int, int);
int dotProduct(vector, vector, int);
int isZeroVectorInListVector(listVector*, int);
