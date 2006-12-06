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

void printVersionInfo();
void printVector(vector, int);
void printListVector(listVector*, int);
void printVectorToFile(FILE*, vector, int);
void printListVectorToFile(char*, listVector*, int);
void printTransposedListVectorToFile(char*, listVector*, int);
void printVectorToFileWithBrackets(FILE*, vector, int);
void printBinomialToFile(FILE*, vector, int, char**);
void printListBinomialsToFile(char*, listVector*, int, char**);
void printMonomialToFile(FILE*, vector, int, char**);
void printListMonomialsAndBinomialsToFile(char*, listVector*, int, char**);
void printMatrix(vector, int, int);
void printVectorToFileMaple(FILE*, vector, int);
void printListVectorMaple(char*, listVector*, int);
void printListVectorMacaulay2(char*, listVector*, int);
void print3wayTables(char*, listVector*, int, int, int, int);
void printL1NormOfListVector(listVector*, int);
void printListVectorWithGivenDegreeToFile(char*, listVector*, int, int);
void printListVectorWithGivenNonzeroEntryToFile(char*, listVector*, int, int);
void writeResult(listVector*, int, char*, char*, int);
void printListRepresentativesToFile(char*, listOrbit*, int);
void printRationalVector(rationalVector*, int);
void printRationalVectorToFileWithoutBrackets(FILE*, rationalVector*, int);
