/*
	slice.h is the open source part of vtemlab v0.0 engine,
	providing slice data structures/operations for vtemlab.

	Copyright (C) 2022  Francis Black Lee (Chinese name: ¿ÓœÕ, Li Xian)

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.

	Email: warner323@outlook.com
*/

#ifndef SLICE_H
#define SLICE_H

// mkl dependence
#include <mkl.h>

/***************************************************************/
// data structure:
// Type 1: all the information of one slice is stored in one single SliceNode:
// SliceNode: a slice containing several types of AtomNode's, a node in the SliceList
// This could be applied to projected potential generation using convolution by simply sorting the data:
typedef struct SliceNode
{
	// SliceNode data:
	double LattConstA, LattConstB, SliceDist;
	int AtomNum, * AtomType, PositionNum, * PositionIndices;
	struct SliceNode* NextSlice;
	double* EleProp, * CoordX, * CoordY;
}SliceList;

// SliceAtomList contains all the atoms after lattice expansion:
typedef struct SliceAtomNode
{
	int AtomType;
	double EleProp, CoordX, CoordY;
	struct SliceAtomNode* NextAtom;
}SliceAtomList;

typedef struct TypeAtomNode
{
	int PositionIdx;
	double EleProp, CoordX, CoordY;
	struct TypeAtomNode* NextAtom;
}TypeAtomList;

typedef struct SliceTypeNode
{
	int AtomType, AtomNum;
	double MassNum, DebyeTemp, MeanSquareDisplace;
	struct SliceTypeNode* NextType;
	struct TypeAtomNode* AtomListHead;
}SliceTypeList;

/****************************************************************/
// Operation:
// Create an slice node and fill in the information:
/*Crystal information file format:
The crystal data is saved in a *.txt file;
Format of the file content:
LattConstA	LattConstB	SliceDist
AtomNum
Type1	EleProp1	CoordX1	CorrdY1
Type2	EleProp2	CoordX2	CoordY2
		...
TypeN	ElePropN	CoordXN	CoordYN

Note: AtomNum = N, LattConst = lattice constant, EleProp = elemental proportion, Coord = coordinate
	The Type-EleProp-CoordX-CoordY list does not have to be sorted with repect to the atomic types.
*/
int CreateSliceNode(char* sliceFileName, SliceList*& destSliceNode);

// Sort the data in one SliceNode by the elemental num in ascending order:
// Z denotes elemental number.
void SortSliceNodeByZ(SliceList* destSliceNode);

// Classify all the atoms in a SliceNode by their positions
void ClassifySliceNodeByPosition(SliceList* destSliceNode);

// Create a SliceList composed of a series of SliceNode's:
/* Rules of naming each slice under the input directory:
	the filename must use a lowercased "s" as its prefix and a number
	starting from 1 in ascending order as the rest part of the name,
	file extension is ".txt", i.e. s1.txt, s2.txt, ...
*/
int CreateSliceList(char* sliceDir, SliceList*& destSliceList, int& sliceNum);

// Scan the SliceList: display all the stored data
void ScanSliceList(SliceList* destList, int sliceNum);

// Free a single slice node
void FreeSliceNode(SliceList*& destSliceNode);

// Free the SliceList:
void FreeSliceList(SliceList*& destSliceList);

// Remove the repeated atoms on the boundary (Linked TypeList, Cartesian coordinates):
int RmvRepAtom_LT_Cart(SliceTypeList*& destTypeList, double Lx, double Ly,
	double distTol, double elePropTol);

// Transform the SliceNode to SliceTypeList, note that all the atoms has been 
// sorted by the elemental number (Cartesian coordinates)
int SliceNodeToTypeList_Cart(SliceList* destSliceNode,
	SliceTypeList*& destTypeList, double distTol, double elePropTol);

// Transform the SliceNode to SliceTypeList, note that all the atoms has been 
// sorted by the elemental number (fractional coordinates)
int SliceNodeToTypeList_Frac(SliceList* destSliceNode,
	SliceTypeList*& destTypeList, double distTol, double elePropTol);

// Display the info in a SliceTypeList:
void ScanSliceTypeList(SliceTypeList* destTypeList);

// Save the info in a SliceTypeList:
int SaveSliceTypeList(SliceTypeList* destTypeList, char* filename);

// Free the SliceTypeList:
void FreeSliceTypeList(SliceTypeList*& destTypeList);

#endif // !SLICE_H
