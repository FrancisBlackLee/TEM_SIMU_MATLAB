/*
	slice.cpp is the open source part of vtemlab v0.0 engine,
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

#include <stdio.h>
#include <stdlib.h>
#include <io.h>
#include <direct.h>
#include <math.h>
#include <string.h>

#include "mathconst.h"
#include "mode.h"
#include "slice.h"
#include "status.h"
#include "vtemlabio.h"


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

Note: AtomNum = N, LattConst = lattice constant, EleProp = elemental proportion, 
Coord = coordinate
The Type-EleProp-CoordX-CoordY list does not have to be sorted with repect to 
the atomic types.
*/
int CreateSliceNode(char* sliceFileName, SliceList*& destSliceNode)
{
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, sliceFileName, "r");
	if (err)
	{
		PrintErrorMsg("CreateSliceNode", VTEMLAB_EFILER);
		return VTEMLAB_EFILER;
	}
	destSliceNode = (SliceList*)malloc(sizeof(SliceList));
	// initializing:
	// two dimensional lattice constants and slice distance to the next slice:
	fscanf_s(destFile, "%lf", &destSliceNode->LattConstA);
	fscanf_s(destFile, "%lf", &destSliceNode->LattConstB);
	fscanf_s(destFile, "%lf", &destSliceNode->SliceDist);

	// read number of atoms:
	fscanf_s(destFile, "%d", &destSliceNode->AtomNum);
	int atomNum = destSliceNode->AtomNum;

	// allocate memory to *AtomType, *EleProp, *CoordX and *CoordY:
	destSliceNode->AtomType = (int*)malloc(atomNum * sizeof(int));
	destSliceNode->PositionIndices = (int*)malloc(atomNum * sizeof(int));
	destSliceNode->EleProp = (double*)mkl_malloc(atomNum * sizeof(double), 64);
	destSliceNode->CoordX = (double*)mkl_malloc(atomNum * sizeof(double), 64);
	destSliceNode->CoordY = (double*)mkl_malloc(atomNum * sizeof(double), 64);
	if ((destSliceNode->AtomType == NULL) || (destSliceNode->EleProp == NULL)
		|| (destSliceNode->CoordX == NULL || (destSliceNode->CoordY == NULL)))
	{
		PrintErrorMsg("CreateSliceNode", VTEMLAB_ENOMEM);
		free(destSliceNode->AtomType);
		mkl_free(destSliceNode->EleProp);
		mkl_free(destSliceNode->CoordX);
		mkl_free(destSliceNode->CoordY);
		free(destSliceNode);
		return VTEMLAB_ENOMEM;
	}

	// Read in the rest of the data:
	for (int atomIdx = 0; atomIdx < atomNum; atomIdx++)
	{
		destSliceNode->PositionIndices[atomIdx] = atomIdx;// initialization
		fscanf_s(destFile, "%d", &destSliceNode->AtomType[atomIdx]);
		fscanf_s(destFile, "%lf", &destSliceNode->EleProp[atomIdx]);
		fscanf_s(destFile, "%lf", &destSliceNode->CoordX[atomIdx]);
		fscanf_s(destFile, "%lf", &destSliceNode->CoordY[atomIdx]);
	}
	fclose(destFile);

	// set the reference of the NextSlice to be NULL:
	destSliceNode->NextSlice = NULL;

	// sort the SliceNode data by elemental number in ascending order:
	SortSliceNodeByZ(destSliceNode);
	ClassifySliceNodeByPosition(destSliceNode);

	return VTEMLAB_SUCCESS;
}


// Sort the data in one SliceNode by the elemental number in ascending order:
// Z denotes elemental number.
void SortSliceNodeByZ(SliceList* destSliceNode)
{
	int atomNum, i, j, tmpAtomZ;
	double tmpEleProp, tmpCoordX, tmpCoordY;
	atomNum = destSliceNode->AtomNum;
	for (i = 0; i < atomNum - 1; i++)
	{
		for (j = i + 1; j < atomNum; j++)
		{
			if (destSliceNode->AtomType[i] > destSliceNode->AtomType[j])
			{
				// switch atomic types:
				tmpAtomZ = destSliceNode->AtomType[i];
				destSliceNode->AtomType[i] = destSliceNode->AtomType[j];
				destSliceNode->AtomType[j] = tmpAtomZ;

				// switch elemental proportions:
				tmpEleProp = destSliceNode->EleProp[i];
				destSliceNode->EleProp[i] = destSliceNode->EleProp[j];
				destSliceNode->EleProp[j] = tmpEleProp;

				// switch atomic x coordinate:
				tmpCoordX = destSliceNode->CoordX[i];
				destSliceNode->CoordX[i] = destSliceNode->CoordX[j];
				destSliceNode->CoordX[j] = tmpCoordX;

				// switch atomic y coordinate:
				tmpCoordY = destSliceNode->CoordY[i];
				destSliceNode->CoordY[i] = destSliceNode->CoordY[j];
				destSliceNode->CoordY[j] = tmpCoordY;
			}
		}
	}
}

// Classify all the atoms in a SliceNode by their positions
void ClassifySliceNodeByPosition(SliceList* destSliceNode)
{
	int currentPosition = 0, i, j, atomNum;
	atomNum = destSliceNode->AtomNum;
	// initialize CheckList as all zeros, 0 for unchecked, 1 for checked:
	int* checkList = (int*)malloc(destSliceNode->AtomNum * sizeof(int));
	for (i = 0; i < atomNum; i++)
		checkList[i] = 0;

	// Start checking:
	double tolerance = 1.0e-4;
	for (i = 0; i < atomNum; i++)
	{
		if (checkList[i] == 0)
		{
			destSliceNode->PositionIndices[i] = currentPosition;
			checkList[i] = 1;
			for (j = 0; j < atomNum; j++)
			{
				if ((checkList[j] == 0)
					&& (fabs(destSliceNode->CoordX[j] - destSliceNode->CoordX[i]) < tolerance)
					&& (fabs(destSliceNode->CoordY[j] - destSliceNode->CoordY[i]) < tolerance))
				{
					destSliceNode->PositionIndices[j] = currentPosition;
					checkList[j] = 1;
				}
			}
			currentPosition++;
		}
	}
	destSliceNode->PositionNum = currentPosition;
	free(checkList);
}

// Create a SliceList composed of a series of SliceNode's:
/* Rules of naming each slice under the input directory:
	the filename must use a lowercased "s" as its prefix and a number
	starting from 1 in ascending order as the rest part of the name,
	file extension is ".txt", i.e. s1.txt, s2.txt, ...
*/
int CreateSliceList(char* sliceDir, SliceList*& destSliceList, int& sliceNum)
{
	// Allocate memory to the SliceListHead and declare a pointer and container for SliceNode's:
	SliceList* prevNode, * nextNode;
	destSliceList = (SliceList*)malloc(sizeof(SliceList));
	prevNode = destSliceList;

	// initializing SliceNum as 0:
	sliceNum = 0;
	errno_t err = false;
	int errorCode; // operation status denoting whether an operation is successful.
	while (!err)
	{
		char filename[FILENAME_MAX];
		snprintf(filename, FILENAME_MAX, "%s\\s%d.txt", sliceDir, sliceNum + 1);

		err = _access_s(filename, 04);
		if (!err)
		{
			errorCode = CreateSliceNode(filename, nextNode);
			if (errorCode)
			{
				PrintErrorMsg("CreateSliceList", VTEMLAB_FAILURE);
				return VTEMLAB_FAILURE;
			}
			prevNode->NextSlice = nextNode;
			prevNode = nextNode;
			sliceNum++;
		}
	}
	prevNode->NextSlice = NULL;

	return VTEMLAB_SUCCESS;
}

// Scan the SliceList: display all the stored data
void ScanSliceList(SliceList* destList, int sliceNum)
{
	SliceList* dispNode = destList->NextSlice;
	int i, j;
	for (i = 0; i < sliceNum; i++)
	{
		printf("slice %d:\n\n", i + 1);
		printf("%.4f\t%.4f\t%.4f\n", dispNode->LattConstA, dispNode->LattConstB, dispNode->SliceDist);
		printf("%d\n", dispNode->AtomNum);
		for (j = 0; j < dispNode->AtomNum; j++)
			printf("%d\t%.4f\t%.4f\t%.4f\n", dispNode->AtomType[j], dispNode->EleProp[j], dispNode->CoordX[j], dispNode->CoordY[j]);
		printf("\n");
		dispNode = dispNode->NextSlice;
	}
}


// Free a single slice node
void FreeSliceNode(SliceList*& destSliceNode)
{
	if (destSliceNode != NULL)
	{
		if (destSliceNode->AtomType != NULL)
			free(destSliceNode->AtomType);
		if (destSliceNode->PositionIndices != NULL)
			free(destSliceNode->PositionIndices);
		if (destSliceNode->EleProp != NULL)
			mkl_free(destSliceNode->EleProp);
		if (destSliceNode->CoordX != NULL)
			mkl_free(destSliceNode->CoordX);
		if (destSliceNode->CoordY != NULL)
			mkl_free(destSliceNode->CoordY);

		free(destSliceNode);
	}
}


// Free the SliceList:
void FreeSliceList(SliceList*& destSliceList)
{
	SliceList* preNode = destSliceList, * postNode = destSliceList->NextSlice;
	while (postNode != NULL)
	{
		free(preNode);
		free(postNode->AtomType);
		free(postNode->PositionIndices);
		mkl_free(postNode->EleProp);
		mkl_free(postNode->CoordX);
		mkl_free(postNode->CoordY);
		preNode = postNode;
		postNode = preNode->NextSlice;
	}
	free(preNode);
}


// Remove the repeated atoms on the boundary (Linked TypeList, Cartesian coordinates):
int RmvRepAtom_LT_Cart(SliceTypeList*& destTypeList, double Lx, double Ly, 
	double distTol, double elePropTol)
{
	if (destTypeList->NextType == NULL)
	{
		PrintErrorMsg("RmvRepAtom_LT_Cart", VTEMLAB_ENULLPTR);
		return VTEMLAB_ENULLPTR;
	}

	/*Note that now the AtomContainer contains the atom found to be deleted,
	and the AtomPointer contains the previous atom of the atom to be deleted,
	an additional AtomNode is used to point to the atom being compared with*/
	SliceTypeList* tmpTypeNode = destTypeList->NextType;
	TypeAtomList* tmpAtomNode, * prevAtomNode, * postAtomNode;

	while (tmpTypeNode != NULL)
	{
		tmpAtomNode = (tmpTypeNode->AtomListHead)->NextAtom;
		prevAtomNode = tmpAtomNode, postAtomNode = tmpAtomNode->NextAtom;

		// Search along x direction:
		while (tmpAtomNode != NULL)
		{
			while (postAtomNode != NULL)
			{
				if ((fabs(postAtomNode->EleProp - tmpAtomNode->EleProp) < elePropTol) &&
					(fabs(postAtomNode->CoordX - tmpAtomNode->CoordX) < distTol) &&
					(fabs(fabs(postAtomNode->CoordY - tmpAtomNode->CoordY) - Ly) < distTol))
				{
					prevAtomNode->NextAtom = postAtomNode->NextAtom;
					free(postAtomNode);
					postAtomNode = prevAtomNode->NextAtom;
				}
				else
				{
					prevAtomNode = postAtomNode;
					postAtomNode = prevAtomNode->NextAtom;
				}
			}
			tmpAtomNode = tmpAtomNode->NextAtom;
			prevAtomNode = tmpAtomNode;
			if (prevAtomNode != NULL)
				postAtomNode = prevAtomNode->NextAtom;
		}

		tmpAtomNode = (tmpTypeNode->AtomListHead)->NextAtom;
		prevAtomNode = tmpAtomNode, postAtomNode = tmpAtomNode->NextAtom;

		// Search along y direction:
		while (tmpAtomNode != NULL)
		{
			while (postAtomNode != NULL)
			{
				if ((fabs(postAtomNode->EleProp - tmpAtomNode->EleProp) < elePropTol) &&
					(fabs(fabs(postAtomNode->CoordX - tmpAtomNode->CoordX) - Lx) < distTol) &&
					(fabs(postAtomNode->CoordY - tmpAtomNode->CoordY) < distTol))
				{
					prevAtomNode->NextAtom = postAtomNode->NextAtom;
					free(postAtomNode);
					postAtomNode = prevAtomNode->NextAtom;
				}
				else
				{
					prevAtomNode = postAtomNode;
					postAtomNode = prevAtomNode->NextAtom;
				}
			}
			tmpAtomNode = tmpAtomNode->NextAtom;
			prevAtomNode = tmpAtomNode;
			if (prevAtomNode != NULL)
				postAtomNode = prevAtomNode->NextAtom;
		}

		tmpAtomNode = (tmpTypeNode->AtomListHead)->NextAtom;
		prevAtomNode = tmpAtomNode, postAtomNode = tmpAtomNode->NextAtom;

		// Search along diagonal direction:
		while (tmpAtomNode != NULL)
		{
			while (postAtomNode != NULL)
			{
				if ((fabs(postAtomNode->EleProp - tmpAtomNode->EleProp) < elePropTol) &&
					(fabs(fabs(postAtomNode->CoordX - tmpAtomNode->CoordX) - Lx) < distTol) &&
					(fabs(fabs(postAtomNode->CoordY - tmpAtomNode->CoordY) - Ly) < distTol))
				{
					prevAtomNode->NextAtom = postAtomNode->NextAtom;
					free(postAtomNode);
					postAtomNode = prevAtomNode->NextAtom;
				}
				else
				{
					prevAtomNode = postAtomNode;
					postAtomNode = prevAtomNode->NextAtom;
				}
			}
			tmpAtomNode = tmpAtomNode->NextAtom;
			prevAtomNode = tmpAtomNode;
			if (prevAtomNode != NULL)
				postAtomNode = prevAtomNode->NextAtom;
		}

		tmpTypeNode = tmpTypeNode->NextType;
	}

	return VTEMLAB_SUCCESS;
}

// Transform the SliceNode to SliceTypeList, note that all the atoms has been 
// sorted by the elemental number (Cartesian coordinates)
int SliceNodeToTypeList_Cart(SliceList* destSliceNode, 
	SliceTypeList*& destTypeList, double distTol, double elePropTol)
{
	destTypeList = (SliceTypeList*)malloc(sizeof(SliceTypeList));
	if (destTypeList == NULL)
	{
		PrintErrorMsg("SliceNodeToTypeList_Cart", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}
	int i, prevType, errorCode;// Previous type.
	SliceTypeList* prevTypeNode, * postTypeNode;
	prevTypeNode = destTypeList;
	postTypeNode = (SliceTypeList*)malloc(sizeof(SliceTypeList));

	TypeAtomList* prevAtomNode, * postAtomNode;
	prevAtomNode = (TypeAtomList*)malloc(sizeof(TypeAtomList));

	postTypeNode->AtomListHead = prevAtomNode;
	prevType = destSliceNode->AtomType[0];
	postTypeNode->AtomType = prevType;
	postTypeNode->DebyeTemp = 645.0;// default
	postTypeNode->MassNum = 28.0;// default
	postTypeNode->MeanSquareDisplace = 1.0e-4;// default
	postTypeNode->AtomNum = 0;
	prevTypeNode->NextType = postTypeNode;
	prevTypeNode = postTypeNode;

	for (i = 0; i < destSliceNode->AtomNum; i++)
	{
		if (destSliceNode->AtomType[i] == prevType)
		{
			postAtomNode = (TypeAtomList*)malloc(sizeof(TypeAtomList));
			postAtomNode->PositionIdx = destSliceNode->PositionIndices[i];
			postAtomNode->EleProp = destSliceNode->EleProp[i];
			postAtomNode->CoordX = destSliceNode->CoordX[i];
			postAtomNode->CoordY = destSliceNode->CoordY[i];
			prevAtomNode->NextAtom = postAtomNode;
			prevAtomNode = postAtomNode;
			postTypeNode->AtomNum++;
		}
		else
		{
			prevAtomNode->NextAtom = NULL;
			// Create a new type:
			postTypeNode = (SliceTypeList*)malloc(sizeof(SliceTypeList));
			prevType = destSliceNode->AtomType[i];
			postTypeNode->AtomType = prevType;
			prevAtomNode = (TypeAtomList*)malloc(sizeof(TypeAtomList));
			postTypeNode->AtomListHead = prevAtomNode;
			postTypeNode->DebyeTemp = 300.0;// default
			postTypeNode->MassNum = 100.0;// default
			postTypeNode->MeanSquareDisplace = 1.0e-4;// default

			postAtomNode = (TypeAtomList*)malloc(sizeof(TypeAtomList));
			postAtomNode->PositionIdx = destSliceNode->PositionIndices[i];
			postAtomNode->EleProp = destSliceNode->EleProp[i];
			postAtomNode->CoordX = destSliceNode->CoordX[i];
			postAtomNode->CoordY = destSliceNode->CoordY[i];
			prevAtomNode->NextAtom = postAtomNode;
			prevAtomNode = postAtomNode;

			postTypeNode->AtomNum = 1;

			prevTypeNode->NextType = postTypeNode;
			prevTypeNode = postTypeNode;
		}
	}
	prevAtomNode->NextAtom = NULL;
	prevTypeNode->NextType = NULL;

	// Remove the repeated atoms on the boundary:
	errorCode = RmvRepAtom_LT_Cart(destTypeList, destSliceNode->LattConstA, destSliceNode->LattConstB, distTol, elePropTol);
	if (errorCode)
	{
		PrintErrorMsg("SliceNodeToTypeList_Cart", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	return VTEMLAB_SUCCESS;
}

// Transform the SliceNode to SliceTypeList, note that all the atoms has been 
// sorted by the elemental number (fractional coordinates)
int SliceNodeToTypeList_Frac(SliceList* destSliceNode, 
	SliceTypeList*& destTypeList, double distTol, double elePropTol)
{
	destTypeList = (SliceTypeList*)malloc(sizeof(SliceTypeList));
	if (destTypeList == NULL)
	{
		PrintErrorMsg("SliceNodeToTypeList_Frac", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}
	int i, prevType, errorCode;// Previous type.
	SliceTypeList* prevTypeNode, * postTypeNode;
	prevTypeNode = destTypeList;
	postTypeNode = (SliceTypeList*)malloc(sizeof(SliceTypeList));

	TypeAtomList* prevAtomNode, * postAtomNode;
	prevAtomNode = (TypeAtomList*)malloc(sizeof(TypeAtomList));

	postTypeNode->AtomListHead = prevAtomNode;
	prevType = destSliceNode->AtomType[0];
	postTypeNode->AtomType = prevType;
	postTypeNode->DebyeTemp = 300.0;// default
	postTypeNode->MassNum = 100.0;// default
	postTypeNode->MeanSquareDisplace = 1.0e-4;// default
	postTypeNode->AtomNum = 0;
	prevTypeNode->NextType = postTypeNode;
	prevTypeNode = postTypeNode;

	for (i = 0; i < destSliceNode->AtomNum; i++)
	{
		if (destSliceNode->AtomType[i] == prevType)
		{
			postAtomNode = (TypeAtomList*)malloc(sizeof(TypeAtomList));
			postAtomNode->PositionIdx = destSliceNode->PositionIndices[i];
			postAtomNode->EleProp = destSliceNode->EleProp[i];
			postAtomNode->CoordX = destSliceNode->CoordX[i] * destSliceNode->LattConstA;
			postAtomNode->CoordY = destSliceNode->CoordY[i] * destSliceNode->LattConstB;
			prevAtomNode->NextAtom = postAtomNode;
			prevAtomNode = postAtomNode;
			postTypeNode->AtomNum++;
		}
		else
		{
			prevAtomNode->NextAtom = NULL;
			// Create a new type:
			postTypeNode = (SliceTypeList*)malloc(sizeof(SliceTypeList));
			prevType = destSliceNode->AtomType[i];
			postTypeNode->AtomType = prevType;
			prevAtomNode = (TypeAtomList*)malloc(sizeof(TypeAtomList));
			postTypeNode->AtomListHead = prevAtomNode;
			postTypeNode->DebyeTemp = 300.0;// default
			postTypeNode->MassNum = 100.0;// default
			postTypeNode->MeanSquareDisplace = 1.0e-4;// default

			postAtomNode = (TypeAtomList*)malloc(sizeof(TypeAtomList));
			postAtomNode->PositionIdx = destSliceNode->PositionIndices[i];
			postAtomNode->EleProp = destSliceNode->EleProp[i];
			postAtomNode->CoordX = destSliceNode->CoordX[i] * destSliceNode->LattConstA;
			postAtomNode->CoordY = destSliceNode->CoordY[i] * destSliceNode->LattConstB;
			prevAtomNode->NextAtom = postAtomNode;
			prevAtomNode = postAtomNode;

			postTypeNode->AtomNum = 1;

			prevTypeNode->NextType = postTypeNode;
			prevTypeNode = postTypeNode;
		}
	}
	prevAtomNode->NextAtom = NULL;
	prevTypeNode->NextType = NULL;

	// Remove the repeated atoms on the boundary:
	errorCode = RmvRepAtom_LT_Cart(destTypeList, destSliceNode->LattConstA, destSliceNode->LattConstB, distTol, elePropTol);
	if (errorCode)
	{
		PrintErrorMsg("SliceNodeToTypeList_Cart", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	return VTEMLAB_SUCCESS;
}

// Display the info in a SliceTypeList:
void ScanSliceTypeList(SliceTypeList* destTypeList)
{
	SliceTypeList* dispTypeNode;
	TypeAtomList* dispAtomNode;
	dispTypeNode = destTypeList->NextType;
	while (dispTypeNode != NULL)
	{
		printf("Z = %d\n u2 = %.4f\n", dispTypeNode->AtomType, dispTypeNode->MeanSquareDisplace);
		dispAtomNode = dispTypeNode->AtomListHead;
		dispAtomNode = dispAtomNode->NextAtom;
		while (dispAtomNode != NULL)
		{
			printf("\t%d\t%.4f\t%.4f\t%.4f\n", dispAtomNode->PositionIdx, dispAtomNode->EleProp, dispAtomNode->CoordX, dispAtomNode->CoordY);
			dispAtomNode = dispAtomNode->NextAtom;
		}
		dispTypeNode = dispTypeNode->NextType;
	}
}

// Save the info in a SliceTypeList:
int SaveSliceTypeList(SliceTypeList* destTypeList, char* filename)
{
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "w+");
	if (err)
	{
		PrintErrorMsg("SaveSliceTypeList", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}

	SliceTypeList* tmpTypeNode;
	TypeAtomList* tmpAtomNode;

	// Start saving:
	tmpTypeNode = destTypeList->NextType;
	while (tmpTypeNode != NULL)
	{
		tmpAtomNode = tmpTypeNode->AtomListHead->NextAtom;
		while (tmpAtomNode != NULL)
		{
			fprintf_s(destFile, "%d\t%lf\t%lf\t%lf\n", tmpTypeNode->AtomType,
				tmpAtomNode->EleProp, tmpAtomNode->CoordX, tmpAtomNode->CoordY);
			tmpAtomNode = tmpAtomNode->NextAtom;
		}
		tmpTypeNode = tmpTypeNode->NextType;
	}
	fclose(destFile);

	return VTEMLAB_SUCCESS;
}

// Free the SliceTypeList:
void FreeSliceTypeList(SliceTypeList*& destTypeList)
{
	SliceTypeList* prevTypeNode, * postTypeNode;
	TypeAtomList* prevAtomNode, * postAtomNode;
	prevTypeNode = destTypeList, postTypeNode = destTypeList->NextType;
	while (postTypeNode != NULL)
	{
		free(prevTypeNode);
		prevAtomNode = postTypeNode->AtomListHead;
		postAtomNode = prevAtomNode->NextAtom;
		while (postAtomNode != NULL)
		{
			free(prevAtomNode);
			prevAtomNode = postAtomNode;
			postAtomNode = prevAtomNode->NextAtom;
		}
		free(prevAtomNode);
		prevTypeNode = postTypeNode;
		postTypeNode = prevTypeNode->NextType;
	}
	free(prevTypeNode);
}