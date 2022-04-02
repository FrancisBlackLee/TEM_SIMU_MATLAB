/*
	thermo.h is the open source part of vtemlab v0.0 engine,
	providing computation of TDS for multislice simulation.

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

#ifndef THERMO_H
#define THERMO_H

#include <mkl.h>

#include "slice.h"

// Supplimentary slice list, containing sliceExpanNum[2], configNum, stackNum,
// posiNum, posiOrientDispArray, rndThermoDispX, rndThermoDispY
typedef struct ThermoSliceNode
{
	int inStackSliceIdx;
	double lattConstA, lattConstB, sliceDist;
	int sliceExpanNum[2];
	int configNum, stackNum;
	int posiNum;
	SliceTypeList* typeListHead;
	double* posiOrientDispArray;
	double** rndThermoDispX, ** rndThermoDispY;

	struct ThermoSliceNode* nextSlice;
}ThermoSliceList;

// Calculate the projected potential for one type of atoms using convolution (atoms of 
// single type): for the case of thermo effect.
// mode == SHIFTED_VEC: shifted ProjPot;
// mode == UNSHIFTED_VEC: unshifted ProjPot;
int MonoAtomProjPot_thermo_conv(MKL_Complex16* projPot, MKL_Complex16* convKer,
	double eleScattParams[4][3], double Lx, double Ly, int Nx, int Ny, int mode);

// Temporary: slice proj pot with thermal vibration (Einstein model)
int SliceComplexProjPot_thermo_conv(ThermoSliceList* thermoSlice,
	int configIdx, int stackIdx, double* scattParams, MKL_Complex16* fxMesh,
	MKL_Complex16* fyMesh, MKL_Complex16* fMesh, MKL_Complex16* convKer,
	MKL_Complex16* projPot, MKL_Complex16* container, int Nx, int Ny,
	int projPotOutMode);

// Bandwidth limit a transmission function:
void BwlSliceTF(MKL_Complex16* sliceTF, double Lx, double Ly, int Nx, int Ny,
	double bwlProp, int tfInMode, int tfOutMode);

// Scan the ThermoSliceList
void ScanThermoSliceList(ThermoSliceList* tsList);

// Free the ThermoSliceList
void FreeThermoSliceList(ThermoSliceList*& tsList);

// Instead of calculating the deviation of thermal displacement with debye temperature,
// it is better to do it with the B factor (or Debye Waller factor), because acquiring an
// accurate debye temperature table is no easier than acquiring a B factor table.
// For more detailed information about B factors, see J.M. Zuo and J.C.H. Spence, Advanced
// Transmission Electron Microscopy.
void FillThermoInfoWithBFac(int typeNum, int* typeArray, double* bFactors_293K,
	SliceTypeList* destTypeList, double simuTemp, int positionNum,
	double* posiOrientThermoDispArray);

// Create the random displacement deviation arrays for both x and y direction
int CreateRandomDispDeviXY(int expanNum[2], int configNum, int stackNum, int posiNum,
	double* posiOrientDispArray, VSLStreamStatePtr stream, double**& rndThermoDispX,
	double**& rndThermoDispY);

// Create ThermoSliceList
int CreateThermoSliceList(SliceList* destSliceList, int coordType, int expanNum[2],
	ThermoSliceList*& tsList, int configNum, int stackNum, int typeNum,
	int* typeArray, double* bFactors_293K, double simuTemp);

// Prepare sample with thermal vibration
int CreateThermoSample(char* folder, int coordType, int expanNum[2], SliceList*& sList,
	ThermoSliceList*& tsList, int configNum, int stackNum, int& sliceNum, int typeNum,
	int* typeArray, double* bFactors_293K, double simuTemp);

#endif // !THERMO_H
