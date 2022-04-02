/*
	multislice.h is the open source part of vtemlab v0.0 engine,
	providing related operations for CBED multislice simulation.

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

#ifndef MULTISLICE_H
#define MULTISLICE_H

#include <mkl.h>

#include "slice.h"
#include "thermo.h"

// Fresnel propagation kernel list
typedef struct PropKerNode
{
	int inStackSliceIdx;
	double Lx, Ly;
	int Nx, Ny;
	double sliceDist;
	MKL_Complex16* propKer;
	struct PropKerNode* nextSlice;
}PropKerList;

typedef struct
{
	// Cs3 -- third spherical aberration (mm);
	// Cs5 -- fifth spherical aberration (mm);
	// defocus (Angstrom)
	double Cs3, Cs5, defocus;
}OtfParamSet;

// Create PropKerList
int CreatePropKerList(SliceList* sList, PropKerList*& pkList,
	double wavLen, int expanNum[2], int Nx, int Ny, int pkMode);

// Create PropKerList for slices with uniform spacing
int CreateUniformPropKerList(PropKerList*& pkList, int sliceNum,
	double sliceSpacing, double wavLen, double Lx, double Ly, int Nx, int Ny,
	int pkMode);

// Automatically choose CreatePropKerList or CreateUniformPropKerList to create 
// PropKerList
int AutoCreatePropKerList(SliceList* sList, PropKerList*& pkList, double wavLen,
	int expanNum[2], int Nx, int Ny, int pkMode, bool& sliceDistEqual);

// Free the PropKerList
void FreePropKerList(PropKerList*& pkList, bool sliceDistEqual);

// Calculate the wavelength (Angstrom) of electron beam given the accelerating voltage (KeV)
double EleWavLen(double KeV);

// Calculate the interaction parameter sigma.
/*Note that the unit of electron wavelength is Angstrom, and the unit of sigma
is (Angstrom * kV)^-1, anyway, not all the units of voltage is kV.*/
double InteractionCoefficient(double KeV);

// Generate the objective transfer function in reciprocal space, only 
// including Cs3, Cs5 and df
int ObjTransFunc(MKL_Complex16* otf, OtfParamSet param, double wavLen,
	double Lx, double Ly, int Nx, int Ny);

// Add a circular aperture to OTF:
void AddCircApert(MKL_Complex16* OTF, double WavLen, double numApert,
	double Lx, double Ly, int Nx, int Ny);

// Add an annular aperture to OTF:
void AddAnnularAperture(MKL_Complex16* otf, int otfMode, double wavLen,
	double innerAngle, double outerAngle, double Lx, double Ly, int Nx, int Ny);

// Generate an electron probe:
void GenerateProbe(MKL_Complex16* otf, MKL_Complex16* probe, double xp,
	double yp, double Lx, double Ly, int Nx, int Ny, int outputMode);

void InitFreqXYMesh(MKL_Complex16* fxMesh, double Lx, int Nx,
	MKL_Complex16* fyMesh, double Ly, int Ny);

// save 3D ADF-STEM image
int Save3dImage(char* filename, double** stemImg, int sliceNum,
	int vecLen);

// Convert a depth list to slice index list
void DepthListToSliceIndexList(SliceList* sList, int sliceNumPerStack,
	int maxStackNum, double* depthList, int* sliceIndexList, int depthLevelNum);

int CbedUpdateSliceTF(MKL_Complex16* tf, int tfMode, double interCoeff,
	ThermoSliceList* tsList, int configIdx, int stackIdx, double bwlProp,
	double* scattFac, MKL_Complex16* fxMesh, MKL_Complex16* fyMesh,
	MKL_Complex16* fMesh, MKL_Complex16* convKer, MKL_Complex16* container,
	int Nx, int Ny);

// Dedicated multislice kernel for CBED simulation with TDS, cbedMat output in
// w_oMode
int CbedMultisliceKernel(MKL_Complex16* wave, int wIMode, int wOMode,
	double voltage, double bwlProp, ThermoSliceList* tsList,
	PropKerList* pkList, int selectSliceNum, int* selectSliceIndices,
	double** cbedMat, int Nx, int Ny);

int CbedTdsKernel(OtfParamSet otfParam, double voltage, double numApert,
	ThermoSliceList* tsList, PropKerList* pkList, double xp, double yp,
	double bwlProp, int Nx, int Ny, int selectSliceNum, int* selectSliceIndices,
	char* cbedFilename);

#endif // !MULTISLICE_H
