/*
	potential.h is the open source part of vtemlab v0.0 engine,
	providing projected potential calculations for vtemlab.

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

#ifndef POTENTIAL_H
#define POTENTIAL_H

// mkl dependence
#include <mkl.h>

#include "slice.h"

// Load the scattering factor table:
int LoadScattFacTable(double*& scattParams);

// Extract elemental scattering parameters from the ScattFacTable:
int ExtractEleScattParam(int atomType, double* scattParams,
	double eleScattParams[4][3]);

// Atomic scattering factor:
double AtomScatteringFactor(double eleScattParams[4][3], double q);

// Calculate the fft of the atomic projected potential (as complex double-
// precision vector) using atomic scattering factors
int AtomProjPotFFT(MKL_Complex16* projPotFFT, double eleScattParams[4][3],
	double Lx, double Ly, double Nx, double Ny);

// Calculate the projected potential using convolution (atoms of single type):
// SHIFTED_VEC: shifted ProjPot;
// UNSHIFTED_VEC: unshifted ProjPot;
int MonoAtomProjPot_conv(MKL_Complex16* projPot, MKL_Complex16* convKer,
	double eleScattParams[4][3], TypeAtomList* destAtomList, int Nx, int Ny,
	int cellNum[2], double lattConsts[2], int mode);

// Calculate the projected potential using convolution (atoms of multiple types):
int MultiAtomProjPot_conv(double* scattParams, double*& projPot,
	SliceTypeList* destTypeList, int Nx, int Ny, int cellNum[2],
	double lattConsts[2]);

// Save projected potential in a special format:
// first line in the file: Lx (double), Ly (double), Nx (int), Ny (int), SliceDist (double)
// second line in the file: projected potential values as a double-precision real vector.
int SaveProjPot(char* filename, double* projPot, double Lx, double Ly,
	int Nx, int Ny, double sliceDist);

// Load projected potential in a special format:
// first line in the file: Lx (double), Ly (double), Nx (int), Ny (int), SliceDist (double)
// second line in the file: projected potential values as a double-precision real vector.
// Note that argument ProjPot needs to be initialized before being passed to the function:
int LoadProjPot(char* filename, double*& projPot, double& Lx, double& Ly,
	int& Nx, int& Ny, double& sliceDist);

// convert text-image-type projected potential file to binary image:
int ConvertProjPotFromTextToBinary(char* textFilename, char* binFilename);

// convert binary projected potential files to text files.
int ConvertProjPotFromBinaryToText(char* binFilename, char* textFilename);

// save binary-type projected potential file
int SaveBinaryProjPot(char* binFilename, double* projPot, double Lx, double Ly,
	int Nx, int Ny, double sliceDist);

// load binary-type projected potential file
int LoadBinaryProjPot(char* binFilename, double*& projPot, double& Lx, double& Ly,
	int& Nx, int& Ny, double& sliceDist);

// Collectively read in the crystal files and write encrypted projected potential files 
// under the same directory (convolutional way):
int CrysFileToProjPot_conv(char* folder, int expanNum[2], int Nx, int Ny,
	int fileCoordType);

int CrysFileToBinProjPot_conv(char* folder, int expanNum[2], int Nx, int Ny,
	int fileCoordType);

// Decrypt the encrypted projected potential files
int DecryptProjPotFiles(char* folder, int sliceNum);

#endif // !POTENTIAL_H
