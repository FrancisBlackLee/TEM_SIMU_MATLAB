/*
	potential.cpp is the open source part of vtemlab v0.0 engine,
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

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "mathconst.h"
#include "mode.h"
#include "potential.h"
#include "slice.h"
#include "status.h"
#include "vtemlabio.h"
#include "vtemlabmath.h"


// Load the scattering factor table:
int LoadScattFacTable(double*& scattParams)
{
	scattParams = (double*)mkl_malloc(103 * 12 * sizeof(double), 64);
	char filename[] = "ScattParams.txt";

	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "r");
	if (err)
	{
		PrintErrorMsg("LoadScattFacTable", VTEMLAB_EFILER);
		mkl_free(scattParams);
		return VTEMLAB_EFILER;
	}

	int i;
	for (i = 0; i < 103 * 12; i++)
		fscanf_s(destFile, "%lf", &scattParams[i]);

	fclose(destFile);

	return VTEMLAB_SUCCESS;
}


// Extract elemental scattering parameters from the ScattFacTable:
int ExtractEleScattParam(int atomType, double* scattParams, 
	double eleScattParams[4][3])
{
	if (scattParams == NULL)
	{
		PrintErrorMsg("ExtractEleScattParam", VTEMLAB_ENULLPTR);
		return VTEMLAB_ENULLPTR;
	}

	int eleStartIdx = 12 * (atomType - 1);
	for (int colIdx = 0; colIdx < 3; colIdx++)
	{
		eleScattParams[0][colIdx] = scattParams[eleStartIdx + 0 + 2 * colIdx];
		eleScattParams[1][colIdx] = scattParams[eleStartIdx + 1 + 2 * colIdx];
		eleScattParams[2][colIdx] = scattParams[eleStartIdx + 6 + 2 * colIdx];
		eleScattParams[3][colIdx] = scattParams[eleStartIdx + 7 + 2 * colIdx];
	}

	return VTEMLAB_SUCCESS;
}


// Atomic scattering factor:
double AtomScatteringFactor(double eleScattParams[4][3], double q)
{
	double sf = 0.0;
	for (int i = 0; i < 3; i++)
	{
		sf += eleScattParams[0][i] / (q * q + eleScattParams[1][i]) +
			eleScattParams[2][i] * exp(-eleScattParams[3][i] * q * q);
	}

	return sf;
}


// Calculate the fft of the atomic projected potential (as complex double-
// precision vector) using atomic scattering factors
int AtomProjPotFFT(MKL_Complex16* projPotFFT, double eleScattParams[4][3],
	double Lx, double Ly, double Nx, double Ny)
{
	if (projPotFFT == NULL)
	{
		PrintErrorMsg("AtomProjPotFFT", VTEMLAB_ENULLPTR);
		return VTEMLAB_ENULLPTR;
	}

	double fxStart = -1.0 / (2 * Lx / Nx);
	double fyStart = -1.0 / (2 * Ly / Ny);
	double dfx = 1.0 / Lx;
	double dfy = 1.0 / Ly;

	double scaleCoeff = 2.0 * MathPI * BohrR * EleCharge;
	for (int yIdx = 0; yIdx < Ny; yIdx++)
	{
		double fy = fyStart + (double)yIdx * dfy;
		for (int xIdx = 0; xIdx < Nx; xIdx++)
		{
			double fx = fxStart + (double)xIdx * dfx;
			double q = sqrt(fx * fx + fy * fy);
			int idx = yIdx * Nx + xIdx;
			double sf = AtomScatteringFactor(eleScattParams, q);
			projPotFFT[idx].real = scaleCoeff * sf;
			projPotFFT[idx].imag = 0.0;
		}
	}

	return VTEMLAB_SUCCESS;
}


// Calculate the projected potential using convolution (atoms of single type):
// SHIFTED_VEC: shifted ProjPot;
// UNSHIFTED_VEC: unshifted ProjPot;
int MonoAtomProjPot_conv(MKL_Complex16* projPot, MKL_Complex16* convKer, 
	double eleScattParams[4][3], TypeAtomList* destAtomList, int Nx, int Ny, 
	int cellNum[2], double lattConsts[2], int mode)
{
	int vecLen;
	vecLen = Ny * Nx;
	// initializing the convolution kernel:
	int i;
	for (i = 0; i < vecLen; i++)
	{
		projPot[i] = { 0.0,0.0 };
		convKer[i] = { 0.0,0.0 };
	}

	double Lx, Ly, x, y, dx, dy, fx, fy, fxStart, fyStart, dfx, dfy, 
		xshift, yshift;
	Lx = (double)cellNum[0] * lattConsts[0];
	Ly = (double)cellNum[1] * lattConsts[1];
	dx = Lx / (double)Nx;
	dy = Ly / (double)Ny;

	xshift = Lx / 2.0;
	yshift = Ly / 2.0;

	fxStart = -1 / (2.0 * dx);
	fyStart = -1 / (2.0 * dy);
	dfx = 1.0 / Lx;
	dfy = 1.0 / Ly;

	int errorCode;
	errorCode = AtomProjPotFFT(projPot, eleScattParams, Lx, Ly, Nx, Ny);
	if (errorCode)
	{
		PrintErrorMsg("MonoAtomProjPot_conv", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}
	fftshift_Z_2D_IP(projPot, Nx, Ny);

	// Constructing fft descriptor:
	DFTI_DESCRIPTOR_HANDLE fftDescHandle;
	MKL_LONG fftStat, matSize[2] = { (MKL_LONG)Ny,(MKL_LONG)Nx };
	fftStat = DftiCreateDescriptor(&fftDescHandle, DFTI_DOUBLE, DFTI_COMPLEX, 
		2, matSize);
	fftStat = DftiCommitDescriptor(fftDescHandle);

	// Constructing the convolution kernel:
	TypeAtomList* tmpAtomNode = destAtomList->NextAtom;
	int ixCell, iyCell, ix, iy;
	double KerFac;
	while (tmpAtomNode != NULL)
	{
		for (ixCell = 0; ixCell < cellNum[0]; ixCell++)
		{
			x = tmpAtomNode->CoordX + (double)ixCell * lattConsts[0] - xshift;
			for (iyCell = 0; iyCell < cellNum[1]; iyCell++)
			{
				y = tmpAtomNode->CoordY + (double)iyCell * lattConsts[1] - yshift;
				for (iy = 0; iy < Ny; iy++)
				{
					fy = fyStart + (double)iy * dfy;
					for (ix = 0; ix < Nx; ix++)
					{
						i = iy * Nx + ix;
						fx = fxStart + (double)ix * dfx;
						KerFac = -2.0 * MathPI * (fx * x + fy * y);
						convKer[i].real += cos(KerFac) * tmpAtomNode->EleProp;
						convKer[i].imag += sin(KerFac) * tmpAtomNode->EleProp;
					}
				}
			}
		}
		tmpAtomNode = tmpAtomNode->NextAtom;
	}
	fftshift_Z_2D_IP(convKer, Nx, Ny);
	vzMul(vecLen, convKer, projPot, projPot);

	fftStat = DftiComputeBackward(fftDescHandle, projPot);

	// Scale ProjPot:
	double fftScaleCoeff = 1.0 / (Lx * Ly);
	cblas_zdscal(vecLen, fftScaleCoeff, projPot, 1);

	if (mode == UNFFTSHIFTED_VEC)
		ifftshift_Z_2D_IP(projPot, Nx, Ny);
	fftStat = DftiFreeDescriptor(&fftDescHandle);

	return VTEMLAB_SUCCESS;
}

// Calculate the projected potential using convolution (atoms of multiple types):
int MultiAtomProjPot_conv(double* scattParams, double*& projPot, 
	SliceTypeList* destTypeList, int Nx, int Ny, int cellNum[2], 
	double lattConsts[2])
{
	int vecLen = Ny * Nx;
	projPot = (double*)mkl_malloc(vecLen * sizeof(double), 64);
	if (projPot == NULL)
	{
		PrintErrorMsg("MultiAtomProjPot_conv", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}

	SliceTypeList* tmpTypeNode = destTypeList->NextType;
	int i, atomType, errorCode;
	double eleScattParams[4][3];
	MKL_Complex16* tmpProjPot = 
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (tmpProjPot == NULL)
	{
		PrintErrorMsg("MultiAtomProjPot_conv", VTEMLAB_ENOMEM);
		mkl_free(projPot);
		return VTEMLAB_ENOMEM;
	}

	MKL_Complex16* convKer = 
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (convKer == NULL)
	{
		PrintErrorMsg("MultiAtomProjPot_conv", VTEMLAB_ENOMEM);
		mkl_free(projPot);
		mkl_free(tmpProjPot);
		return VTEMLAB_ENOMEM;
	}

	// initializing
	for (i = 0; i < vecLen; i++)
	{
		projPot[i] = 0.0;
		tmpProjPot[i] = { 0.0,0.0 };
		convKer[i] = { 0.0,0.0 };
	}

	while (tmpTypeNode != NULL)
	{
		atomType = tmpTypeNode->AtomType;
		errorCode = ExtractEleScattParam(atomType, scattParams, eleScattParams);
		if (errorCode)
		{
			PrintErrorMsg("MultiAtomProjPot_conv", VTEMLAB_FAILURE);
			mkl_free(projPot);
			mkl_free(tmpProjPot);
			mkl_free(convKer);
			return VTEMLAB_FAILURE;
		}

		errorCode = MonoAtomProjPot_conv(tmpProjPot, convKer, eleScattParams, tmpTypeNode->AtomListHead, Nx, Ny, cellNum, lattConsts, FFTSHIFTED_VEC);
		if (errorCode)
		{
			PrintErrorMsg("MultiAtomProjPot_conv", VTEMLAB_FAILURE);
			mkl_free(projPot);
			mkl_free(tmpProjPot);
			mkl_free(convKer);
			return VTEMLAB_FAILURE;
		}
		for (i = 0; i < vecLen; i++)
			projPot[i] += tmpProjPot[i].real;

		tmpTypeNode = tmpTypeNode->NextType;
	}
	ifftshift_D_2D_IP(projPot, Nx, Ny);
	mkl_free(tmpProjPot);
	mkl_free(convKer);

	return VTEMLAB_SUCCESS;
}

// Save projected potential in a special format:
// first line in the file: Lx (double), Ly (double), Nx (int), Ny (int), SliceDist (double)
// second line in the file: projected potential values as a double-precision real vector.
int SaveProjPot(char* filename, double* projPot, double Lx, double Ly, 
	int Nx, int Ny, double sliceDist)
{
	if ((Nx == 0) || (Ny == 0))
	{
		PrintErrorMsg("SaveProjPot", VTEMLAB_EINVAL);
		return VTEMLAB_EINVAL;
	}

	int ix, iy;
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "w+");
	if (err)
	{
		PrintErrorMsg("SaveProjPot", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}
	fprintf_s(destFile, "%lf\t%lf\t%d\t%d\t%lf\t", Lx, Ly, Nx, Ny, sliceDist);
	for (ix = 5; ix < Nx; ix++)
	{
		if (ix != Nx - 1)
			fprintf_s(destFile, "%d\t", 0);
		else
			fprintf_s(destFile, "%d\n", 0);
	}

	for (iy = 0; iy < Ny; iy++)
	{
		for (ix = 0; ix < Nx - 1; ix++)
		{
			fprintf_s(destFile, "%lf\t", projPot[iy * Nx + ix]);
		}
		if (iy != Ny - 1)
			fprintf_s(destFile, "%lf\n", projPot[iy * Nx + ix]);
		else
			fprintf_s(destFile, "%lf", projPot[iy * Nx + ix]);
	}

	fclose(destFile);

	return VTEMLAB_SUCCESS;
}

// Load projected potential in a special format:
// first line in the file: Lx (double), Ly (double), Nx (int), Ny (int), SliceDist (double)
// second line in the file: projected potential values as a double-precision real vector.
// Note that argument ProjPot needs to be initialized before being passed to the function:
int LoadProjPot(char* filename, double*& projPot, double& Lx, double& Ly, 
	int& Nx, int& Ny, double& sliceDist)
{
	int i, vecLen, itemp;
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "r");
	if (err)
	{
		PrintErrorMsg("LoadProjPot", VTEMLAB_EFILER);
		return VTEMLAB_EFILER;
	}
	fscanf_s(destFile, "%lf\t%lf\t%d\t%d\t%lf", &Lx, &Ly, &Nx, &Ny, &sliceDist);
	for (i = 5; i < Nx; i++)
	{
		fscanf_s(destFile, "%d", &itemp);
	}
	vecLen = Ny * Nx;
	projPot = (double*)mkl_malloc(vecLen * sizeof(double), 64);
	if (projPot == NULL)
	{
		PrintErrorMsg("LoadProjPot", VTEMLAB_ENOMEM);
		fclose(destFile);
		return VTEMLAB_ENOMEM;
	}
	for (i = 0; i < vecLen; i++)
		fscanf_s(destFile, "%lf", &projPot[i]);
	fclose(destFile);

	return VTEMLAB_SUCCESS;
}


// convert text-image-type projected potential file to binary image:
int ConvertProjPotFromTextToBinary(char* textFilename, char* binFilename)
{
	int Nx, Ny;
	double Lx, Ly, sliceDist;
	double* projPot;
	int errorCode = LoadProjPot(textFilename, projPot, Lx, Ly, Nx, Ny, sliceDist);
	if (errorCode)
	{
		PrintErrorMsg("ConvertProjPotFromTextToBinary", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	errorCode = SaveBinaryProjPot(binFilename, projPot, Lx, Ly, Nx, Ny, sliceDist);
	if (errorCode)
	{
		PrintErrorMsg("ConvertProjPotFromTextToBinary", VTEMLAB_FAILURE);
		mkl_free(projPot);
		return VTEMLAB_FAILURE;
	}

	mkl_free(projPot);

	return VTEMLAB_SUCCESS;
}


// convert binary projected potential files to text files.
int ConvertProjPotFromBinaryToText(char* binFilename, char* textFilename)
{
	int Nx, Ny;
	double Lx, Ly, sliceDist;
	double* projPot;
	int opStat = LoadBinaryProjPot(binFilename, projPot, Lx, Ly, Nx, Ny, sliceDist);
	if (opStat)
	{
		PrintErrorMsg("ConvertProjPotFromBinaryToText", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	opStat = SaveProjPot(textFilename, projPot, Lx, Ly, Nx, Ny, sliceDist);
	if (opStat)
	{
		PrintErrorMsg("ConvertProjPotFromBinaryToText", VTEMLAB_FAILURE);
		mkl_free(projPot);
		return VTEMLAB_FAILURE;
	}

	mkl_free(projPot);

	return VTEMLAB_SUCCESS;
}


// save binary-type projected potential file
int SaveBinaryProjPot(char* binFilename, double* projPot, double Lx, double Ly,
	int Nx, int Ny, double sliceDist)
{
	FILE* destFile;
	if (fopen_s(&destFile, binFilename, "wb") == 0)
	{
		int numWritten = 0;
		numWritten += fwrite(&Lx, sizeof(double), 1, destFile);
		numWritten += fwrite(&Ly, sizeof(double), 1, destFile);
		numWritten += fwrite(&Nx, sizeof(int), 1, destFile);
		numWritten += fwrite(&Ny, sizeof(int), 1, destFile);
		numWritten += fwrite(&sliceDist, sizeof(double), 1, destFile);
		numWritten += fwrite(projPot, sizeof(double), Ny * Nx, destFile);
		printf("Binary file save as\n%s\nWrote %d items.\n\n", 
			binFilename, numWritten);
		fclose(destFile);
	}
	else
	{
		PrintErrorMsg("SaveBinaryProjPot", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}

	return VTEMLAB_SUCCESS;
}


// load binary-type projected potential file
int LoadBinaryProjPot(char* binFilename, double*& projPot, double& Lx, double& Ly,
	int& Nx, int& Ny, double& sliceDist)
{
	FILE* destFile;
	if (fopen_s(&destFile, binFilename, "rb") == 0)
	{
		int numRead = 0;
		numRead += fread(&Lx, sizeof(double), 1, destFile);
		numRead += fread(&Ly, sizeof(double), 1, destFile);
		numRead += fread(&Nx, sizeof(int), 1, destFile);
		numRead += fread(&Ny, sizeof(int), 1, destFile);
		numRead += fread(&sliceDist, sizeof(double), 1, destFile);

		projPot = (double*)mkl_malloc(Ny * Nx * sizeof(double), 64);
		if (projPot == NULL)
		{
			PrintErrorMsg("LoadBinaryProjPot", VTEMLAB_ENOMEM);
			fclose(destFile);
			return VTEMLAB_ENOMEM;
		}

		numRead += fread(projPot, sizeof(double), Ny * Nx, destFile);
		printf("Read %d items.\n\n", numRead);
		fclose(destFile);
	}
	else
	{
		PrintErrorMsg("LoadBinaryProjPot", VTEMLAB_EFILER);
		return VTEMLAB_EFILER;
	}

	return VTEMLAB_SUCCESS;
}


// Collectively read in the crystal files and write encrypted projected potential files 
// under the same directory (convolutional way):
int CrysFileToProjPot_conv(char* folder, int expanNum[2], int Nx, int Ny,
	int fileCoordType)
{
	int sliceNum = 0;
	SliceList* slices;
	int errorCode = CreateSliceList(folder, slices, sliceNum);
	if (errorCode)
	{
		PrintErrorMsg("CrysFileToProjPot_conv", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	double* scattParams;
	errorCode = LoadScattFacTable(scattParams);
	if (errorCode)
	{
		PrintErrorMsg("CrysFileToProjPot_conv", VTEMLAB_FAILURE);
		FreeSliceList(slices);
		return VTEMLAB_FAILURE;
	}

	SliceList* tmpSlice = slices->NextSlice;
	double lattConsts[2] = { tmpSlice->LattConstA,tmpSlice->LattConstB };
	double Lx = (double)expanNum[0] * lattConsts[0];
	double Ly = (double)expanNum[1] * lattConsts[1];

	int sliceCount = 0;
	while (tmpSlice != NULL)
	{
		sliceCount++;
		SliceTypeList* tmpTypeList;
		if (fileCoordType == FRAC_COORD)
			errorCode = SliceNodeToTypeList_Frac(tmpSlice, tmpTypeList, 1.0e-4, 1.0e-4);
		else if (fileCoordType == CART_COORD)
			errorCode = SliceNodeToTypeList_Cart(tmpSlice, tmpTypeList, 1.0e-4, 1.0e-4);
		else
		{
			PrintErrorMsg("CrysFileToProjPot_conv", VTEMLAB_EINVAL);
			FreeSliceList(slices);
			mkl_free(scattParams);
			return VTEMLAB_EINVAL;
		}

		if (errorCode)
		{
			PrintErrorMsg("CrysFileToProjPot_conv", VTEMLAB_FAILURE);
			FreeSliceList(slices);
			mkl_free(scattParams);
			return VTEMLAB_FAILURE;
		}

		double* projPot;
		errorCode = MultiAtomProjPot_conv(scattParams, projPot, tmpTypeList, 
			Nx, Ny, expanNum, lattConsts);
		if (errorCode)
		{
			PrintErrorMsg("CrysFileToProjPot_conv", VTEMLAB_FAILURE);
			FreeSliceList(slices);
			mkl_free(scattParams);
			FreeSliceTypeList(tmpTypeList);
			return VTEMLAB_FAILURE;
		}

		char filename[FILENAME_MAX];
		snprintf(filename, FILENAME_MAX, "%s\\p%d.txt", folder, sliceCount);

		errorCode = SaveProjPot(filename, projPot, Lx, Ly, Nx, Ny, tmpSlice->SliceDist);
		if (errorCode)
		{
			PrintErrorMsg("CrysFileToProjPot_conv", VTEMLAB_FAILURE);
			FreeSliceList(slices);
			mkl_free(scattParams);
			FreeSliceTypeList(tmpTypeList);
			mkl_free(projPot);
			return VTEMLAB_FAILURE;
		}

		FreeSliceTypeList(tmpTypeList);
		mkl_free(projPot);

		tmpSlice = tmpSlice->NextSlice;
	}

	FreeSliceList(slices);
	mkl_free(scattParams);

	return VTEMLAB_SUCCESS;
}


int CrysFileToBinProjPot_conv(char* folder, int expanNum[2], int Nx, int Ny,
	int fileCoordType)
{
	int sliceNum = 0;
	SliceList* slices;
	int errorCode = CreateSliceList(folder, slices, sliceNum);
	if (errorCode)
	{
		PrintErrorMsg("CrysFileToBinProjPot_conv", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	double* scattParams;
	errorCode = LoadScattFacTable(scattParams);
	if (errorCode)
	{
		PrintErrorMsg("CrysFileToBinProjPot_conv", VTEMLAB_FAILURE);
		FreeSliceList(slices);
		return VTEMLAB_FAILURE;
	}

	SliceList* tmpSlice = slices->NextSlice;
	double lattConsts[2] = { tmpSlice->LattConstA,tmpSlice->LattConstB };
	double Lx = (double)expanNum[0] * lattConsts[0];
	double Ly = (double)expanNum[1] * lattConsts[1];

	int sliceCount = 0;
	while (tmpSlice != NULL)
	{
		sliceCount++;
		SliceTypeList* tmpTypeList;
		if (fileCoordType == FRAC_COORD)
			errorCode = SliceNodeToTypeList_Frac(tmpSlice, tmpTypeList, 1.0e-4, 1.0e-4);
		else if (fileCoordType == CART_COORD)
			errorCode = SliceNodeToTypeList_Cart(tmpSlice, tmpTypeList, 1.0e-4, 1.0e-4);
		else
		{
			PrintErrorMsg("CrysFileToBinProjPot_conv", VTEMLAB_EINVAL);
			FreeSliceList(slices);
			mkl_free(scattParams);
			return VTEMLAB_EINVAL;
		}

		if (errorCode)
		{
			PrintErrorMsg("CrysFileToBinProjPot_conv", VTEMLAB_FAILURE);
			FreeSliceList(slices);
			mkl_free(scattParams);
			return VTEMLAB_FAILURE;
		}

		double* projPot;
		errorCode = MultiAtomProjPot_conv(scattParams, projPot, tmpTypeList, Nx, Ny,
			expanNum, lattConsts);
		if (errorCode)
		{
			PrintErrorMsg("CrysFileToBinProjPot_conv", VTEMLAB_FAILURE);
			FreeSliceList(slices);
			mkl_free(scattParams);
			FreeSliceTypeList(tmpTypeList);
			return VTEMLAB_FAILURE;
		}

		char filename[FILENAME_MAX];
		snprintf(filename, FILENAME_MAX, "%s\\p%d.bin", folder, sliceCount);

		errorCode = SaveBinaryProjPot(filename, projPot, Lx, Ly, Nx, Ny, tmpSlice->SliceDist);
		if (errorCode)
		{
			PrintErrorMsg("CrysFileToBinProjPot_conv", VTEMLAB_FAILURE);
			FreeSliceList(slices);
			mkl_free(scattParams);
			FreeSliceTypeList(tmpTypeList);
			mkl_free(projPot);
			return VTEMLAB_FAILURE;
		}

		FreeSliceTypeList(tmpTypeList);
		mkl_free(projPot);

		tmpSlice = tmpSlice->NextSlice;
	}

	FreeSliceList(slices);
	mkl_free(scattParams);

	return VTEMLAB_SUCCESS;
}


// Decrypt the encrypted projected potential files
int DecryptProjPotFiles(char* folder, int sliceNum)
{
	double* projPot, Lx, Ly, sliceDist;
	int Nx, Ny, errorCode;
	for (int sliceIdx = 0; sliceIdx < sliceNum; sliceIdx++)
	{
		char filename[FILENAME_MAX];
		snprintf(filename, FILENAME_MAX, "%s\\p%d.txt", folder, sliceIdx + 1);

		errorCode = LoadProjPot(filename, projPot, Lx, Ly, Nx, Ny, sliceDist);
		if (errorCode)
		{
			PrintErrorMsg("DecryptProjPotFiles", VTEMLAB_FAILURE);
			return VTEMLAB_FAILURE;
		}
		
		memset(filename, 0, FILENAME_MAX);
		snprintf(filename, FILENAME_MAX, "%s\\P_dec%d.txt", 
			folder, sliceIdx + 1);

		errorCode = SaveDvecAsMat(Nx, Ny, projPot, filename);
		if (errorCode)
		{
			PrintErrorMsg("DecryptProjPotFiles", VTEMLAB_FAILURE);
			mkl_free(projPot);
			return VTEMLAB_FAILURE;
		}

		mkl_free(projPot);
	}
	
	return VTEMLAB_SUCCESS;
}