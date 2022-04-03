/*
	thermo.cpp is the open source part of vtemlab v0.0 engine,
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

#include <math.h>
#include <stdio.h>

#include "mathconst.h"
#include "mode.h"
#include "potential.h"
#include "slice.h"
#include "status.h"
#include "thermo.h"
#include "vtemlabio.h"
#include "vtemlabmath.h"

// Calculate the projected potential for one type of atoms using convolution (atoms of 
// single type): for the case of thermo effect.
// mode == SHIFTED_VEC: shifted ProjPot;
// mode == UNSHIFTED_VEC: unshifted ProjPot;
int MonoAtomProjPot_thermo_conv(MKL_Complex16* projPot, MKL_Complex16* convKer,
	double eleScattParams[4][3], double Lx, double Ly, int Nx, int Ny, int mode)
{
	int vecLen, i, errorCode;
	vecLen = Nx * Ny;
	// initializing ProjPot:
	for (i = 0; i < vecLen; i++)
		projPot[i] = { 0.0,0.0 };

	// initializing convolution seed:
	errorCode = AtomProjPotFFT(projPot, eleScattParams, Lx, Ly, Nx, Ny);
	if (errorCode)
	{
		PrintErrorMsg("MonoAtomProjPot_thermo_conv", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}
	fftshift_Z_2D_IP(projPot, Nx, Ny);

	// Constructing fft descriptor:
	DFTI_DESCRIPTOR_HANDLE fftDescHandle;
	MKL_LONG fftStat, matSize[2] = { (MKL_LONG)Ny,(MKL_LONG)Nx };
	fftStat = DftiCreateDescriptor(&fftDescHandle, DFTI_DOUBLE, DFTI_COMPLEX, 
		2, matSize);
	fftStat = DftiCommitDescriptor(fftDescHandle);

	fftshift_Z_2D_IP(convKer, Nx, Ny);
	vzMul(vecLen, convKer, projPot, projPot);
	fftStat = DftiComputeBackward(fftDescHandle, projPot);

	// Scale ProjPot:
	double fftScaleCoeff = 1 / (Lx * Ly);
	cblas_zdscal(vecLen, fftScaleCoeff, projPot, 1);

	if (mode == UNFFTSHIFTED_VEC)
		ifftshift_Z_2D_IP(projPot, Nx, Ny);
	fftStat = DftiFreeDescriptor(&fftDescHandle);

	return VTEMLAB_SUCCESS;
}


// Temporary: slice proj pot with thermal vibration (Einstein model)
int SliceComplexProjPot_thermo_conv(ThermoSliceList* thermoSlice,
	int configIdx, int stackIdx, double* scattParams, MKL_Complex16* fxMesh,
	MKL_Complex16* fyMesh, MKL_Complex16* fMesh, MKL_Complex16* convKer,
	MKL_Complex16* projPot, MKL_Complex16* container, int Nx, int Ny,
	int projPotOutMode)
{
	bool anyNullPtr = (thermoSlice == NULL) ||
		(thermoSlice->typeListHead == NULL) ||
		(thermoSlice->posiOrientDispArray == NULL) ||
		(thermoSlice->rndThermoDispX == NULL) ||
		(thermoSlice->rndThermoDispY == NULL) ||
		(scattParams == NULL) || (convKer == NULL) ||
		(projPot == NULL) || (container == NULL);
	if (anyNullPtr)
	{
		PrintErrorMsg("SliceComplexProjPot_thermo_conv", VTEMLAB_ENULLPTR);
		return VTEMLAB_ENULLPTR;
	}

	int vecLen = Ny * Nx;
	int expanNum[2] = { thermoSlice->sliceExpanNum[0],thermoSlice->sliceExpanNum[1] };
	double lattConst[2] = { thermoSlice->lattConstA,thermoSlice->lattConstB };
	double Lx = (double)(expanNum[0]) * lattConst[0];
	double Ly = (double)(expanNum[1]) * lattConst[1];

	double xshift = Lx / 2.0;
	double yshift = Ly / 2.0;

	double dx = Lx / (double)Nx;
	double dy = Ly / (double)Ny;

	int rndNumPerStack = expanNum[1] * expanNum[0];
	int rndNumPerConfig = thermoSlice->stackNum * rndNumPerStack;

	// initialize projPot
	for (int i = 0; i < vecLen; i++)
		projPot[i] = { 0.0,0.0 };

	SliceTypeList* tmpTypeNode = thermoSlice->typeListHead->NextType;
	double kerFacScal = -2.0 * MathPI;
	int errorCode;
	while (tmpTypeNode != NULL)
	{
		int tmpType = tmpTypeNode->AtomType;

		// extract elemental scattering factors:
		double eleScattParams[4][3];
		errorCode = ExtractEleScattParam(tmpType, scattParams, eleScattParams);
		if (errorCode)
		{
			PrintErrorMsg("SliceComplexProjPot_thermo_conv", VTEMLAB_FAILURE);
			return VTEMLAB_FAILURE;
		}

		// initialize convKer
		for (int i = 0; i < vecLen; i++)
			convKer[i] = { 0.0,0.0 };

		// build convKer
		TypeAtomList* tmpAtomNode = tmpTypeNode->AtomListHead->NextAtom;
		while (tmpAtomNode != NULL)
		{
			int tmpPosiIdx = tmpAtomNode->PositionIdx;
			double tmpEleProp = tmpAtomNode->EleProp;

			for (int yExpanIdx = 0; yExpanIdx < expanNum[1]; yExpanIdx++)
			{
				double y0 = tmpAtomNode->CoordY +
					(double)yExpanIdx * lattConst[1] - yshift;

				for (int xExpanIdx = 0; xExpanIdx < expanNum[0]; xExpanIdx++)
				{
					double x0 = tmpAtomNode->CoordX +
						(double)xExpanIdx * lattConst[0] - xshift;

					int rndIdx = configIdx * rndNumPerConfig + stackIdx * rndNumPerStack +
						yExpanIdx * expanNum[0] + xExpanIdx;

					double x = x0 + thermoSlice->rndThermoDispX[tmpPosiIdx][rndIdx];

					if (x > xshift)
						x -= Lx;
					else if (x < -xshift)
						x += Lx;

					double y = y0 + 
						thermoSlice->rndThermoDispY[tmpPosiIdx][rndIdx];
					if (y > yshift)
						y -= Ly;
					else if (y < -yshift)
						y += Ly;

					// write the current atom into convKer:
					MKL_Complex16 complexY[1] = { y,0.0 };

					cblas_zcopy(vecLen, fxMesh, 1, fMesh, 1);
					cblas_zdscal(vecLen, x, fMesh, 1);
					cblas_zaxpy(vecLen, complexY, fyMesh, 1, fMesh, 1);
					cblas_zdscal(vecLen, kerFacScal, fMesh, 1);
					vzExp(vecLen, fMesh, fMesh);
					cblas_zdscal(vecLen, tmpEleProp, fMesh, 1);
					vzAdd(vecLen, convKer, fMesh, convKer);
				}
			}

			tmpAtomNode = tmpAtomNode->NextAtom;
		}
		int opStat = MonoAtomProjPot_thermo_conv(container, convKer, eleScattParams, Lx, Ly,
			Nx, Ny, FFTSHIFTED_VEC);
		if (opStat)
		{
			PrintErrorMsg("SliceComplexProjPot_thermo_conv", VTEMLAB_FAILURE);
			return VTEMLAB_FAILURE;
		}

		for (int idx = 0; idx < vecLen; idx++)
			projPot[idx].imag += container[idx].real;

		tmpTypeNode = tmpTypeNode->NextType;
	}

	if (projPotOutMode == UNFFTSHIFTED_VEC)
		ifftshift_Z_2D_IP(projPot, Nx, Ny);

	return VTEMLAB_SUCCESS;
}


// Bandwidth limit a transmission function:
void BwlSliceTF(MKL_Complex16* sliceTF, double Lx, double Ly, int Nx, int Ny,
	double bwlProp, int tfInMode, int tfOutMode)
{
	double dx = Lx / (double)Nx;
	double dy = Ly / (double)Ny;

	double fxStart = -1.0 / (2.0 * dx);
	double fyStart = -1.0 / (2.0 * dy);

	double dfx = 1.0 / Lx;
	double dfy = 1.0 / Ly;

	double minR = (dx > dy) ? (1.0 / (2.0 * dx)) : (1.0 / (2.0 * dy));
	double apertRadiusSqr = bwlProp * minR;
	apertRadiusSqr = apertRadiusSqr * apertRadiusSqr;

	// Construct a fft handle
	MKL_LONG fftStat, matSize[2] = { (MKL_LONG)Ny,(MKL_LONG)Nx };
	DFTI_DESCRIPTOR_HANDLE descHandle_1;
	fftStat = DftiCreateDescriptor(&descHandle_1, DFTI_DOUBLE, DFTI_COMPLEX, 2, matSize);
	fftStat = DftiCommitDescriptor(descHandle_1);

	if (tfInMode == UNFFTSHIFTED_VEC)
		fftshift_Z_2D_IP(sliceTF, Nx, Ny);

	fftStat = DftiComputeForward(descHandle_1, sliceTF);

	// Add aperture
	ifftshift_Z_2D_IP(sliceTF, Nx, Ny);

	for (int iy = 0; iy < Ny; iy++)
	{
		double fy = fyStart + (double)iy * dfy;
		for (int ix = 0; ix < Nx; ix++)
		{
			double fx = fxStart + (double)ix * dfx;
			if ((fx * fx + fy * fy) > apertRadiusSqr)
				sliceTF[iy * Nx + ix] = { 0.0,0.0 };
		}
	}

	fftshift_Z_2D_IP(sliceTF, Nx, Ny);

	fftStat = DftiComputeBackward(descHandle_1, sliceTF);
	fftStat = DftiFreeDescriptor(&descHandle_1);

	int vecLen = Ny * Nx;
	for (int idx = 0; idx < vecLen; idx++)
	{
		double normCoeff = sqrt(sliceTF[idx].real * sliceTF[idx].real +
			sliceTF[idx].imag * sliceTF[idx].imag);
		sliceTF[idx].real /= normCoeff;
		sliceTF[idx].imag /= normCoeff;
	}

	if (tfOutMode == UNFFTSHIFTED_VEC)
		ifftshift_Z_2D_IP(sliceTF, Nx, Ny);
}


// Scan the ThermoSliceList
void ScanThermoSliceList(ThermoSliceList* tsList)
{
	ThermoSliceList* tsNode = tsList->nextSlice;
	while (tsNode != NULL)
	{
		printf(">> Slice %d:\n\n", tsNode->inStackSliceIdx + 1);

		printf("> Lattice constants A and B:\n%.5f\t%.5f\n\n",
			tsNode->lattConstA, tsNode->lattConstB);
		printf("> Slice distance: %.5f\n\n", tsNode->sliceDist);
		printf("> Lattice expansion along x and y:\n%d\t%d\n\n",
			tsNode->sliceExpanNum[0], tsNode->sliceExpanNum[1]);
		printf("> Config num: %d\n> Stack num: %d\n\n",
			tsNode->configNum, tsNode->stackNum);
		printf("> %d independent positions\n\n", tsNode->posiNum);

		int rndNumPerPosi = 0;
		rndNumPerPosi = tsNode->sliceExpanNum[0] *
			tsNode->sliceExpanNum[1] * tsNode->configNum *
			tsNode->stackNum;

		printf("> %d random displacements per position per dimension\n\n", rndNumPerPosi);

		// check if typeListHead is attached
		if (tsNode->typeListHead == NULL)
			printf("> WARNING: No slice type list is attached!\n\n");
		else
		{
			printf("> Slice type list attached.\n\n");
			ScanSliceTypeList(tsNode->typeListHead);
		}

		// check if posiOrientDispArray is initialized
		if (tsNode->posiOrientDispArray == NULL)
			printf("> NULL ARRAY: posiOrientDispArray uninitialized!\n\n");
		else
			printf("> posiOrientDispArray initialized.\n\n");

		// check if rndThermoDispX is initialized
		if (tsNode->rndThermoDispX == NULL)
			printf("> NULL ARRAY: rndThermoDispX uninitialized!\n\n");
		else
			printf("> rndThermoDispX initialized.\n\n");

		// check if rndThermoDispY is initialized
		if (tsNode->rndThermoDispY == NULL)
			printf("> NULL ARRAY: rndThermoDispY uninitialized!\n\n");
		else
			printf("> rndThermoDispY initialized.\n\n");

		tsNode = tsNode->nextSlice;
	}
}


// Free the ThermoSliceList
void FreeThermoSliceList(ThermoSliceList*& tsList)
{
	ThermoSliceList* prevNode = tsList;
	ThermoSliceList* postNode = tsList->nextSlice;
	while (postNode != NULL)
	{
		free(prevNode);

		if (postNode->typeListHead != NULL)
			FreeSliceTypeList(postNode->typeListHead);

		mkl_free(postNode->posiOrientDispArray);
		for (int posiIdx = 0; posiIdx < postNode->posiNum; posiIdx++)
		{
			mkl_free(postNode->rndThermoDispX[posiIdx]);
			mkl_free(postNode->rndThermoDispY[posiIdx]);
		}
		mkl_free(postNode->rndThermoDispX);
		mkl_free(postNode->rndThermoDispY);

		prevNode = postNode;
		postNode = postNode->nextSlice;
	}
	free(prevNode);
}


// Calculate isotropic atomic vibration amplitude deviation
double IsoAtomMeanSqrVibAmp(double simuTemp, double bFactor_293K)
{
	double isoU = 0.0;
	isoU = bFactor_293K / (8 * MathPI * MathPI);
	isoU *= (simuTemp / 293);
	return isoU;
}


// Instead of calculating the deviation of thermal displacement with debye temperature,
// it is better to do it with the B factor (or Debye Waller factor), because acquiring an
// accurate debye temperature table is no easier than acquiring a B factor table.
// For more detailed information about B factors, see J.M. Zuo and J.C.H. Spence, Advanced
// Transmission Electron Microscopy.
void FillThermoInfoWithBFac(int typeNum, int* typeArray, double* bFactors_293K,
	SliceTypeList* destTypeList, double simuTemp, int positionNum,
	double* posiOrientThermoDispArray)
{
	// Initialize posiOrientThermoDispArray as all zeros
	for (int posiIdx = 0; posiIdx < positionNum; posiIdx++)
		posiOrientThermoDispArray[posiIdx] = 0.0;

	double tmpMassNum = 12.0, tmpBFactor = 0.0, tmpSqrThermoDisp = 0.0;
	int tmpType;
	SliceTypeList* tmpTypeNode = destTypeList->NextType;
	TypeAtomList* tmpTypeAtomNode;
	while (tmpTypeNode != NULL)
	{
		tmpType = tmpTypeNode->AtomType;

		for (int typeIdx = 0; typeIdx < typeNum; typeIdx++)
		{
			if (tmpType == typeArray[typeIdx])
			{
				tmpBFactor = bFactors_293K[typeIdx];
				tmpSqrThermoDisp = IsoAtomMeanSqrVibAmp(simuTemp, tmpBFactor);
				tmpTypeNode->MeanSquareDisplace = tmpSqrThermoDisp;
			}
		}

		tmpTypeAtomNode = tmpTypeNode->AtomListHead->NextAtom;
		while (tmpTypeAtomNode != NULL)
		{
			posiOrientThermoDispArray[tmpTypeAtomNode->PositionIdx] +=
				tmpTypeAtomNode->EleProp * tmpSqrThermoDisp;

			tmpTypeAtomNode = tmpTypeAtomNode->NextAtom;
		}

		tmpTypeNode = tmpTypeNode->NextType;
	}
	vdSqrt(positionNum, posiOrientThermoDispArray, posiOrientThermoDispArray);
}


// Create the random displacement deviation arrays for both x and y direction
int CreateRandomDispDeviXY(int expanNum[2], int configNum, int stackNum, int posiNum,
	double* posiOrientDispArray, VSLStreamStatePtr stream, double**& rndThermoDispX,
	double**& rndThermoDispY)
{
	if (posiOrientDispArray == NULL)
	{
		PrintErrorMsg("CreateRandomDispDeviXY", VTEMLAB_ENULLPTR);
		return VTEMLAB_ENULLPTR;
	}

	rndThermoDispX = (double**)mkl_malloc(posiNum * sizeof(double*), 64);
	if (rndThermoDispX == NULL)
	{
		PrintErrorMsg("CreateRandomDispDeviXY", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}

	rndThermoDispY = (double**)mkl_malloc(posiNum * sizeof(double*), 64);
	if (rndThermoDispY == NULL)
	{
		PrintErrorMsg("CreateRandomDispDeviXY", VTEMLAB_ENOMEM);
		free(rndThermoDispX);
		return VTEMLAB_ENOMEM;
	}

	int rndNumPerPosi = configNum * stackNum * expanNum[1] * expanNum[0];
	for (int posiIdx = 0; posiIdx < posiNum; posiIdx++)
	{
		rndThermoDispX[posiIdx] = (double*)mkl_malloc(rndNumPerPosi * sizeof(double), 64);
		if (rndThermoDispX[posiIdx] == NULL)
		{
			for (int freeIdx = 0; freeIdx < posiIdx; freeIdx++)
			{
				mkl_free(rndThermoDispX[freeIdx]);
				mkl_free(rndThermoDispY[freeIdx]);
			}
			mkl_free(rndThermoDispX);
			mkl_free(rndThermoDispY);

			PrintErrorMsg("CreateRandomDispDeviXY", VTEMLAB_ENOMEM);
			return VTEMLAB_ENOMEM;
		}

		rndThermoDispY[posiIdx] = (double*)mkl_malloc(rndNumPerPosi * sizeof(double), 64);
		if (rndThermoDispY[posiIdx] == NULL)
		{
			for (int freeIdx = 0; freeIdx < posiIdx; freeIdx++)
			{
				mkl_free(rndThermoDispX[freeIdx]);
				mkl_free(rndThermoDispY[freeIdx]);
			}
			mkl_free(rndThermoDispX[posiIdx]);
			mkl_free(rndThermoDispX);
			mkl_free(rndThermoDispY);

			PrintErrorMsg("CreateRandomDispDeviXY", VTEMLAB_ENOMEM);
			return VTEMLAB_ENOMEM;
		}

		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, rndNumPerPosi,
			rndThermoDispX[posiIdx], 0.0, posiOrientDispArray[posiIdx]);

		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, rndNumPerPosi,
			rndThermoDispY[posiIdx], 0.0, posiOrientDispArray[posiIdx]);
	}

	return VTEMLAB_SUCCESS;
}


// Create ThermoSliceList
int CreateThermoSliceList(SliceList* destSliceList, int coordType, int expanNum[2],
	ThermoSliceList*& tsList, int configNum, int stackNum, int typeNum,
	int* typeArray, double* bFactors_293K, double simuTemp)
{
	tsList = (ThermoSliceList*)malloc(sizeof(ThermoSliceList));
	ThermoSliceList* prevTsNode = tsList;
	ThermoSliceList* postTsNode;

	// initialize the random number generator
	VSLStreamStatePtr stream;
	vslNewStream(&stream, VSL_BRNG_MT19937, 777);

	int errorCode;
	int sliceCount = 0;
	SliceList* tmpSliceNode = destSliceList->NextSlice;
	while (tmpSliceNode != NULL)
	{
		postTsNode = (ThermoSliceList*)malloc(sizeof(ThermoSliceList));

		postTsNode->inStackSliceIdx = sliceCount;

		postTsNode->lattConstA = tmpSliceNode->LattConstA;
		postTsNode->lattConstB = tmpSliceNode->LattConstB;
		postTsNode->sliceDist = tmpSliceNode->SliceDist;

		postTsNode->sliceExpanNum[0] = expanNum[0];
		postTsNode->sliceExpanNum[1] = expanNum[1];

		postTsNode->configNum = configNum;
		postTsNode->stackNum = stackNum;

		postTsNode->posiNum = tmpSliceNode->PositionNum;

		// create inner SliceTypeList
		SliceTypeList* tmpTypeList;
		if (coordType == CART_COORD)
			errorCode = SliceNodeToTypeList_Cart(tmpSliceNode, tmpTypeList, 1.0e-5, 1.0e-5);
		else if (coordType == FRAC_COORD)
			errorCode = SliceNodeToTypeList_Frac(tmpSliceNode, tmpTypeList, 1.0e-5, 1.0e-5);
		else
		{
			PrintErrorMsg("CreateThermoSliceList", VTEMLAB_EINVAL);
			vslDeleteStream(&stream);
			FreeThermoSliceList(tsList);

			return VTEMLAB_EINVAL;
		}

		if (errorCode)
		{
			PrintErrorMsg("CreateThermoSliceList", VTEMLAB_FAILURE);
			vslDeleteStream(&stream);
			FreeThermoSliceList(tsList);

			return VTEMLAB_FAILURE;
		}

		postTsNode->typeListHead = tmpTypeList;

		// create posiOrientDispArray:
		postTsNode->posiOrientDispArray =
			(double*)mkl_malloc(postTsNode->posiNum * sizeof(double), 64);
		FillThermoInfoWithBFac(typeNum, typeArray, bFactors_293K,
			postTsNode->typeListHead, simuTemp, postTsNode->posiNum,
			postTsNode->posiOrientDispArray);

		// create rndThermoDispX and rndThermoDispY
		errorCode = CreateRandomDispDeviXY(expanNum, configNum, stackNum,
			postTsNode->posiNum, postTsNode->posiOrientDispArray,
			stream, postTsNode->rndThermoDispX, postTsNode->rndThermoDispY);
		if (errorCode)
		{
			PrintErrorMsg("CreateThermoSliceList", VTEMLAB_FAILURE);
			vslDeleteStream(&stream);
			FreeThermoSliceList(tsList);

			return VTEMLAB_FAILURE;
		}

		prevTsNode->nextSlice = postTsNode;
		prevTsNode = postTsNode;

		tmpSliceNode = tmpSliceNode->NextSlice;
		sliceCount++;
	}
	prevTsNode->nextSlice = NULL;

	vslDeleteStream(&stream);

	return VTEMLAB_SUCCESS;
}


// Prepare sample with thermal vibration
int CreateThermoSample(char* folder, int coordType, int expanNum[2], SliceList*& sList,
	ThermoSliceList*& tsList, int configNum, int stackNum, int& sliceNum, int typeNum,
	int* typeArray, double* bFactors_293K, double simuTemp)
{
	int errorCode = CreateSliceList(folder, sList, sliceNum);
	if (errorCode)
	{
		PrintErrorMsg("CreateThermoSample", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	errorCode = CreateThermoSliceList(sList, coordType, expanNum, tsList, configNum,
		stackNum, typeNum, typeArray, bFactors_293K, simuTemp);
	if (errorCode)
	{
		FreeSliceList(sList);
		PrintErrorMsg("CreateThermoSample", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	return VTEMLAB_SUCCESS;
}