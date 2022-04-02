/*
	multislice.cpp is the open source part of vtemlab v0.0 engine,
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

#include <math.h>
#include <stdio.h>

#include "mathconst.h"
#include "mode.h"
#include "multislice.h"
#include "optics.h"
#include "potential.h"
#include "slice.h"
#include "status.h"
#include "thermo.h"
#include "vtemlabio.h"
#include "vtemlabmath.h"

// Create PropKerList
int CreatePropKerList(SliceList* sList, PropKerList*& pkList,
	double wavLen, int expanNum[2], int Nx, int Ny, int pkMode)
{
	SliceList* tmpSNode = sList->NextSlice;
	PropKerList* prevPkNode, * postPkNode;
	pkList = (PropKerList*)malloc(sizeof(PropKerList));
	prevPkNode = pkList;

	int errorCode;
	int sliceIdx = 0;
	int vecLen = Ny * Nx;
	while (tmpSNode != NULL)
	{
		double Lx = (double)expanNum[0] * tmpSNode->LattConstA;
		double Ly = (double)expanNum[1] * tmpSNode->LattConstB;

		postPkNode = (PropKerList*)malloc(sizeof(PropKerList));

		postPkNode->inStackSliceIdx = sliceIdx;
		postPkNode->Lx = Lx;
		postPkNode->Ly = Ly;
		postPkNode->Nx = Nx;
		postPkNode->Ny = Ny;
		postPkNode->sliceDist = tmpSNode->SliceDist;
		postPkNode->propKer =
			(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
		if (postPkNode->propKer == NULL)
		{
			PrintErrorMsg("CreatePropKerList_SL", VTEMLAB_ENOMEM);
			// Left blank to add code freeing the PropKerList:
			FreePropKerList(pkList, false);

			return VTEMLAB_ENOMEM;
		}

		errorCode = FresPropKer(postPkNode->propKer, wavLen, tmpSNode->SliceDist,
			Lx, Ly, Nx, Ny);
		if (errorCode)
		{
			PrintErrorMsg("CreatePropKerList_SL", VTEMLAB_FAILURE);
			// Left blank to add code freeing the PropKerList:
			FreePropKerList(pkList, false);

			return VTEMLAB_FAILURE;
		}

		if (pkMode == FFTSHIFTED_VEC)
			fftshift_Z_2D_IP(postPkNode->propKer, Nx, Ny);

		prevPkNode->nextSlice = postPkNode;
		prevPkNode = postPkNode;

		tmpSNode = tmpSNode->NextSlice;
		sliceIdx++;
	}
	prevPkNode->nextSlice = NULL;

	return VTEMLAB_SUCCESS;
}


// Create PropKerList for slices with uniform spacing
int CreateUniformPropKerList(PropKerList*& pkList, int sliceNum, 
	double sliceSpacing, double wavLen, double Lx, double Ly, int Nx, int Ny, 
	int pkMode)
{
	pkList = (PropKerList*)malloc(sizeof(PropKerList));
	PropKerList* prevNode = pkList;
	PropKerList* postNode;

	int vecLen = Ny * Nx;
	MKL_Complex16* sharedPropKer =
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (sharedPropKer == NULL)
	{
		PrintErrorMsg("CreateUniformPropKerList", VTEMLAB_ENOMEM);
		FreePropKerList(pkList, true);

		return VTEMLAB_ENOMEM;
	}

	int errorCode = FresPropKer(sharedPropKer, wavLen, sliceSpacing, Lx, Ly, Nx, Ny);
	if (errorCode)
	{
		PrintErrorMsg("CreateUniformPropKerList", VTEMLAB_FAILURE);
		FreePropKerList(pkList, true);
		if (sharedPropKer != NULL)
			mkl_free(sharedPropKer);

		return VTEMLAB_FAILURE;
	}

	if (pkMode == FFTSHIFTED_VEC)
		fftshift_Z_2D_IP(sharedPropKer, Nx, Ny);

	for (int sliceIdx = 0; sliceIdx < sliceNum; sliceIdx++)
	{
		postNode = (PropKerList*)malloc(sizeof(PropKerList));

		postNode->inStackSliceIdx = sliceIdx;
		postNode->Lx = Lx;
		postNode->Ly = Ly;
		postNode->Nx = Nx;
		postNode->Ny = Ny;
		postNode->sliceDist = sliceSpacing;
		postNode->propKer = sharedPropKer;

		prevNode->nextSlice = postNode;
		prevNode = postNode;
	}

	prevNode->nextSlice = NULL;

	return VTEMLAB_SUCCESS;
}


// Automatically choose CreatePropKerList or CreateUniformPropKerList to create PropKerList
int AutoCreatePropKerList(SliceList* sList, PropKerList*& pkList, double wavLen,
	int expanNum[2], int Nx, int Ny, int pkMode, bool& sliceDistEqual)
{
	// check if sliceDist of each slice is identical (within an acceptant error)
	SliceList* tmpSlice = sList->NextSlice;
	double refSliceDist = tmpSlice->SliceDist;
	double Lx = (double)expanNum[0] * tmpSlice->LattConstA;
	double Ly = (double)expanNum[1] * tmpSlice->LattConstB;
	double acceptError = 1.0e-5;
	int sliceCount = 0;
	sliceDistEqual = true;
	while (tmpSlice != NULL)
	{
		sliceCount++;
		double sliceDistDiff = fabs(tmpSlice->SliceDist - refSliceDist);
		if (sliceDistDiff > acceptError)
			sliceDistEqual = false;

		tmpSlice = tmpSlice->NextSlice;
	}

	int errorCode;
	if (sliceDistEqual)
	{
		errorCode = CreateUniformPropKerList(pkList, sliceCount, refSliceDist, wavLen, Lx, Ly,
			Nx, Ny, pkMode);
	}
	else
	{
		errorCode = CreatePropKerList(sList, pkList, wavLen, expanNum, Nx, Ny,
			pkMode);
	}

	if (errorCode)
	{
		PrintErrorMsg("AutoCreatePropKerList_SL", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	return VTEMLAB_SUCCESS;
}


// Free the PropKerList
void FreePropKerList(PropKerList*& pkList, bool sliceDistEqual)
{
	PropKerList* prevNode = pkList;
	PropKerList* postNode = pkList->nextSlice;
	if (sliceDistEqual)
		mkl_free(postNode->propKer);
	while (postNode != NULL)
	{
		free(prevNode);
		if (!sliceDistEqual)
			mkl_free(postNode->propKer);

		prevNode = postNode;
		postNode = postNode->nextSlice;
	}
	free(prevNode);
}


// Calculate the wavelength (Angstrom) of electron beam given the accelerating voltage (KeV)
double EleWavLen(double KeV)
{
	return 12.3986 / sqrt((2.0 * 511.0 + KeV) * KeV);
}

// Calculate the interaction parameter sigma.
/*Note that the unit of electron wavelength is Angstrom, and the unit of sigma
is (Angstrom * kV)^-1, anyway, not all the units of voltage is kV.*/
double InteractionCoefficient(double KeV)
{
	double wavLen = EleWavLen(KeV);
	double sigma = 2.0 * MathPI / (wavLen * KeV) * 
		(511.0 + KeV) / (2.0 * 511.0 + KeV);
	return sigma;
}

// Generate the objective transfer function in reciprocal space, only 
// including Cs3, Cs5 and df
int ObjTransFunc(MKL_Complex16* otf, OtfParamSet param, double wavLen,
	double Lx, double Ly, int Nx, int Ny)
{
	if (otf == NULL)
	{
		PrintErrorMsg("ObjTransFunc", VTEMLAB_ENULLPTR);
		mkl_free(otf);
		return VTEMLAB_ENULLPTR;
	}
	double dx, dy, fxStart, fyStart, dfx, dfy, fx, fy, 
		angFreqPwr2, angFreqPwr4, angFreqPwr6, otfPhase;
	dx = Lx / (double)Nx;
	dy = Ly / (double)Ny;
	fxStart = -1.0 / (2.0 * dx);
	fyStart = -1.0 / (2.0 * dy);
	dfx = 1.0 / Lx;
	dfy = 1.0 / Ly;

	double df, Cs3, Cs5;
	df = param.defocus;
	Cs3 = param.Cs3 * 1.0e7;
	Cs5 = param.Cs5 * 1.0e7;

	int ix, iy;
	for (iy = 0; iy < Ny; iy++)
	{
		fy = fyStart + (double)iy * dfy;
		for (ix = 0; ix < Nx; ix++)
		{
			fx = fxStart + (double)ix * dfx;
			angFreqPwr2 = (fx * fx + fy * fy) * wavLen * wavLen;
			angFreqPwr4 = angFreqPwr2 * angFreqPwr2;
			angFreqPwr6 = angFreqPwr4 * angFreqPwr2;
			otfPhase = -2.0 * MathPI / wavLen * (0.5 * df * angFreqPwr2 +
				0.25 * Cs3 * angFreqPwr4 + 1.0 / 6.0 * Cs5 * angFreqPwr6);
			otf[iy * Nx + ix] = { cos(otfPhase),sin(otfPhase) };
		}
	}

	return 0;
}

// Add a circular aperture to OTF:
void AddCircApert(MKL_Complex16* OTF, double WavLen, double numApert, 
	double Lx, double Ly, int Nx, int Ny)
{
	double dx, dy, fxStart, fyStart, dfx, dfy, fx, fy, angFreqSqr, 
		scaledNumApertSqr;
	dx = Lx / (double)Nx;
	dy = Ly / (double)Ny;
	fxStart = -1.0 / (2.0 * dx);
	fyStart = -1.0 / (2.0 * dy);
	dfx = 1.0 / Lx;
	dfy = 1.0 / Ly;
	scaledNumApertSqr = 1.0e-3 * numApert;
	scaledNumApertSqr = scaledNumApertSqr * scaledNumApertSqr;

	int ix, iy;
	for (iy = 0; iy < Ny; iy++)
	{
		fy = fyStart + (double)iy * dfy;
		for (ix = 0; ix < Nx; ix++)
		{
			fx = fxStart + (double)ix * dfx;
			angFreqSqr = (fx * fx + fy * fy) * WavLen * WavLen;
			if (angFreqSqr > scaledNumApertSqr)
				OTF[iy * Nx + ix] = { 0.0,0.0 };
		}
	}
}


// Add an annular aperture to OTF:
void AddAnnularAperture(MKL_Complex16* otf, int otfMode, double wavLen, 
	double innerAngle, double outerAngle, double Lx, double Ly, int Nx, int Ny)
{
	double dx = Lx / (double)Nx;
	double dy = Ly / (double)Ny;
	double fxStart = -1.0 / (2.0 * dx);
	double fyStart = -1.0 / (2.0 * dy);
	double dfx = 1.0 / Lx;
	double dfy = 1.0 / Ly;

	double innerAngleFreqSqr = 1.0e-3 * innerAngle / wavLen;
	innerAngleFreqSqr *= innerAngleFreqSqr;

	double outerAngleFreqSqr = 1.0e-3 * outerAngle / wavLen;
	outerAngleFreqSqr *= outerAngleFreqSqr;

	for (int yIdx = 0; yIdx < Ny; yIdx++)
	{
		double fy = fyStart + (double)yIdx * dfy;
		for (int xIdx = 0; xIdx < Nx; xIdx++)
		{
			double fx = fxStart + (double)xIdx * dfx;
			double freqSqr = fx * fx + fy * fy;
			if ((freqSqr > (outerAngleFreqSqr + POSZERO)) ||
				(freqSqr < (innerAngleFreqSqr + NEGZERO)))
			{
				int i = yIdx * Nx + xIdx;
				if (otfMode == FFTSHIFTED_VEC)
					i = fftshift_Index_exchange_2D(i, Nx, Ny);
				otf[i] = { 0.0,0.0 };
			}
		}
	}
}


// Generate an electron probe:
void GenerateProbe(MKL_Complex16* otf, MKL_Complex16* probe, double xp, 
	double yp, double Lx, double Ly, int Nx, int Ny, int outputMode)
{
	double dx, dy, fx, fy, fxStart, fyStart, dfx, dfy, shiftFac, normCoeff;
	int ix, iy, i, VecLen = Ny * Nx;

	dx = Lx / (double)Nx;
	dy = Ly / (double)Ny;
	fxStart = -1.0 / (2.0 * dx);
	fyStart = -1.0 / (2.0 * dy);
	dfx = 1.0 / Lx;
	dfy = 1.0 / Ly;

	for (iy = 0; iy < Ny; iy++)
	{
		fy = fyStart + (double)iy * dfy;
		for (ix = 0; ix < Nx; ix++)
		{
			fx = fxStart + (double)ix * dfx;
			shiftFac = -2.0 * MathPI * (fx * xp + fy * yp);
			probe[iy * Nx + ix] = { cos(shiftFac),sin(shiftFac) };
		}
	}
	vzMul(VecLen, otf, probe, probe);
	fftshift_Z_2D_IP(probe, Nx, Ny);

	// Construct a fft handle:
	DFTI_DESCRIPTOR_HANDLE descHandle;
	MKL_LONG fftStat, matSize[2] = { (MKL_LONG)Ny,(MKL_LONG)Nx };
	fftStat = DftiCreateDescriptor(&descHandle, DFTI_DOUBLE, DFTI_COMPLEX, 2,
		matSize);
	fftStat = DftiCommitDescriptor(descHandle);
	fftStat = DftiComputeBackward(descHandle, probe);
	fftStat = DftiFreeDescriptor(&descHandle);

	normCoeff = 0.0;
	for (i = 0; i < VecLen; i++)
		normCoeff += probe[i].real * probe[i].real + probe[i].imag * probe[i].imag;
	normCoeff *= (dx * dy);
	normCoeff = sqrt(normCoeff);

	for (i = 0; i < VecLen; i++)
	{
		probe[i].real /= normCoeff;
		probe[i].imag /= normCoeff;
	}

	if (outputMode == UNFFTSHIFTED_VEC)
		ifftshift_Z_2D_IP(probe, Nx, Ny);
}


void InitFreqXYMesh(MKL_Complex16* fxMesh, double Lx, int Nx,
	MKL_Complex16* fyMesh, double Ly, int Ny)
{
	double dx = Lx / (double)Nx;
	double dy = Ly / (double)Ny;
	double fxStart = -1.0 / (2.0 * dx);
	double fyStart = -1.0 / (2.0 * dy);
	double dfx = 1.0 / Lx;
	double dfy = 1.0 / Ly;

	for (int yIdx = 0; yIdx < Ny; yIdx++)
	{
		double fy = fyStart + (double)yIdx * dfy;
		for (int xIdx = 0; xIdx < Nx; xIdx++)
		{
			double fx = fxStart + (double)xIdx * dfx;
			int i = yIdx * Nx + xIdx;
			fxMesh[i] = { 0.0,fx };
			fyMesh[i] = { 0.0,fy };
		}
	}
}


// save 3D ADF-STEM image
int Save3dImage(char* filename, double** stemImg, int sliceNum, 
	int vecLen)
{
	if ((sliceNum == 0) || (vecLen == 0))
	{
		PrintErrorMsg("Save3dImage", VTEMLAB_EINVAL);
		return VTEMLAB_EINVAL;
	}

	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "w+");
	if (err)
	{
		PrintErrorMsg("Save3dImage", VTEMLAB_EINVAL);
		return VTEMLAB_EINVAL;
	}

	for (int sliceIdx = 0; sliceIdx < sliceNum; sliceIdx++)
	{
		for (int i = 0; i < vecLen; i++)
			fprintf_s(destFile, "%lf\t", stemImg[sliceIdx][i]);
		if (sliceIdx != sliceNum - 1)
			fprintf_s(destFile, "\n");
	}

	fclose(destFile);

	return VTEMLAB_SUCCESS;
}


// Convert a depth list to slice index list
void DepthListToSliceIndexList(SliceList* sList, int sliceNumPerStack,
	int maxStackNum, double* depthList, int* sliceIndexList, int depthLevelNum)
{
	int sliceNum = sliceNumPerStack * maxStackNum;
	double* sliceDepthList = (double*)mkl_malloc(sliceNum * sizeof(double), 64);

	double tmpDepth = 0.0;
	int sliceIdx = 0;
	// Initialize sliceDepthList
	for (int stackIdx = 0; stackIdx < maxStackNum; stackIdx++)
	{
		SliceList* tmpSliceNode = sList->NextSlice;
		for (int inStackSliceIdx = 0; inStackSliceIdx < sliceNumPerStack; inStackSliceIdx++)
		{
			sliceIdx = stackIdx * sliceNumPerStack + inStackSliceIdx;
			tmpDepth += tmpSliceNode->SliceDist;
			sliceDepthList[sliceIdx] = tmpDepth;
			tmpSliceNode = tmpSliceNode->NextSlice;
		}
	}

	sliceIdx = 0;
	// Find the slice indices for each thickness level:
	for (int thickIdx = 0; thickIdx < depthLevelNum; thickIdx++)
	{
		double tmpThickness = depthList[thickIdx];
		while (sliceDepthList[sliceIdx] < tmpThickness)
		{
			if (sliceIdx == sliceNum - 1)
				break;

			sliceIdx++;
		}
		if (sliceIdx == 0)
			sliceIndexList[thickIdx] = sliceIdx + 1;
		else
		{
			double prevDist = fabs(tmpThickness - sliceDepthList[sliceIdx - 1]);
			double postDist = fabs(tmpThickness - sliceDepthList[sliceIdx]);
			sliceIndexList[thickIdx] = (prevDist < postDist) ? sliceIdx : (sliceIdx + 1);
		}
	}

	mkl_free(sliceDepthList);
}


int CbedUpdateSliceTF(MKL_Complex16* tf, int tfMode, double interCoeff,
	ThermoSliceList* tsList, int configIdx, int stackIdx, double bwlProp,
	double* scattFac, MKL_Complex16* fxMesh, MKL_Complex16* fyMesh,
	MKL_Complex16* fMesh, MKL_Complex16* convKer, MKL_Complex16* container,
	int Nx, int Ny)
{
	bool anyNullPtr = (tf == NULL) || (scattFac == NULL) || (fxMesh == NULL) ||
		(fyMesh == NULL) || (fMesh == NULL) || (convKer == NULL) ||
		(container == NULL);
	if (anyNullPtr)
	{
		PrintErrorMsg("CbedUpdateSliceTF", VTEMLAB_ENULLPTR);
		return VTEMLAB_ENULLPTR;
	}

	double scaledInterCoeff = 1.0e-3 * interCoeff;
	int vecLen = Ny * Nx;
	double Lx = (double)tsList->sliceExpanNum[0] * tsList->lattConstA;
	double Ly = (double)tsList->sliceExpanNum[1] * tsList->lattConstB;

	int errorCode = SliceComplexProjPot_thermo_conv(tsList, configIdx, stackIdx,
		scattFac, fxMesh, fyMesh, fMesh, convKer, tf, container, Nx, Ny, tfMode);
	if (errorCode)
	{
		PrintErrorMsg("CbedUpdateSliceTF", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	cblas_zdscal(vecLen, scaledInterCoeff, tf, 1);
	vzExp(vecLen, tf, tf);
	BwlSliceTF(tf, Lx, Ly, Nx, Ny, bwlProp, tfMode, tfMode);

	return VTEMLAB_SUCCESS;
}


// Dedicated multislice kernel for CBED simulation with TDS, cbedMat output in
// w_oMode
int CbedMultisliceKernel(MKL_Complex16* wave, int wIMode, int wOMode,
	double voltage, double bwlProp, ThermoSliceList* tsList, 
	PropKerList* pkList, int selectSliceNum, int* selectSliceIndices, 
	double** cbedMat, int Nx, int Ny)
{
	bool anyNullPtr = (wave == NULL) || (selectSliceIndices == NULL) ||
		(cbedMat == NULL);
	if (anyNullPtr)
	{
		PrintErrorMsg("CbedMultisliceKernel", VTEMLAB_ENULLPTR);
		return VTEMLAB_ENULLPTR;
	}

	if (wIMode == UNFFTSHIFTED_VEC)
		fftshift_Z_2D_IP(wave, Nx, Ny);

	int vecLen = Ny * Nx;
	double invVecLen = 1.0 / (double)vecLen;
	double interCoeff = InteractionCoefficient(voltage);
	ThermoSliceList* tmpTsList = tsList->nextSlice;
	double Lx = (double)tmpTsList->sliceExpanNum[0] * tmpTsList->lattConstA;
	double Ly = (double)tmpTsList->sliceExpanNum[1] * tmpTsList->lattConstB;

	// Initialize cbedMat as all zeros:
	for (int sliceIdx = 0; sliceIdx < selectSliceNum; sliceIdx++)
	{
		for (int i = 0; i < vecLen; i++)
			cbedMat[sliceIdx][i] = 0.0;
	}

	// initialize tf
	MKL_Complex16* tf =
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (tf == NULL)
	{
		PrintErrorMsg("CbedMultisliceKernel", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}

	// Initialize fxMesh (complex16*), fyMesh (complex16*), fMesh (complex16*)
	MKL_Complex16* fxMesh =
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (fxMesh == NULL)
	{
		mkl_free(tf);
		PrintErrorMsg("CbedMultisliceKernel", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}

	MKL_Complex16* fyMesh =
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (fyMesh == NULL)
	{
		mkl_free(tf);
		mkl_free(fxMesh);
		PrintErrorMsg("CbedMultisliceKernel", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}

	InitFreqXYMesh(fxMesh, Lx, Nx, fyMesh, Ly, Ny);

	MKL_Complex16* fMesh =
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (fMesh == NULL)
	{
		mkl_free(tf);
		mkl_free(fxMesh);
		mkl_free(fyMesh);
		PrintErrorMsg("CbedMultisliceKernel", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}

	// initialize convKer (complex16*), tmpProjPot (complex16*)
	MKL_Complex16* convKer =
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (convKer == NULL)
	{
		mkl_free(tf);
		mkl_free(fxMesh);
		mkl_free(fyMesh);
		mkl_free(fMesh);
		PrintErrorMsg("CbedMultisliceKernel", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}

	MKL_Complex16* tmpProjPot =
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (tmpProjPot == NULL)
	{
		mkl_free(tf);
		mkl_free(fxMesh);
		mkl_free(fyMesh);
		mkl_free(fMesh);
		mkl_free(convKer);
		PrintErrorMsg("CbedMultisliceKernel", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}

	// initialize waveI (complex16)
	double* waveI = (double*)mkl_malloc(vecLen * sizeof(double), 64);
	if (waveI == NULL)
	{
		mkl_free(tf);
		mkl_free(fxMesh);
		mkl_free(fyMesh);
		mkl_free(fMesh);
		mkl_free(convKer);
		mkl_free(tmpProjPot);
		PrintErrorMsg("CbedMultisliceKernel", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}

	// initialize a initWave to store the initial incident wave:
	MKL_Complex16* initWave =
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (initWave == NULL)
	{
		mkl_free(tf);
		mkl_free(fxMesh);
		mkl_free(fyMesh);
		mkl_free(fMesh);
		mkl_free(convKer);
		mkl_free(tmpProjPot);
		mkl_free(waveI);
		PrintErrorMsg("CbedMultisliceKernel", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}
	cblas_zcopy(vecLen, wave, 1, initWave, 1);

	// Load scattFac
	double* scattFac;
	int errorCode = LoadScattFacTable(scattFac);
	if (errorCode)
	{
		mkl_free(tf);
		mkl_free(fxMesh);
		mkl_free(fyMesh);
		mkl_free(fMesh);
		mkl_free(convKer);
		mkl_free(tmpProjPot);
		mkl_free(waveI);
		mkl_free(initWave);
		PrintErrorMsg("CbedMultisliceKernel", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	// Constructing fft descriptor:
	DFTI_DESCRIPTOR_HANDLE fftDescHandle;
	MKL_LONG fftStat, matSize[2] = { (MKL_LONG)Ny,(MKL_LONG)Nx };
	fftStat = DftiCreateDescriptor(&fftDescHandle, DFTI_DOUBLE, DFTI_COMPLEX,
		2, matSize);
	fftStat = DftiCommitDescriptor(fftDescHandle);

	// loop over all configs and stacks:
	int configNum = tmpTsList->configNum;
	int stackNum = tmpTsList->stackNum;

	MKL_Complex16* recordWave =
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (recordWave == NULL)
	{
		mkl_free(tf);
		mkl_free(fxMesh);
		mkl_free(fyMesh);
		mkl_free(fMesh);
		mkl_free(convKer);
		mkl_free(tmpProjPot);
		mkl_free(waveI);
		mkl_free(scattFac);
		mkl_free(initWave);
		PrintErrorMsg("CbedMultisliceKernel", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}

	int sliceCount = 0;
	int tmpRecordSliceIdx = 0;
	for (int configIdx = 0; configIdx < configNum; configIdx++)
	{
		UpdateProcess("CbedMultisliceKernel", configIdx, configNum);
		cblas_zcopy(vecLen, initWave, 1, wave, 1);
		for (int stackIdx = 0; stackIdx < stackNum; stackIdx++)
		{
			tmpTsList = tsList->nextSlice;
			PropKerList* tmpPkList = pkList->nextSlice;
			while (tmpTsList != NULL)
			{
				errorCode = CbedUpdateSliceTF(tf, FFTSHIFTED_VEC, interCoeff,
					tmpTsList, configIdx, stackIdx, bwlProp, scattFac, fxMesh,
					fyMesh, fMesh, convKer, tmpProjPot, Nx, Ny);
				if (errorCode)
				{
					mkl_free(tf);
					mkl_free(fxMesh);
					mkl_free(fyMesh);
					mkl_free(fMesh);
					mkl_free(convKer);
					mkl_free(tmpProjPot);
					mkl_free(waveI);
					mkl_free(scattFac);
					mkl_free(initWave);
					mkl_free(recordWave);
					PrintErrorMsg("CbedMultisliceKernel", VTEMLAB_FAILURE);
					return VTEMLAB_FAILURE;
				}

				vzMul(vecLen, wave, tf, wave);
				fftStat = DftiComputeForward(fftDescHandle, wave);
				vzMul(vecLen, tmpPkList->propKer, wave, wave);
				fftStat = DftiComputeBackward(fftDescHandle, wave);
				cblas_zdscal(vecLen, invVecLen, wave, 1);

				tmpTsList = tmpTsList->nextSlice;
				tmpPkList = tmpPkList->nextSlice;

				sliceCount++;
				if (sliceCount == selectSliceIndices[tmpRecordSliceIdx])
				{
					cblas_zcopy(vecLen, wave, 1, recordWave, 1);
					fftStat = DftiComputeForward(fftDescHandle, recordWave);
					GetWaveInten(waveI, recordWave, vecLen);
					vdAdd(vecLen, waveI, cbedMat[tmpRecordSliceIdx],
						cbedMat[tmpRecordSliceIdx]);

					tmpRecordSliceIdx++;
				}
			}
		}

		UpdateProcess("CbedMultisliceKernel", configIdx + 1, configNum);
	}
	PrintErrorMsg("CbedMultisliceKernel", VTEMLAB_SUCCESS);

	double invConfigNum = 1.0 / (double)configNum;
	for (int sliceIdx = 0; sliceIdx < selectSliceNum; sliceIdx++)
	{
		cblas_dscal(vecLen, invConfigNum, cbedMat[sliceIdx], 1);
		if (wOMode == UNFFTSHIFTED_VEC)
			ifftshift_D_2D_IP(cbedMat[sliceIdx], Nx, Ny);
	}

	fftStat = DftiFreeDescriptor(&fftDescHandle);

	if (wOMode == UNFFTSHIFTED_VEC)
		ifftshift_Z_2D_IP(wave, Nx, Ny);

	mkl_free(tf);
	mkl_free(fxMesh);
	mkl_free(fyMesh);
	mkl_free(fMesh);
	mkl_free(convKer);
	mkl_free(tmpProjPot);
	mkl_free(waveI);
	mkl_free(scattFac);
	mkl_free(initWave);
	mkl_free(recordWave);

	return VTEMLAB_SUCCESS;
}


int CbedTdsKernel(OtfParamSet otfParam, double voltage, double numApert,
	ThermoSliceList* tsList, PropKerList* pkList, double xp, double yp,
	double bwlProp, int Nx, int Ny, int selectSliceNum, int* selectSliceIndices,
	char* cbedFilename)
{
	double wavLen = EleWavLen(voltage);
	double interCoeff = InteractionCoefficient(voltage);

	int vecLen = Ny * Nx;
	ThermoSliceList* tmpTsList = tsList->nextSlice;
	double Lx = (double)tmpTsList->sliceExpanNum[0] * tmpTsList->lattConstA;
	double Ly = (double)tmpTsList->sliceExpanNum[1] * tmpTsList->lattConstB;

	// Initialize OTF
	MKL_Complex16* otf =
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (otf == NULL)
	{
		PrintErrorMsg("CbedTdsKernel", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}
	int errorCode = ObjTransFunc(otf, otfParam, wavLen, Lx, Ly, Nx, Ny);
	if (errorCode)
	{
		mkl_free(otf);
		PrintErrorMsg("CbedTdsKernel", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	AddCircApert(otf, wavLen, numApert, Lx, Ly, Nx, Ny);

	// initialize wave (complex16*), waveI (double*)
	MKL_Complex16* wave =
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (wave == NULL)
	{
		mkl_free(otf);
		PrintErrorMsg("CbedTdsKernel", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}

	GenerateProbe(otf, wave, xp, yp, Lx, Ly, Nx, Ny, FFTSHIFTED_VEC);

	// initialize cbedMat
	double** cbedMat =
		(double**)mkl_malloc(selectSliceNum * sizeof(double*), 64);
	if (cbedMat == NULL)
	{
		mkl_free(otf);
		mkl_free(wave);
		PrintErrorMsg("CbedTdsKernel", VTEMLAB_ENOMEM);
		return VTEMLAB_ENOMEM;
	}

	for (int sliceIdx = 0; sliceIdx < selectSliceNum; sliceIdx++)
	{
		cbedMat[sliceIdx] = (double*)mkl_malloc(vecLen * sizeof(double), 64);
		if (cbedMat[sliceIdx] == NULL)
		{
			for (int freeSliceIdx = 0; freeSliceIdx < sliceIdx; freeSliceIdx++)
				mkl_free(cbedMat[freeSliceIdx]);
			mkl_free(cbedMat);
			mkl_free(otf);
			mkl_free(wave);
			PrintErrorMsg("CbedTdsKernel", VTEMLAB_ENOMEM);
			return VTEMLAB_ENOMEM;
		}
	}

	errorCode = CbedMultisliceKernel(wave, FFTSHIFTED_VEC, UNFFTSHIFTED_VEC,
		voltage, bwlProp, tsList, pkList, selectSliceNum, selectSliceIndices,
		cbedMat, Nx, Ny);
	if (errorCode)
	{
		for (int sliceIdx = 0; sliceIdx < selectSliceNum; sliceIdx++)
			mkl_free(cbedMat[sliceIdx]);
		mkl_free(cbedMat);
		mkl_free(otf);
		mkl_free(wave);
		PrintErrorMsg("CbedTdsKernel", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	errorCode = Save3dImage(cbedFilename, cbedMat, selectSliceNum, 
		vecLen);
	if (errorCode)
	{
		for (int sliceIdx = 0; sliceIdx < selectSliceNum; sliceIdx++)
			mkl_free(cbedMat[sliceIdx]);
		mkl_free(cbedMat);
		mkl_free(otf);
		mkl_free(wave);
		PrintErrorMsg("CbedTdsKernel", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	// free memory:
	for (int sliceIdx = 0; sliceIdx < selectSliceNum; sliceIdx++)
		mkl_free(cbedMat[sliceIdx]);
	mkl_free(cbedMat);
	mkl_free(otf);
	mkl_free(wave);

	return VTEMLAB_SUCCESS;
}