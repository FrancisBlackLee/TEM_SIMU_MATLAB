/*
	test.cpp is the open source part of vtemlab v0.0 engine,
	providing unit tests for vtemlab.

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
#include <Windows.h>
#include <math.h>
#include <time.h>
#include <direct.h>

#include "interact.h"
#include "mathconst.h"
#include "mode.h"
#include "multislice.h"
#include "optics.h"
#include "potential.h"
#include "slice.h"
#include "status.h"
#include "test.h"
#include "thermo.h"
#include "vtemlabio.h"
#include "vtemlabmath.h"

// Test PrintErrorMsg
void TestPrintErrorMsg()
{
	PrintErrorMsg("TestPrintErrorMsg", VTEMLAB_SUCCESS);
	PrintErrorMsg("TestPrintErrorMsg", VTEMLAB_FAILURE);
	PrintErrorMsg("TestPrintErrorMsg", VTEMLAB_ENULLPTR);
	PrintErrorMsg("TestPrintErrorMsg", VTEMLAB_ENOMEM);
	PrintErrorMsg("TestPrintErrorMsg", VTEMLAB_ELOSS);
	PrintErrorMsg("TestPrintErrorMsg", VTEMLAB_EFILEW);
	PrintErrorMsg("TestPrintErrorMsg", VTEMLAB_EFILER);
	PrintErrorMsg("TestPrintErrorMsg", VTEMLAB_ENOENT);
	PrintErrorMsg("TestPrintErrorMsg", VTEMLAB_EINVAL);
	PrintErrorMsg("TestPrintErrorMsg", -13);
}


// Test UpdateProcess
void TestUpdateProcess()
{
	int totalNum = 10;
	for (int doneNum = 0; doneNum < totalNum; doneNum++)
	{
		UpdateProcess("Test process", doneNum, totalNum);
		Sleep(1000);
	}

	UpdateProcess("Test process", totalNum, totalNum);
	printf("\n\n");
}


// Simple test of MKL:
int TestMklZvecAdd()
{
	int vecLen = 10;
	MKL_Complex16* x1 = 
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (x1 == NULL)
	{
		PrintErrorMsg("TestMklZvecAdd", VTEMLAB_ENULLPTR);
		return VTEMLAB_ENULLPTR;
	}

	MKL_Complex16* x2 = 
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (x2 == NULL)
	{
		PrintErrorMsg("TestMklZvecAdd", VTEMLAB_ENULLPTR);
		mkl_free(x1);
		return VTEMLAB_ENULLPTR;
	}

	MKL_Complex16* y = 
		(MKL_Complex16*)mkl_malloc(vecLen * sizeof(MKL_Complex16), 64);
	if (y == NULL)
	{
		PrintErrorMsg("TestMklZvecAdd", VTEMLAB_ENULLPTR);
		mkl_free(x1);
		mkl_free(x2);
		return VTEMLAB_ENULLPTR;
	}

	for (int i = 0; i < vecLen; i++)
	{
		x1[i] = { (double)i,(double)(i + 1) };
		x2[i] = { (double)i / 2.0,(double)i / 3.0 };
	}

	vzAdd(vecLen, x1, x2, y);
	double sigma = 0.0;
	for (int i = 0; i < vecLen; i++)
	{
		double realDiff = y[i].real - ((double)i + (double)i / 2.0);
		double imagDiff = y[i].imag - ((double)(i + 1) + (double)i / 3.0);
		sigma += realDiff * realDiff + imagDiff * imagDiff;
	}

	sigma /= (double)vecLen;
	sigma = sqrt(sigma);
	double tolerance = 1.0e-8;
	if (sigma > tolerance)
	{
		PrintErrorMsg("TestMklZvecAdd", VTEMLAB_ELOSS);
		mkl_free(x1);
		mkl_free(x2);
		mkl_free(y);
		return VTEMLAB_ELOSS;
	}

	PrintErrorMsg("TestMklZvecAdd", VTEMLAB_SUCCESS);
	mkl_free(x1);
	mkl_free(x2);
	mkl_free(y);

	return VTEMLAB_SUCCESS;
}


// Test the CreateSliceList(): part1 generates test data, part2 performs the test using data generated in part1
int TestCreateSliceList_p1()
{
	int sliceNum, atomNum, sliceIdx, atomIdx, tmpInt;
	double tmpDouble;
	char destDir[] = "test_slices";

	srand((unsigned int)time(NULL));
	sliceNum = 5;

	errno_t err = _mkdir(destDir);
	if (err == ENOENT)
	{
		PrintErrorMsg("TestCreateSliceList_p1", VTEMLAB_ENOENT);
		return VTEMLAB_ENOENT;
	}
	else
	{
		for (sliceIdx = 0; sliceIdx < sliceNum; sliceIdx++)
		{
			char filename[FILENAME_MAX];
			snprintf(filename, FILENAME_MAX, "%s\\s%d.txt", destDir, sliceIdx + 1);
			puts(filename);

			FILE* destFile;
			err = fopen_s(&destFile, filename, "w+");
			if (err)
			{
				PrintErrorMsg("TestCreateSliceList_p1", VTEMLAB_EFILEW);
			}

			// Generate random LattConstA, LattConstB and SliceDist:
			tmpDouble = RandomDouble(1, 5);
			fprintf_s(destFile, "%.4f\t", tmpDouble);
			tmpDouble = RandomDouble(1, 5);
			fprintf_s(destFile, "%.4f\t", tmpDouble);
			tmpDouble = RandomDouble(1, 5);
			fprintf_s(destFile, "%.4f\n", tmpDouble);

			// Generate random AtomNum:
			atomNum = RandomInt(5, 20);
			fprintf_s(destFile, "%d\n", atomNum);
			// assign random values to AtomType, EleProp, CoordX and CoordY for each atom:
			for (atomIdx = 0; atomIdx < atomNum; atomIdx++)
			{
				tmpInt = RandomInt(10, 50);
				fprintf_s(destFile, "%d\t", tmpInt);
				tmpDouble = RandomDouble(0.5, 1);
				fprintf_s(destFile, "%.4f\t", tmpDouble);
				tmpDouble = RandomDouble(0, 1);
				fprintf_s(destFile, "%.4f\t", tmpDouble);
				tmpDouble = RandomDouble(0, 1);
				fprintf_s(destFile, "%.4f\n", tmpDouble);
			}
			fclose(destFile);
		}
	}

	return VTEMLAB_SUCCESS;
}

int TestCreateSliceList_p2()
{
	int errorCode, sliceNum = 0;
	char destDir[] = "test_slices";

	// Start test:
	SliceList* testSliceList;
	errorCode = CreateSliceList(destDir, testSliceList, sliceNum);

	if (errorCode)
	{
		PrintErrorMsg("TestCreateSliceList_p2", VTEMLAB_FAILURE);
		FreeSliceList(testSliceList);
		return VTEMLAB_FAILURE;
	}
	printf("%d slices detected:\n\n", sliceNum);
	ScanSliceList(testSliceList, sliceNum);

	FreeSliceList(testSliceList);

	return VTEMLAB_SUCCESS;
}

// Test removing repeated atoms on the boundary and transforming SliceNode 
// to SliceTypeList:
int TestSliceNodeToTypeList()
{
	char filename[] = "test_slices\\s1.txt";
	SliceList* testSliceNode;
	SliceTypeList* testTypeList;

	// Initializing SliceNode:
	int errorCode;
	errorCode = CreateSliceNode(filename, testSliceNode);
	if (errorCode)
	{
		PrintErrorMsg("TestSliceNodeToTypeList", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	// Transfer data from the SliceNode to TypeList:
	errorCode = SliceNodeToTypeList_Frac(testSliceNode, testTypeList, 1.0e-4, 1.0e-4);
	if (errorCode)
	{
		PrintErrorMsg("TestSliceNodeToTypeList", VTEMLAB_FAILURE);
		free(testSliceNode);
		return VTEMLAB_FAILURE;
	}
	ScanSliceTypeList(testTypeList);

	free(testSliceNode);
	FreeSliceTypeList(testTypeList);

	return VTEMLAB_SUCCESS;
}


// Test the loading of scattering factors:
int TestLoadScattFac()
{
	double* scattParams;
	int i, errorCode = LoadScattFacTable(scattParams);
	if (errorCode)
	{
		PrintErrorMsg("TestLoadScattFac", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	printf("Elemental scattering factors (EJK):\n");

	for (i = 0; i < 103 * 12; i++)
	{
		if (i % 12 == 0)
			printf("**************************************************************\n"
				"Z = %d\n", i / 12 + 1);
		printf("%f\t", scattParams[i]);
		if (i % 4 == 3)
			printf("\n");
	}

	mkl_free(scattParams);

	return VTEMLAB_SUCCESS;
}


// Test CrysFileToProjPot()
int TestCrysFileToProjPot_conv_Frac()
{
	char folder[] = "si_110";
	int cellNum[2] = { 3,2 }, errorCode;
	int Nx = 512, Ny = 512;

	errorCode = CrysFileToProjPot_conv(folder, cellNum, Nx, Ny, FRAC_COORD);
	if (errorCode)
	{
		PrintErrorMsg("Test_CrysFileToProjPot_conv_Frac", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	errorCode = DecryptProjPotFiles(folder, 2);
	if (errorCode)
	{
		PrintErrorMsg("Test_CrysFileToProjPot_conv_Frac", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	PrintErrorMsg("Test_CrysFileToProjPot_conv_Frac", errorCode);

	return VTEMLAB_SUCCESS;
}


// Test CbedTdsKernel
int TestCbedTdsKernel()
{
	char folder[] = "si_111";

	int expanNum[2] = { 12,7 };
	int configNum = 4;
	int stackNum = 107;
	int Nx = 512, Ny = 512;
	int typeNum = 1;
	int typeList[1] = { 14 };
	double bFactors_293K[1] = { 0.4668 };
	int coordType = FRAC_COORD;
	double simuTemp = 300.0;

	// prepare thermo sample:
	SliceList* testSliceList;
	ThermoSliceList* testThermoSliceList;
	int sliceNum = 0;
	int errorCode = CreateThermoSample(folder, coordType, expanNum, testSliceList,
		testThermoSliceList, configNum, stackNum, sliceNum, typeNum, typeList,
		bFactors_293K, simuTemp);
	if (errorCode)
	{
		PrintErrorMsg("TestCbedTdsKernel", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	// create fftshifted PropKerList:
	double voltage = 100.0;
	double wavLen = EleWavLen(voltage);
	PropKerList* testPropKerList;
	bool sliceDistEqual = true;
	errorCode = AutoCreatePropKerList(testSliceList, testPropKerList, wavLen, expanNum,
		Nx, Ny, FFTSHIFTED_VEC, sliceDistEqual);
	if (errorCode)
	{
		FreeSliceList(testSliceList);
		FreeThermoSliceList(testThermoSliceList);
		PrintErrorMsg("TestCbedTdsKernel", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	// CBED:
	OtfParamSet otfParam;
	otfParam.Cs3 = 0.0;
	otfParam.Cs5 = 0.0;
	otfParam.defocus = 0.0;

	double numApert = 8;
	double bwlProp = 0.67;

	int depthNum = 51;
	double* depthList = (double*)mkl_malloc(depthNum * sizeof(double), 64);
	int* sliceIndexList = (int*)malloc(depthNum * sizeof(int));
	for (int depthIdx = 0; depthIdx < depthNum; depthIdx++)
	{
		depthList[depthIdx] = 500.0 + (double)depthIdx * 10.0;
		sliceIndexList[depthIdx] = 0;
	}

	DepthListToSliceIndexList(testSliceList, sliceNum, stackNum, depthList,
		sliceIndexList, depthNum);

	printf("> select data output slice:\n");
	for (int sliceIdx = 0; sliceIdx < depthNum; sliceIdx++)
	{
		printf("%d\t(%.4f Angs.)\n", sliceIndexList[sliceIdx], depthList[sliceIdx]);
	}
	printf("\n\n");

	char cbedFilename[FILENAME_MAX];
	snprintf(cbedFilename, FILENAME_MAX, "%s\\test_cbed.txt", folder);

	errorCode = CbedTdsKernel(otfParam, voltage, numApert, testThermoSliceList,
		testPropKerList, 0.0, 0.0, bwlProp, Nx, Ny, depthNum, sliceIndexList,
		cbedFilename);
	if (errorCode)
	{
		FreeSliceList(testSliceList);
		FreeThermoSliceList(testThermoSliceList);
		FreePropKerList(testPropKerList, sliceDistEqual);
		mkl_free(depthList);
		free(sliceIndexList);
		PrintErrorMsg("TestCbedTdsKernel", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}

	FreeSliceList(testSliceList);
	FreeThermoSliceList(testThermoSliceList);
	FreePropKerList(testPropKerList, sliceDistEqual);
	mkl_free(depthList);
	free(sliceIndexList);

	return VTEMLAB_SUCCESS;
}


int TestInteract()
{
	PrintLicense();

	int testInt;
	AssignInt("integer", testInt);
	printf("val = %d\n", testInt);

	double testDouble;
	AssignDouble("double", testDouble);
	printf("val = %.4f\n", testDouble);

	char filename[FILENAME_MAX];
	AssignString("filename", filename);
	printf("val = %s\n", filename);

	int vecLen = 10;
	double* testDoubleArray = (double*)mkl_malloc(vecLen * sizeof(double), 64);
	AssignDoubleArray("double array", testDoubleArray, vecLen);
	for (int i = 0; i < vecLen; i++)
		printf("val[%d] = %.4f\n", i, testDoubleArray[i]);

	mkl_free(testDoubleArray);

	return VTEMLAB_SUCCESS;
}