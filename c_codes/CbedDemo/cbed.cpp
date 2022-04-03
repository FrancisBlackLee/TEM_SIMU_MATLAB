/*
	cbed.cpp is the open source part of vtemlab v0.0 engine,
	providing CBED console or tests for vtemlab.

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
#include <string.h>
#include <math.h>

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

#define RUN_BASIC_TEST
#define RUN_KERNEL_TEST
#define RUN_CONSOLE

int CbedTdsConsole();

int main()
{
#ifdef RUN_CONSOLE
	CbedTdsConsole();
#else

#ifdef RUN_BASIC_TEST
	// run basic test cases
	TestPrintErrorMsg();
	TestUpdateProcess();
	TestMklZvecAdd();
	TestCreateSliceList_p1();
	TestCreateSliceList_p2();
	TestSliceNodeToTypeList();
	TestLoadScattFac();
	TestCrysFileToProjPot_conv_Frac();
	TestInteract();
	// blabla
#endif // RUN_BASIC_TEST

#ifdef RUN_KERNEL_TEST
	TestCbedTdsKernel();
#endif // RUN_KERNEL_TEST


#endif // RUN_CONSOLE

	return 0;
}


int CbedTdsConsole()
{
	PrintLicense();
	
	char folder[FILENAME_MAX];
	AssignString("slice folder", folder);
	SliceList* testSList;
	int sliceNum;
	int errorCode = CreateSliceList(folder, testSList, sliceNum);
	if (errorCode)
	{
		PrintErrorMsg("CbedTdsConsole", VTEMLAB_FAILURE);
		return VTEMLAB_FAILURE;
	}
	double stackThickness = GetStackThickness(testSList);
	ScanSliceList(testSList, sliceNum);

	int expanNum[2] = { 1,1 };
	AssignInt("unit cell number along x", expanNum[0]);
	AssignInt("unit cell number along y", expanNum[1]);

	double Lx = (double)expanNum[0] * testSList->NextSlice->LattConstA;
	double Ly = (double)expanNum[1] * testSList->NextSlice->LattConstB;
	printf("Specimen size: Lx = %.4f, Ly = %.4f\n\n", Lx, Ly);

	int depthNum = 0;
	AssignInt("depth number", depthNum);
	double* depthList = (double*)mkl_malloc(depthNum * sizeof(double), 64);
	if (depthList == NULL)
	{
		PrintErrorMsg("CbedTdsConsole", VTEMLAB_ENOMEM);
		FreeSliceList(testSList);
		return VTEMLAB_ENOMEM;
	}
	AssignDoubleArray("depths (in Angstrom, ascending order)", 
		depthList, depthNum);

	int stackNum = (int)ceil(depthList[depthNum - 1] / stackThickness);
	printf("Unit cell tiled (z) by %d\n\n", stackNum);

	int* sliceIndexList = (int*)malloc(depthNum * sizeof(int));
	if (sliceIndexList == NULL)
	{
		PrintErrorMsg("CbedTdsConsole", VTEMLAB_ENOMEM);
		FreeSliceList(testSList);
		mkl_free(depthList);
		return VTEMLAB_ENOMEM;
	}
	DepthListToSliceIndexList(testSList, sliceNum, stackNum, depthList,
		sliceIndexList, depthNum);

	printf("> Output slice positions:\n");
	for (int i = 0; i < depthNum; i++)
	{
		printf("%d\t(%.4f Angs.)\n", sliceIndexList[i], depthList[i]);
	}
	printf("\n\n");

	int Nx = 512, Ny = 512;
	AssignInt("pixel number along x", Nx);
	AssignInt("pixel number along y", Ny);

	int configNum = 4;
	AssignInt("number of frozen configurations", configNum);

	int typeNum = 1;
	AssignInt("number of element types", typeNum);
	int* typeList = (int*)malloc(typeNum * sizeof(int));
	if (typeList == NULL)
	{
		PrintErrorMsg("CbedTdsConsole", VTEMLAB_ENOMEM);
		FreeSliceList(testSList);
		mkl_free(depthList);
		free(sliceIndexList);
		return VTEMLAB_ENOMEM;
	}

	AssignIntArray("element types (Z)", typeList, typeNum);
	double* bFactors = (double*)mkl_malloc(typeNum * sizeof(double), 64);
	if (bFactors == NULL)
	{
		PrintErrorMsg("CbedTdsConsole", VTEMLAB_ENOMEM);
		FreeSliceList(testSList);
		mkl_free(depthList);
		free(sliceIndexList);
		free(typeList);
		return VTEMLAB_ENOMEM;
	}

	AssignBFactors(typeList, bFactors, typeNum);

	double simuTemp = 293.0;
	AssignDouble("simulation temperature", simuTemp);

	// prepare specimen with TDS:
	SliceList* sList;
	ThermoSliceList* tsList;
	errorCode = CreateThermoSample(folder, FRAC_COORD, expanNum, sList, tsList,
		configNum, stackNum, sliceNum, typeNum, typeList, bFactors, simuTemp);
	if (errorCode)
	{
		PrintErrorMsg("CbedTdsConsole", VTEMLAB_FAILURE);
		FreeSliceList(testSList);
		mkl_free(depthList);
		free(sliceIndexList);
		free(typeList);
		mkl_free(bFactors);
		return VTEMLAB_FAILURE;
	}

	double voltage = 300.0;
	AssignDouble("voltage (kV)", voltage);

	double wavLen = EleWavLen(voltage);
	PropKerList* pkList;
	bool sliceDistEqual = true;
	errorCode = AutoCreatePropKerList(sList, pkList, wavLen, expanNum, Nx, Ny, 
		FFTSHIFTED_VEC, sliceDistEqual);
	if (errorCode)
	{
		PrintErrorMsg("CbedTdsConsole", VTEMLAB_FAILURE);
		FreeSliceList(testSList);
		mkl_free(depthList);
		free(sliceIndexList);
		free(typeList);
		mkl_free(bFactors);
		FreeSliceList(sList);
		FreeThermoSliceList(tsList);
		return VTEMLAB_FAILURE;
	}

	OtfParamSet otfParam;
	otfParam.Cs3 = 0.0;
	otfParam.Cs5 = 0.0;
	otfParam.defocus = 0.0;
	AssignDouble("Cs3", otfParam.Cs3);
	AssignDouble("Cs5", otfParam.Cs5);
	AssignDouble("defocus", otfParam.defocus);

	double numApert = 7.5;
	AssignDouble("convergent semi-angle", numApert);
	double bwlProp = 0.67;

	double xp = 0.0, yp = 0.0;
	AssignDouble("probe position (x, cartesian)", xp);
	AssignDouble("probe position (x, cartesian)", yp);

	char filename[FILENAME_MAX];
	AssignString("output filename", filename);

	errorCode = CbedTdsKernel(otfParam, voltage, numApert, tsList, pkList, 
		xp, yp, bwlProp, Nx, Ny, depthNum, sliceIndexList, filename);
	if (errorCode)
	{
		PrintErrorMsg("CbedTdsConsole", VTEMLAB_FAILURE);
		FreeSliceList(testSList);
		mkl_free(depthList);
		free(sliceIndexList);
		free(typeList);
		mkl_free(bFactors);
		FreeSliceList(sList);
		FreeThermoSliceList(tsList);
		FreePropKerList(pkList, sliceDistEqual);
		return VTEMLAB_FAILURE;
	}

	// free memory:
	FreeSliceList(testSList);
	mkl_free(depthList);
	free(sliceIndexList);
	free(typeList);
	mkl_free(bFactors);
	FreeSliceList(sList);
	FreeThermoSliceList(tsList);
	FreePropKerList(pkList, sliceDistEqual);

	return VTEMLAB_SUCCESS;
}