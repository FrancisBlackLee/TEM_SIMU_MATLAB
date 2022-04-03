/*
	test.h is the open source part of vtemlab v0.0 engine,
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

#ifndef TEST_H
#define TEST_H

// mkl dependence
#include <mkl.h>

// Test PrintErrorMsg
void TestPrintErrorMsg();

// Test UpdateProcess
void TestUpdateProcess();

// Simple test of MKL:
int TestMklZvecAdd();

// Test the CreateSliceList(): part1 generates test data, 
// part2 performs the test using data generated in part1
int TestCreateSliceList_p1();
int TestCreateSliceList_p2();

// Test removing repeated atoms on the boundary and transforming SliceNode 
// to SliceTypeList:
int TestSliceNodeToTypeList();

// Test the loading of scattering factors:
int TestLoadScattFac();

// Test CrysFileToProjPot()
int TestCrysFileToProjPot_conv_Frac();

// Test CbedTdsKernel
int TestCbedTdsKernel();

int TestInteract();

#endif // !TEST_H
