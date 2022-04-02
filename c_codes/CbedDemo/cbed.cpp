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

#include "status.h"
#include "test.h"

#define RUN_TEST

int main()
{
	// TEST CODE:
#ifdef RUN_TEST
	// run test cases
	TestPrintErrorMsg();
	TestUpdateProcess();
	TestMklZvecAdd();
	TestCreateSliceList_p1();
	TestCreateSliceList_p2();
	TestSliceNodeToTypeList();
	TestLoadScattFac();
	TestCrysFileToProjPot_conv_Frac();
	TestCbedTdsKernel();
	// blabla
#else
	// run CBED module
#endif

	return 0;
}
