/*
	interact.h is the open source part of vtemlab v0.0 engine,
	providing interactive operations for CbedTdsConsole.

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

#ifndef INTERACT_H
#define INTERACT_H

void PrintLicense();

void AssignInt(const char* item, int& val);

void AssignDouble(const char* item, double& val);

void AssignString(const char* item, char* val);

void AssignIntArray(const char* item, int* val, int vecLen);

void AssignDoubleArray(const char* item, double* val, int vecLen);

// Load the name of elements:
const char* LoadElementName(int atomType);

void AssignBFactors(int* atomType, double* bFactors, int vecLen);

#endif // !INTERACT_H
