/*
	interact.cpp is the open source part of vtemlab v0.0 engine,
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

#include <stdio.h>

#include "interact.h"
#include "slice.h"
#include "status.h"

void PrintLicense()
{
	printf("CbedDemo (vtemlab module)  Copyright (C) 2022  "
		"Francis Black Lee (Li Xian)\n"
		"This program comes with ABSOLUTELY NO WARRANTY;\n"
		"This is free software, and you are welcome to \n"
		"redistribute it under certain conditions.\n\n");
}

void AssignInt(const char* item, int& val)
{
	printf("> Enter %s\n", item);
	scanf_s("%d", &val);
	getc(stdin);
}

void AssignDouble(const char* item, double& val)
{
	printf("> Enter %s\n", item);
	scanf_s("%lf", &val);
	getc(stdin);
}

void AssignString(const char* item, char* val)
{
	printf("> Enter %s\n", item);
	gets_s(val, FILENAME_MAX);
}

void AssignIntArray(const char* item, int* val, int vecLen)
{
	printf("> Enter %d items for %s\n", vecLen, item);
	for (int i = 0; i < vecLen; i++)
	{
		printf("item %d: ", i);
		scanf_s("%d", &val[i]);
	}
	getc(stdin);
}

void AssignDoubleArray(const char* item, double* val, int vecLen)
{
	printf("> Enter %d items for %s\n", vecLen, item);
	for (int i = 0; i < vecLen; i++)
	{
		printf("item %d: ", i);
		scanf_s("%lf", &val[i]);
	}
	getc(stdin);
}

// Load the name of elements:
const char* LoadElementName(int atomType)
{
	switch (atomType)
	{
	case 1:
		return "H"; break;
	case 2:
		return "He"; break;
	case 3:
		return "Li"; break;
	case 4:
		return "Be"; break;
	case 5:
		return "B"; break;
	case 6:
		return "C"; break;
	case 7:
		return "N"; break;
	case 8:
		return "O"; break;
	case 9:
		return "F"; break;
	case 10:
		return "Ne"; break;
	case 11:
		return "Na"; break;
	case 12:
		return "Mg"; break;
	case 13:
		return "Al"; break;
	case 14:
		return "Si"; break;
	case 15:
		return "P"; break;
	case 16:
		return "S"; break;
	case 17:
		return "Cl"; break;
	case 18:
		return "Ar"; break;
	case 19:
		return "K"; break;
	case 20:
		return "Ca"; break;
	case 21:
		return "Sc"; break;
	case 22:
		return "Ti"; break;
	case 23:
		return "V"; break;
	case 24:
		return "Cr"; break;
	case 25:
		return "Mn"; break;
	case 26:
		return "Fe"; break;
	case 27:
		return "Co"; break;
	case 28:
		return "Ni"; break;
	case 29:
		return "Cu"; break;
	case 30:
		return "Zn"; break;
	case 31:
		return "Ga"; break;
	case 32:
		return "Ge"; break;
	case 33:
		return "As"; break;
	case 34:
		return "Se"; break;
	case 35:
		return "Br"; break;
	case 36:
		return "Kr"; break;
	case 37:
		return "Rb"; break;
	case 38:
		return "Sr"; break;
	case 39:
		return "Y"; break;
	case 40:
		return "Zr"; break;
	case 41:
		return "Nb"; break;
	case 42:
		return "Mo"; break;
	case 43:
		return "Tc"; break;
	case 44:
		return "Ru"; break;
	case 45:
		return "Rh"; break;
	case 46:
		return "Pd"; break;
	case 47:
		return "Ag"; break;
	case 48:
		return "Cd"; break;
	case 49:
		return "In"; break;
	case 50:
		return "Sn"; break;
	case 51:
		return "Sb"; break;
	case 52:
		return "Te"; break;
	case 53:
		return "I"; break;
	case 54:
		return "Xe"; break;
	case 55:
		return "Cs"; break;
	case 56:
		return "Ba"; break;
	case 57:
		return "La"; break;
	case 58:
		return "Ce"; break;
	case 59:
		return "Pr"; break;
	case 60:
		return "Nd"; break;
	case 61:
		return "Pm"; break;
	case 62:
		return "Sm"; break;
	case 63:
		return "Eu"; break;
	case 64:
		return "Gd"; break;
	case 65:
		return "Tb"; break;
	case 66:
		return "Dy"; break;
	case 67:
		return "Ho"; break;
	case 68:
		return "Er"; break;
	case 69:
		return "Tm"; break;
	case 70:
		return "Yb"; break;
	case 71:
		return "Lu"; break;
	case 72:
		return "Hf"; break;
	case 73:
		return "Ta"; break;
	case 74:
		return "W"; break;
	case 75:
		return "Re"; break;
	case 76:
		return "Os"; break;
	case 77:
		return "Ir"; break;
	case 78:
		return "Pt"; break;
	case 79:
		return "Au"; break;
	case 80:
		return "Hg"; break;
	case 81:
		return "Tl"; break;
	case 82:
		return "Pb"; break;
	case 83:
		return "Bi"; break;
	case 84:
		return "Po"; break;
	case 85:
		return "At"; break;
	case 86:
		return "Rn"; break;
	case 87:
		return "Fr"; break;
	case 88:
		return "Ra"; break;
	case 89:
		return "Ac"; break;
	case 90:
		return "Th"; break;
	case 91:
		return "Pa"; break;
	case 92:
		return "U"; break;
	case 93:
		return "Np"; break;
	case 94:
		return "Pu"; break;
	case 95:
		return "Am"; break;
	case 96:
		return "Cm"; break;
	case 97:
		return "Bk"; break;
	case 98:
		return "Cf"; break;
	case 99:
		return "Es"; break;
	case 100:
		return "Fm"; break;
	case 101:
		return "Md"; break;
	case 102:
		return "No"; break;
	case 103:
		return "Lr"; break;
	default:
		return "**"; 
		PrintErrorMsg("LoadElementName",VTEMLAB_EINVAL); 
		break;
	}
}

void AssignBFactors(int* atomType, double* bFactors, int vecLen)
{
	printf("> Enter B factors (293K) for %d elements\n", vecLen);
	for (int i = 0; i < vecLen; i++)
	{
		printf("%s (Z = %d): ", LoadElementName(atomType[i]), atomType[i]);
		scanf_s("%lf", &bFactors[i]);
	}
	getc(stdin);
}