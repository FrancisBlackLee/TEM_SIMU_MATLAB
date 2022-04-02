/*
	vtemlabio.cpp is the open source part of vtemlab v0.0 engine,
	providing file I/O operations for vtemlab.

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
#include <io.h>
#include <direct.h>

#include "vtemlabio.h"
#include "status.h"
#include "mode.h"
#include "mathconst.h"

// save single-precision vector:
int SaveSvec(int vecLen, float* x, char* filename)
{
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "w+");
	if (err || (vecLen == 0))
	{
		PrintErrorMsg("SaveSvec", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}
	int i;
	for (i = 0; i < vecLen - 1; i++)
	{
		fprintf_s(destFile, "%f\t", x[i]);
	}
	fprintf_s(destFile, "%f", x[i]);
	fclose(destFile);

	return VTEMLAB_SUCCESS;
}

// save double-precision vector:
int SaveDvec(int vecLen, double* x, char* filename)
{
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "w+");
	if (err || (vecLen == 0))
	{
		PrintErrorMsg("SaveDvec", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}
	int i;
	for (i = 0; i < vecLen - 1; i++)
		fprintf_s(destFile, "%lf\t", x[i]);
	fprintf_s(destFile, "%lf", x[i]);
	fclose(destFile);

	return VTEMLAB_SUCCESS;
}

// save single-precision complex vector, 1st line real, 2nd line imag:
int SaveCvec(int vecLen, MKL_Complex8* x, char* filename)
{
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "w+");
	if (err || (vecLen == 0))
	{
		PrintErrorMsg("SaveCvec", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}
	int i;
	for (i = 0; i < vecLen - 1; i++)
		fprintf_s(destFile, "%f\t", x[i].real);
	fprintf_s(destFile, "%f\n", x[i].real);
	for (i = 0; i < vecLen - 1; i++)
		fprintf_s(destFile, "%f\t", x[i].imag);
	fprintf_s(destFile, "%f", x[i].imag);
	fclose(destFile);
	printf("Vector saved as\n");
	puts(filename);
	printf("1st line real and 2nd line imaginary.\n\n");

	return VTEMLAB_SUCCESS;
}

// save double-precision complex vector, 1st line real, 2nd line imag:
int SaveZvec(int vecLen, MKL_Complex16* x, char* filename)
{
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "w+");
	if (err || (vecLen == 0))
	{
		PrintErrorMsg("SaveZvec", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}
	int i;
	for (i = 0; i < vecLen - 1; i++)
		fprintf_s(destFile, "%lf\t", x[i].real);
	fprintf_s(destFile, "%lf\n", x[i].real);
	for (i = 0; i < vecLen - 1; i++)
		fprintf_s(destFile, "%lf\t", x[i].imag);
	fprintf_s(destFile, "%lf", x[i].imag);
	fclose(destFile);
	printf("Vector saved as\n");
	puts(filename);
	printf("1st line real and 2nd line imaginary.\n\n");

	return VTEMLAB_SUCCESS;
}

// save single-precision vector as a matrix:
int SaveSvecAsMat(int Nx, int Ny, float* x, char* filename)
{
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "w+");
	if (err || (Ny * Nx == 0))
	{
		PrintErrorMsg("SaveSvecAsMat", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}
	int ix, iy;
	for (iy = 0; iy < Ny; iy++)
	{
		for (ix = 0; ix < Nx - 1; ix++)
			fprintf_s(destFile, "%f\t", x[iy * Nx + ix]);
		if (iy == Ny - 1)
			fprintf_s(destFile, "%f", x[iy * Nx + ix]);
		else
			fprintf_s(destFile, "%f\n", x[iy * Nx + ix]);
	}
	fclose(destFile);

	return VTEMLAB_SUCCESS;
}

//save double-precision vector as a matrix:
int SaveDvecAsMat(int Nx, int Ny, double* x, char* filename)
{
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "w+");
	if (err || (Ny * Nx == 0))
	{
		PrintErrorMsg("SaveDvecAsMat", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}
	int ix, iy;
	for (iy = 0; iy < Ny; iy++)
	{
		for (ix = 0; ix < Nx - 1; ix++)
			fprintf_s(destFile, "%lf\t", x[iy * Nx + ix]);
		if (iy == Ny - 1)
			fprintf_s(destFile, "%lf", x[iy * Nx + ix]);
		else
			fprintf_s(destFile, "%lf\n", x[iy * Nx + ix]);
	}
	fclose(destFile);

	return VTEMLAB_SUCCESS;
}

// save single-precision complex vector as a matrix:
//(Ny*Nx) sized array saved as 2Ny by Nx matrix, 1st Ny by Nx block is real, 
// 2nd Ny by Nx block is imag:
int SaveCvecAsMat(int Nx, int Ny, MKL_Complex8* x, char* filename)
{
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "w+");
	if (err || (Ny * Nx == 0))
	{
		PrintErrorMsg("SaveCvecAsMat", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}
	int ix, iy;
	// Saving real component:
	for (iy = 0; iy < Ny; iy++)
	{
		for (ix = 0; ix < Nx - 1; ix++)
			fprintf_s(destFile, "%f\t", x[iy * Nx + ix].real);
		fprintf_s(destFile, "%f\n", x[iy * Nx + ix].real);
	}
	// Saving imag component:
	for (iy = 0; iy < Ny; iy++)
	{
		for (ix = 0; ix < Nx - 1; ix++)
			fprintf_s(destFile, "%f\t", x[iy * Nx + ix].imag);
		if (iy == Ny - 1)
			fprintf_s(destFile, "%f", x[iy * Nx + ix].imag);
		else
			fprintf_s(destFile, "%f\n", x[iy * Nx + ix].imag);
	}
	fclose(destFile);
	printf("Vector saved as\n");
	puts(filename);
	printf("Real component saved as the 1st %d by %d block\n"
		"Imaginary component saved as the 2nd %d by %d block\n", 
		Ny, Nx, Ny, Nx);

	return VTEMLAB_SUCCESS;
}

// save double-precision complex vector as a matrix:
//(Ny*Nx) sized array saved as 2Ny by Nx matrix, 1st Ny by Nx block is real, 
// 2nd Ny by Nx block is imag:
int SaveZvecAsMat(int Nx, int Ny, MKL_Complex16* x, char* filename)
{
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "w+");
	if (err || (Ny * Nx == 0))
	{
		PrintErrorMsg("SaveZvecAsMat", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}
	int ix, iy;
	// Saving real component:
	for (iy = 0; iy < Ny; iy++)
	{
		for (ix = 0; ix < Nx - 1; ix++)
			fprintf_s(destFile, "%lf\t", x[iy * Nx + ix].real);
		fprintf_s(destFile, "%lf\n", x[iy * Nx + ix].real);
	}
	// Saving imag component:
	for (iy = 0; iy < Ny; iy++)
	{
		for (ix = 0; ix < Nx - 1; ix++)
			fprintf_s(destFile, "%lf\t", x[iy * Nx + ix].imag);
		if (iy == Ny - 1)
			fprintf_s(destFile, "%lf", x[iy * Nx + ix].imag);
		else
			fprintf_s(destFile, "%lf\n", x[iy * Nx + ix].imag);
	}
	fclose(destFile);
	printf("Vector saved as\n");
	puts(filename);
	printf("Real component saved as the 1st %d by %d block\n"
		"Imaginary component saved as the 2nd %d by %d block\n", 
		Ny, Nx, Ny, Nx);

	return VTEMLAB_SUCCESS;
}


// load integer vector;
int LoadIntVec(int vecLen, int* x, char* filename)
{
	if (x == NULL)
	{
		PrintErrorMsg("LoadIntVec", VTEMLAB_ENULLPTR);
		return VTEMLAB_ENULLPTR;
	}
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "r");
	if (err || (vecLen == 0))
	{
		PrintErrorMsg("LoadIntVec", VTEMLAB_EFILER);
		return VTEMLAB_EFILER;
	}

	for (int i = 0; i < vecLen; i++)
		fscanf_s(destFile, "%d", &x[i]);
	fclose(destFile);

	return VTEMLAB_SUCCESS;
}


// Load single-precision vector:
int LoadSvec(int vecLen, float* x, char* filename)
{
	if (x == NULL)
	{
		PrintErrorMsg("LoadSvec", VTEMLAB_ENULLPTR);
		return VTEMLAB_ENULLPTR;
	}
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "r");
	if (err || (vecLen == 0))
	{
		PrintErrorMsg("LoadSvec", VTEMLAB_EFILER);
		return VTEMLAB_EFILER;
	}

	for (int i = 0; i < vecLen; i++)
		fscanf_s(destFile, "%f", &x[i]);
	fclose(destFile);

	return VTEMLAB_SUCCESS;
}

// Load double-precision vector:
int LoadDvec(int vecLen, double* x, char* filename)
{
	if (x == NULL)
	{
		PrintErrorMsg("LoadDvec", VTEMLAB_ENULLPTR);
		return VTEMLAB_ENULLPTR;
	}
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "r");
	if (err || (vecLen == 0))
	{
		PrintErrorMsg("LoadDvec", VTEMLAB_EFILER);
		return VTEMLAB_EFILER;
	}

	for (int i = 0; i < vecLen; i++)
		fscanf_s(destFile, "%lf", &x[i]);
	fclose(destFile);

	return VTEMLAB_SUCCESS;
}

// Load single-precision complex vector:
// length = N, file contains 2N elements, 1st ~ Nth for real component; 
// (N+1)th ~ (2N)th for imag component.
int LoadCvec(int vecLen, MKL_Complex8* x, char* filename)
{
	if (x == NULL)
	{
		PrintErrorMsg("LoadCvec", VTEMLAB_ENULLPTR);
		return VTEMLAB_ENULLPTR;
	}
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "r");
	if (err || (vecLen == 0))
	{
		PrintErrorMsg("LoadCvec", VTEMLAB_EFILER);
		return VTEMLAB_EFILER;
	}

	for (int i = 0; i < vecLen; i++)
		fscanf_s(destFile, "%f", &x[i].real);
	for (int i = 0; i < vecLen; i++)
		fscanf_s(destFile, "%f", &x[i].imag);
	fclose(destFile);

	return VTEMLAB_SUCCESS;
}

// Load double-precision complex vector:
// length = N, file contains 2N elements, 1st ~ Nth for real component; 
// (N+1)th ~ (2N)th for imag component.
int LoadZvec(int vecLen, MKL_Complex16* x, char* filename)
{
	if (x == NULL)
	{
		PrintErrorMsg("LoadZvec", VTEMLAB_ENULLPTR);
		return VTEMLAB_ENULLPTR;
	}
	FILE* destFile;
	errno_t err;
	err = fopen_s(&destFile, filename, "r");
	if (err || (vecLen == 0))
	{
		PrintErrorMsg("LoadZvec", VTEMLAB_EFILER);
		return VTEMLAB_EFILER;
	}

	for (int i = 0; i < vecLen; i++)
		fscanf_s(destFile, "%lf", &x[i].real);
	for (int i = 0; i < vecLen; i++)
		fscanf_s(destFile, "%lf", &x[i].imag);
	fclose(destFile);

	return VTEMLAB_SUCCESS;
}


// save vector as binary file:
// single-precision floating-point vector
int SaveSvecAsBinaryFile(char* filename, float* x, int vecLen)
{
	FILE* destFile;
	if (fopen_s(&destFile, filename, "wb") == 0)
	{
		int count = 0;
		count += fwrite(x, sizeof(float), vecLen, destFile);
		fclose(destFile);
	}
	else
	{
		PrintErrorMsg("SaveSvecAsBinaryFile", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}

	return VTEMLAB_SUCCESS;
}


// double-precision floating-point vector
int SaveDvecAsBinaryFile(char* filename, double* x, int vecLen)
{
	FILE* destFile;
	if (fopen_s(&destFile, filename, "wb") == 0)
	{
		int count = 0;
		count += fwrite(x, sizeof(double), vecLen, destFile);
		fclose(destFile);
	}
	else
	{
		PrintErrorMsg("SaveDvecAsBinaryFile", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}

	return VTEMLAB_SUCCESS;
}

// single-precision floating-point complex vector
// every two saved values (real + imag) are one element in the complex vector.
int SaveCvecAsBinaryFile(char* filename, MKL_Complex8* x, int vecLen)
{
	FILE* destFile;
	if (fopen_s(&destFile, filename, "wb") == 0)
	{
		int count = 0;
		for (int i = 0; i < vecLen; i++)
		{
			count += fwrite(&x[i].real, sizeof(float), 1, destFile);
			count += fwrite(&x[i].imag, sizeof(float), 1, destFile);
		}
		fclose(destFile);
	}
	else
	{
		PrintErrorMsg("SaveCvecAsBinaryFile", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}

	return VTEMLAB_SUCCESS;
}

// double-precision floating-point complex vector
// every two saved values (real + imag) are one element in the complex vector.
int SaveZvecAsBinaryFile(char* filename, MKL_Complex16* x, int vecLen)
{
	FILE* destFile;
	if (fopen_s(&destFile, filename, "wb") == 0)
	{
		int count = 0;
		for (int i = 0; i < vecLen; i++)
		{
			count += fwrite(&x[i].real, sizeof(double), 1, destFile);
			count += fwrite(&x[i].imag, sizeof(double), 1, destFile);
		}
		fclose(destFile);
	}
	else
	{
		PrintErrorMsg("SaveZvecAsBinaryFile", VTEMLAB_EFILEW);
		return VTEMLAB_EFILEW;
	}

	return VTEMLAB_SUCCESS;
}


// Load vector from binary file:
// single-precision floating-point vector
int LoadSvecFromBinaryFile(char* filename, float* x, int vecLen)
{
	FILE* destFile;
	if (fopen_s(&destFile, filename, "rb") == 0)
	{
		int count = 0;
		count += fread(x, sizeof(float), vecLen, destFile);
		fclose(destFile);
	}
	else
	{
		PrintErrorMsg("LoadSvecFromBinaryFile", VTEMLAB_EFILER);
		return VTEMLAB_EFILER;
	}

	return VTEMLAB_SUCCESS;
}


// double-precision floating-point vector
int LoadDvecFromBinaryFile(char* filename, double* x, int vecLen)
{
	FILE* destFile;
	if (fopen_s(&destFile, filename, "rb") == 0)
	{
		int count = 0;
		count += fread(x, sizeof(double), vecLen, destFile);
		fclose(destFile);
	}
	else
	{
		PrintErrorMsg("LoadDvecFromBinaryFile", VTEMLAB_EFILER);
		return VTEMLAB_EFILER;
	}

	return VTEMLAB_SUCCESS;
}


// single-precision floating-point complex vector
// every two loaded values (real + imag) are one element in the complex vector.
int LoadCvecFromBinaryFile(char* filename, MKL_Complex8* x, int vecLen)
{
	FILE* destFile;
	if (fopen_s(&destFile, filename, "rb") == 0)
	{
		int count = 0;
		for (int i = 0; i < vecLen; i++)
		{
			count += fread(&x[i].real, sizeof(float), 1, destFile);
			count += fread(&x[i].imag, sizeof(float), 1, destFile);
		}
		fclose(destFile);
	}
	else
	{
		PrintErrorMsg("LoadCvecFromBinaryFile", VTEMLAB_EFILER);
		return VTEMLAB_EFILER;
	}

	return VTEMLAB_SUCCESS;
}


// double-precision floating-point complex vector
// every two loaded values (real + imag) are one element in the complex vector.
int LoadZvecFromBinaryFile(char* filename, MKL_Complex16* x, int vecLen)
{
	FILE* destFile;
	if (fopen_s(&destFile, filename, "rb") == 0)
	{
		int count = 0;
		for (int i = 0; i < vecLen; i++)
		{
			count += fread(&x[i].real, sizeof(double), 1, destFile);
			count += fread(&x[i].imag, sizeof(double), 1, destFile);
		}
		fclose(destFile);
	}
	else
	{
		PrintErrorMsg("LoadZvecFromBinaryFile", VTEMLAB_EFILER);
		return VTEMLAB_EFILER;
	}

	return VTEMLAB_SUCCESS;
}