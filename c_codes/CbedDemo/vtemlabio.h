/*
	vtemlabio.h is the open source part of vtemlab v0.0 engine,
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

#ifndef VTEMLABIO_H
#define VTEMLABIO_H

// mkl dependence
#include <mkl.h>

// save single-precision vector:
int SaveSvec(int vecLen, float* x, char* filename);

// save double-precision vector:
int SaveDvec(int vecLen, double* x, char* filename);

// save single-precision complex vector, 1st line real, 2nd line imag:
int SaveCvec(int vecLen, MKL_Complex8* x, char* filename);

// save double-precision complex vector, 1st line real, 2nd line imag:
int SaveZvec(int vecLen, MKL_Complex16* x, char* filename);

// save single-precision vector as a matrix:
int SaveSvecAsMat(int Nx, int Ny, float* x, char* filename);

//save double-precision vector as a matrix:
int SaveDvecAsMat(int Nx, int Ny, double* x, char* filename);

// save single-precision complex vector as a matrix:
// (Ny*Nx) sized array saved as 2Ny by Nx matrix, 1st Ny by Nx block is real, 
// 2nd Ny by Nx block is imag:
int SaveCvecAsMat(int Nx, int Ny, MKL_Complex8* x, char* filename);

// save double-precision complex vector as a matrix:
// (Ny*Nx) sized array saved as 2Ny by Nx matrix, 1st Ny by Nx block is real, 
// 2nd Ny by Nx block is imag:
int SaveZvecAsMat(int Nx, int Ny, MKL_Complex16* x, char* filename);

// load integer vector;
int LoadIntVec(int vecLen, int* x, char* filename);

// Load single-precision vector:
int LoadSvec(int vecLen, float* x, char* filename);

// Load double-precision vector:
int LoadDvec(int vecLen, double* x, char* filename);

// Load single-precision complex vector:
// length = N, file contains 2N elements, 1st ~ Nth for real component; 
// (N+1)th ~ (2N)th for imag component.
int LoadCvec(int vecLen, MKL_Complex8* x, char* filename);

// Load double-precision complex vector:
// length = N, file contains 2N elements, 1st ~ Nth for real component; 
// (N+1)th ~ (2N)th for imag component.
int LoadZvec(int vecLen, MKL_Complex16* x, char* filename);

// save vector as binary file:
// single-precision floating-point vector
int SaveSvecAsBinaryFile(char* filename, float* x, int vecLen);

// double-precision floating-point vector
int SaveDvecAsBinaryFile(char* filename, double* x, int vecLen);

// single-precision floating-point complex vector
// every two saved values (real + imag) are one element in the complex vector.
int SaveCvecAsBinaryFile(char* filename, MKL_Complex8* x, int vecLen);

// double-precision floating-point complex vector
// every two saved values (real + imag) are one element in the complex vector.
int SaveZvecAsBinaryFile(char* filename, MKL_Complex16* x, int vecLen);

// Load vector from binary file:
// single-precision floating-point vector
int LoadSvecFromBinaryFile(char* filename, float* x, int vecLen);

// double-precision floating-point vector
int LoadDvecFromBinaryFile(char* filename, double* x, int vecLen);

// single-precision floating-point complex vector
// every two loaded values (real + imag) are one element in the complex vector.
int LoadCvecFromBinaryFile(char* filename, MKL_Complex8* x, int vecLen);

// double-precision floating-point complex vector
// every two loaded values (real + imag) are one element in the complex vector.
int LoadZvecFromBinaryFile(char* filename, MKL_Complex16* x, int vecLen);

#endif // !VTEMLABIO_H
