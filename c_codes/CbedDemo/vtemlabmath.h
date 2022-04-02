/*
    vtemlabmath.h is the open source part of vtemlab v0.0 engine,
    providing supplementary mathematical operations for vtemlab.

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

#ifndef VTEMLABMATH_H
#define VTEMLABMATH_H

// mkl dependence
#include <mkl.h>

// Generate random double:
double RandomDouble(double lower, double upper);

// Generate random int:
int RandomInt(int lower, int upper);

/*****************************************************************************/
// bessi0 and bessk0 are codes from EJK's computem
// Modified Bessel function I0(x)
double bessi0(double x);
// Modified Bessel function K0(x)
double bessk0(double x);
/*****************************************************************************/

// calculate the exchange of indices in one dimension after fftshift
int fftshift_Index_exchange_1D(int oldIdx, int dimLength); // DimLength -- dimension length
// calculate the exchange of indices in one dimension after ifftshift
int ifftshift_Index_exchange_1D(int oldIdx, int dimLength); // DimLength -- dimension length

/******************************************************************************************************/
// Indices exchange for 2D fftshift, and the 2D matrices are stored in a Ny * Nx sized row vector
int fftshift_Index_exchange_2D(int oldIdx, int Nx, int Ny);
// Indices exchange for 2D ifftshift, and the 2D matrices are stored in a Ny * Nx sized row vector
int ifftshift_Index_exchange_2D(int oldIdx, int Nx, int Ny);

/******************************************************************************************************/
// out-of-place 2D fftshift of single-precision real vector
void fftshift_S_2D_OP(float* unshiftedVec, float* shiftedVec, int Nx, int Ny);
// out-of-place 2D fftshift of double-precision real vector
void fftshift_D_2D_OP(double* unshiftedVec, double* shiftedVec, int Nx, int Ny);
// out-of-place 2D fftshift of single-precision complex vector
void fftshift_C_2D_OP(MKL_Complex8* unshiftedVec, MKL_Complex8* shiftedVec, int Nx, int Ny);
// out-of-place 2D fftshift of double-precision complex vector
void fftshift_Z_2D_OP(MKL_Complex16* unshiftedVec, MKL_Complex16* shiftedVec, int Nx, int Ny);

/******************************************************************************************************/
// out-of-place 2D ifftshift of single-precision real vector
void ifftshift_S_2D_OP(float* unshiftedVec, float* shiftedVec, int Nx, int Ny);
// out-of-place 2D ifftshift of double-precision real vector
void ifftshift_D_2D_OP(double* unshiftedVec, double* shiftedVec, int Nx, int Ny);
// out-of-place 2D ifftshift of single-precision complex vector
void ifftshift_C_2D_OP(MKL_Complex8* unshiftedVec, MKL_Complex8* shiftedVec, int Nx, int Ny);
// out-of-place 2D ifftshift of double-precision complex vector
void ifftshift_Z_2D_OP(MKL_Complex16* unshiftedVec, MKL_Complex16* shiftedVec, int Nx, int Ny);

/******************************************************************************************************/
// SBCF: supplementary basic calculation functions, GCD: greatest commom divisor
int gcd(int a, int b);

// left circular shift of both dimension (single-precision real vector):
void LeftCircShift_S_2d(float* mat, int sx, int sy, int Nx, int Ny);

// left circular shift of both dimension (double-precision real vector):
void LeftCircShift_D_2d(double* mat, int sx, int sy, int Nx, int Ny);

// left circular shift of both dimension (single-precision complex vector):
void LeftCircShift_C_2d(MKL_Complex8* mat, int sx, int sy, int Nx, int Ny);

// left circular shift of both dimension (double-precision complex vector):
void LeftCircShift_Z_2d(MKL_Complex16* mat, int sx, int sy, int Nx, int Ny);

// Right circular shift of both dimension (single-precision real vector):
void RightCircShift_S_2d(float* mat, int sx, int sy, int Nx, int Ny);

// Right circular shift of both dimension (double-precision real vector):
void RightCircShift_D_2d(double* mat, int sx, int sy, int Nx, int Ny);

// Right circular shift of both dimension (single-precision complex vector):
void RightCircShift_C_2d(MKL_Complex8* mat, int sx, int sy, int Nx, int Ny);

// Right circular shift of both dimension (double-precision complex vector):
void RightCircShift_Z_2d(MKL_Complex16* mat, int sx, int sy, int Nx, int Ny);

// matrix circular shift (single-precision real vector, sx or sy, positive for left shift, negative for right shift):
void CircShift_S_2d(float* mat, int sx, int sy, int Nx, int Ny);

// matrix circular shift (double-precision real vector, sx or sy, positive for left shift, negative for right shift):
void CircShift_D_2d(double* mat, int sx, int sy, int Nx, int Ny);

// matrix circular shift (single-precision complex vector, sx or sy, positive for left shift, negative for right shift):
void CircShift_C_2d(MKL_Complex8* mat, int sx, int sy, int Nx, int Ny);

// matrix circular shift (single-precision complex vector, sx or sy, positive for left shift, negative for right shift):
void CircShift_Z_2d(MKL_Complex16* mat, int sx, int sy, int Nx, int Ny);

// in-place 2D fftshift (single-precision real vector)
void fftshift_S_2D_IP(float* x, int Nx, int Ny);

// in-place 2D fftshift (double-precision real vector)
void fftshift_D_2D_IP(double* x, int Nx, int Ny);

// in-place 2D fftshift (single-precision complex vector)
void fftshift_C_2D_IP(MKL_Complex8* x, int Nx, int Ny);

// in-place 2D fftshift (double-precision complex vector)
void fftshift_Z_2D_IP(MKL_Complex16* x, int Nx, int Ny);

// in-place 2D ifftshift (single-precision real vector)
void ifftshift_S_2D_IP(float* x, int Nx, int Ny);

// in-place 2D ifftshift (double-precision real vector)
void ifftshift_D_2D_IP(double* x, int Nx, int Ny);

// in-place 2D ifftshift (single-precision complex vector)
void ifftshift_C_2D_IP(MKL_Complex8* x, int Nx, int Ny);

// in-place 2D ifftshift (double-precision complex vector)
void ifftshift_Z_2D_IP(MKL_Complex16* x, int Nx, int Ny);

#endif // !VTEMLABMATH_H
