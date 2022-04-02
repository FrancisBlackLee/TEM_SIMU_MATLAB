/*
    vtemlabmath.cpp is the open source part of vtemlab v0.0 engine,
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

#include <math.h>

#include "vtemlabmath.h"
#include "mathconst.h"
#include "mode.h"
#include "status.h"


// Generate random double:
double RandomDouble(double lower, double upper)
{
	return (double)rand() / RAND_MAX * (upper - lower) + lower;
}

// Generate random int:
int RandomInt(int lower, int upper)
{
	return rand() % (upper - lower + 1) + lower;
}


/*****************************************************************************/
// bessi0 and bessk0 are codes from EJK's computem
// Modified Bessel function I0(x)
double bessi0(double x)
{
	int i;
	double ax, sum, t;

	static double i0a[] = { 1.0, 3.5156229, 3.0899424, 1.2067492,
		0.2659732, 0.0360768, 0.0045813 };

	static double i0b[] = { 0.39894228, 0.01328592, 0.00225319,
		-0.00157565, 0.00916281, -0.02057706, 0.02635537,
		-0.01647633, 0.00392377 };

	ax = fabs(x);
	if (ax <= 3.75)
	{
		t = x / 3.75;
		t = t * t;
		sum = i0a[6];
		for (i = 5; i >= 0; i--)
			sum = sum * t + i0a[i];
	}
	else
	{
		t = 3.75 / ax;
		sum = i0b[8];
		for (i = 7; i >= 0; i--)
			sum = sum * t + i0b[i];
		sum = exp(ax) * sum / sqrt(ax);
	}
	return(sum);
}

// Modified Bessel function K0(x)
double bessk0(double x)
{
	double bessi0(double);

	int i;
	double ax, x2, sum;
	static double k0a[] = { -0.57721566, 0.42278420, 0.23069756,
		 0.03488590, 0.00262698, 0.00010750, 0.00000740 };

	static double k0b[] = { 1.25331414, -0.07832358, 0.02189568,
		 -0.01062446, 0.00587872, -0.00251540, 0.00053208 };

	ax = fabs(x);
	if ((ax > 0.0) && (ax <= 2.0))
	{
		x2 = ax / 2.0;
		x2 = x2 * x2;
		sum = k0a[6];
		for (i = 5; i >= 0; i--)
			sum = sum * x2 + k0a[i];
		sum = -log(ax / 2.0) * bessi0(x) + sum;
	}
	else if (ax > 2.0)
	{
		x2 = 2.0 / ax;
		sum = k0b[6];
		for (i = 5; i >= 0; i--)
			sum = sum * x2 + k0b[i];
		sum = exp(-ax) * sum / sqrt(ax);
	}
	else sum = 1.0e20;
	return (sum);
}
/*****************************************************************************/


// calculate the exchange of indices in one dimension after fftshift
int fftshift_Index_exchange_1D(int oldIdx, int dimLength) // DimLength -- dimension length
{
	return (oldIdx + dimLength / 2) % dimLength;
}

// calculate the exchange of indices in one dimension after ifftshift
int ifftshift_Index_exchange_1D(int oldIdx, int dimLength) // DimLength -- dimension length
{
	// return (OldIndex + (DimLength + 1) / 2) % DimLength;
	return (oldIdx - dimLength / 2 + dimLength) % dimLength;
}

// Indices exchange for 2D fftshift, and the 2D matrices are stored in a Ny * Nx sized row vector
int fftshift_Index_exchange_2D(int oldIdx, int Nx, int Ny)
{
	int xIdx, yIdx, xIdx_s, yIdx_s;// the suffix "_s" denotes shifted index
	xIdx = oldIdx % Nx;
	yIdx = oldIdx / Nx;
	xIdx_s = fftshift_Index_exchange_1D(xIdx, Nx);
	yIdx_s = fftshift_Index_exchange_1D(yIdx, Ny);
	return yIdx_s * Nx + xIdx_s;
}

// Indices exchange for 2D ifftshift, and the 2D matrices are stored in a Ny * Nx sized row vector
int ifftshift_Index_exchange_2D(int oldIdx, int Nx, int Ny)
{
	int xIdx, yIdx, xIdx_s, yIdx_s;// the suffix "_s" denotes shifted index
	xIdx = oldIdx % Nx;
	yIdx = oldIdx / Nx;
	xIdx_s = ifftshift_Index_exchange_1D(xIdx, Nx);
	yIdx_s = ifftshift_Index_exchange_1D(yIdx, Ny);
	return yIdx_s * Nx + xIdx_s;
}

// out-of-place 2D fftshift of single-precision real vector
void fftshift_S_2D_OP(float* unshiftedVec, float* shiftedVec, int Nx, int Ny)
{
	int oldIdx, newIdx;
	for (oldIdx = 0; oldIdx < Ny * Nx; oldIdx++)
	{
		newIdx = fftshift_Index_exchange_2D(oldIdx, Nx, Ny);
		shiftedVec[newIdx] = unshiftedVec[oldIdx];
	}
}

// out-of-place 2D fftshift of double-precision real vector
void fftshift_D_2D_OP(double* unshiftedVec, double* shiftedVec, int Nx, int Ny)
{
	int oldIdx, newIdx;
	for (oldIdx = 0; oldIdx < Ny * Nx; oldIdx++)
	{
		newIdx = fftshift_Index_exchange_2D(oldIdx, Nx, Ny);
		shiftedVec[newIdx] = unshiftedVec[oldIdx];
	}
}

// out-of-place 2D fftshift of single-precision complex vector
void fftshift_C_2D_OP(MKL_Complex8* unshiftedVec, MKL_Complex8* shiftedVec, int Nx, int Ny)
{
	int oldIdx, newIdx;
	for (oldIdx = 0; oldIdx < Ny * Nx; oldIdx++)
	{
		newIdx = fftshift_Index_exchange_2D(oldIdx, Nx, Ny);
		shiftedVec[newIdx] = unshiftedVec[oldIdx];
	}
}

// out-of-place 2D fftshift of double-precision complex vector
void fftshift_Z_2D_OP(MKL_Complex16* unshiftedVec, MKL_Complex16* shiftedVec, int Nx, int Ny)
{
	int oldIdx, newIdx;
	for (oldIdx = 0; oldIdx < Ny * Nx; oldIdx++)
	{
		newIdx = fftshift_Index_exchange_2D(oldIdx, Nx, Ny);
		shiftedVec[newIdx] = unshiftedVec[oldIdx];
	}
}

// out-of-place 2D ifftshift of single-precision real vector
void ifftshift_S_2D_OP(float* unshiftedVec, float* shiftedVec, int Nx, int Ny)
{
	int oldIdx, newIdx;
	for (oldIdx = 0; oldIdx < Ny * Nx; oldIdx++)
	{
		newIdx = ifftshift_Index_exchange_2D(oldIdx, Nx, Ny);
		shiftedVec[newIdx] = unshiftedVec[oldIdx];
	}
}

// out-of-place 2D ifftshift of double-precision real vector
void ifftshift_D_2D_OP(double* unshiftedVec, double* shiftedVec, int Nx, int Ny)
{
	int oldIdx, newIdx;
	for (oldIdx = 0; oldIdx < Ny * Nx; oldIdx++)
	{
		newIdx = ifftshift_Index_exchange_2D(oldIdx, Nx, Ny);
		shiftedVec[newIdx] = unshiftedVec[oldIdx];
	}
}

// out-of-place 2D ifftshift of single-precision complex vector
void ifftshift_C_2D_OP(MKL_Complex8* unshiftedVec, MKL_Complex8* shiftedVec, int Nx, int Ny)
{
	int oldIdx, newIdx;
	for (oldIdx = 0; oldIdx < Ny * Nx; oldIdx++)
	{
		newIdx = ifftshift_Index_exchange_2D(oldIdx, Nx, Ny);
		shiftedVec[newIdx] = unshiftedVec[oldIdx];
	}
}

// out-of-place 2D ifftshift of double-precision complex vector
void ifftshift_Z_2D_OP(MKL_Complex16* unshiftedVec, MKL_Complex16* shiftedVec, int Nx, int Ny)
{
	int oldIdx, newIdx;
	for (oldIdx = 0; oldIdx < Ny * Nx; oldIdx++)
	{
		newIdx = ifftshift_Index_exchange_2D(oldIdx, Nx, Ny);
		shiftedVec[newIdx] = unshiftedVec[oldIdx];
	}
}

// GCD: greatest commom divisor
int gcd(int a, int b)
{
	if (b == 0)
		return a;
	else
		return gcd(b, a % b);
}

// left circular shift of both dimension (single-precision real vector):
void LeftCircShift_S_2d(float* mat, int sx, int sy, int Nx, int Ny)
{
	int ix, iy, jx, jy, kx, ky, gcdx, gcdy;
	float temp;
	gcdx = gcd(sx, Nx);
	gcdy = gcd(sy, Ny);
	for (ix = 0; ix < gcdx; ix++)
	{
		for (iy = 0; iy < Ny; iy++)
		{
			temp = mat[iy * Nx + ix];
			jx = ix;
			while (1)
			{
				kx = (jx + sx) % Nx;
				if (kx == ix)
					break;
				mat[iy * Nx + jx] = mat[iy * Nx + kx];
				jx = kx;
			}
			mat[iy * Nx + jx] = temp;
		}
	}
	for (iy = 0; iy < gcdy; iy++)
	{
		for (ix = 0; ix < Nx; ix++)
		{
			temp = mat[iy * Nx + ix];
			jy = iy;
			while (1)
			{
				ky = (jy + sy) % Ny;
				if (ky == iy)
					break;
				mat[jy * Nx + ix] = mat[ky * Nx + ix];
				jy = ky;
			}
			mat[jy * Nx + ix] = temp;
		}
	}
}

// left circular shift of both dimension (double-precision real vector):
void LeftCircShift_D_2d(double* mat, int sx, int sy, int Nx, int Ny)
{
	int ix, iy, jx, jy, kx, ky, gcdx, gcdy;
	double temp;
	gcdx = gcd(sx, Nx);
	gcdy = gcd(sy, Ny);
	for (ix = 0; ix < gcdx; ix++)
	{
		for (iy = 0; iy < Ny; iy++)
		{
			temp = mat[iy * Nx + ix];
			jx = ix;
			while (1)
			{
				kx = (jx + sx) % Nx;
				if (kx == ix)
					break;
				mat[iy * Nx + jx] = mat[iy * Nx + kx];
				jx = kx;
			}
			mat[iy * Nx + jx] = temp;
		}
	}
	for (iy = 0; iy < gcdy; iy++)
	{
		for (ix = 0; ix < Nx; ix++)
		{
			temp = mat[iy * Nx + ix];
			jy = iy;
			while (1)
			{
				ky = (jy + sy) % Ny;
				if (ky == iy)
					break;
				mat[jy * Nx + ix] = mat[ky * Nx + ix];
				jy = ky;
			}
			mat[jy * Nx + ix] = temp;
		}
	}
}

// left circular shift of both dimension (single-precision complex vector):
void LeftCircShift_C_2d(MKL_Complex8* mat, int sx, int sy, int Nx, int Ny)
{
	int ix, iy, jx, jy, kx, ky, gcdx, gcdy;
	MKL_Complex8 temp;
	gcdx = gcd(sx, Nx);
	gcdy = gcd(sy, Ny);
	for (ix = 0; ix < gcdx; ix++)
	{
		for (iy = 0; iy < Ny; iy++)
		{
			temp = mat[iy * Nx + ix];
			jx = ix;
			while (1)
			{
				kx = (jx + sx) % Nx;
				if (kx == ix)
					break;
				mat[iy * Nx + jx] = mat[iy * Nx + kx];
				jx = kx;
			}
			mat[iy * Nx + jx] = temp;
		}
	}
	for (iy = 0; iy < gcdy; iy++)
	{
		for (ix = 0; ix < Nx; ix++)
		{
			temp = mat[iy * Nx + ix];
			jy = iy;
			while (1)
			{
				ky = (jy + sy) % Ny;
				if (ky == iy)
					break;
				mat[jy * Nx + ix] = mat[ky * Nx + ix];
				jy = ky;
			}
			mat[jy * Nx + ix] = temp;
		}
	}
}

// left circular shift of both dimension (double-precision complex vector):
void LeftCircShift_Z_2d(MKL_Complex16* mat, int sx, int sy, int Nx, int Ny)
{
	int ix, iy, jx, jy, kx, ky, gcdx, gcdy;
	MKL_Complex16 temp;
	gcdx = gcd(sx, Nx);
	gcdy = gcd(sy, Ny);
	for (ix = 0; ix < gcdx; ix++)
	{
		for (iy = 0; iy < Ny; iy++)
		{
			temp = mat[iy * Nx + ix];
			jx = ix;
			while (1)
			{
				kx = (jx + sx) % Nx;
				if (kx == ix)
					break;
				mat[iy * Nx + jx] = mat[iy * Nx + kx];
				jx = kx;
			}
			mat[iy * Nx + jx] = temp;
		}
	}
	for (iy = 0; iy < gcdy; iy++)
	{
		for (ix = 0; ix < Nx; ix++)
		{
			temp = mat[iy * Nx + ix];
			jy = iy;
			while (1)
			{
				ky = (jy + sy) % Ny;
				if (ky == iy)
					break;
				mat[jy * Nx + ix] = mat[ky * Nx + ix];
				jy = ky;
			}
			mat[jy * Nx + ix] = temp;
		}
	}
}

// Right circular shift of both dimension (single-precision real vector):
void RightCircShift_S_2d(float* mat, int sx, int sy, int Nx, int Ny)
{
	int ix, iy, jx, jy, kx, ky, gcdx, gcdy;
	float temp;
	gcdx = gcd(sx, Nx);
	gcdy = gcd(sy, Ny);
	for (ix = 0; ix < gcdx; ix++)
	{
		for (iy = 0; iy < Ny; iy++)
		{
			temp = mat[iy * Nx + ix];
			jx = ix;
			while (1)
			{
				kx = (jx - sx + Nx) % Nx;
				if (kx == ix)
					break;
				mat[iy * Nx + jx] = mat[iy * Nx + kx];
				jx = kx;
			}
			mat[iy * Nx + jx] = temp;
		}
	}
	for (iy = 0; iy < gcdy; iy++)
	{
		for (ix = 0; ix < Nx; ix++)
		{
			temp = mat[iy * Nx + ix];
			jy = iy;
			while (1)
			{
				ky = (jy - sy + Ny) % Ny;
				if (ky == iy)
					break;
				mat[jy * Nx + ix] = mat[ky * Nx + ix];
				jy = ky;
			}
			mat[jy * Nx + ix] = temp;
		}
	}
}

// right circular shift of both dimension (double-precision real vector):
void RightCircShift_D_2d(double* mat, int sx, int sy, int Nx, int Ny)
{
	int ix, iy, jx, jy, kx, ky, gcdx, gcdy;
	double temp;
	gcdx = gcd(sx, Nx);
	gcdy = gcd(sy, Ny);
	for (ix = 0; ix < gcdx; ix++)
	{
		for (iy = 0; iy < Ny; iy++)
		{
			temp = mat[iy * Nx + ix];
			jx = ix;
			while (1)
			{
				kx = (jx - sx + Nx) % Nx;
				if (kx == ix)
					break;
				mat[iy * Nx + jx] = mat[iy * Nx + kx];
				jx = kx;
			}
			mat[iy * Nx + jx] = temp;
		}
	}
	for (iy = 0; iy < gcdy; iy++)
	{
		for (ix = 0; ix < Nx; ix++)
		{
			temp = mat[iy * Nx + ix];
			jy = iy;
			while (1)
			{
				ky = (jy - sy + Ny) % Ny;
				if (ky == iy)
					break;
				mat[jy * Nx + ix] = mat[ky * Nx + ix];
				jy = ky;
			}
			mat[jy * Nx + ix] = temp;
		}
	}
}

// right circular shift of both dimension (single-precision complex vector):
void RightCircShift_C_2d(MKL_Complex8* mat, int sx, int sy, int Nx, int Ny)
{
	int ix, iy, jx, jy, kx, ky, gcdx, gcdy;
	MKL_Complex8 temp;
	gcdx = gcd(sx, Nx);
	gcdy = gcd(sy, Ny);
	for (ix = 0; ix < gcdx; ix++)
	{
		for (iy = 0; iy < Ny; iy++)
		{
			temp = mat[iy * Nx + ix];
			jx = ix;
			while (1)
			{
				kx = (jx - sx + Nx) % Nx;
				if (kx == ix)
					break;
				mat[iy * Nx + jx] = mat[iy * Nx + kx];
				jx = kx;
			}
			mat[iy * Nx + jx] = temp;
		}
	}
	for (iy = 0; iy < gcdy; iy++)
	{
		for (ix = 0; ix < Nx; ix++)
		{
			temp = mat[iy * Nx + ix];
			jy = iy;
			while (1)
			{
				ky = (jy - sy + Ny) % Ny;
				if (ky == iy)
					break;
				mat[jy * Nx + ix] = mat[ky * Nx + ix];
				jy = ky;
			}
			mat[jy * Nx + ix] = temp;
		}
	}
}

// right circular shift of both dimension (double-precision complex vector):
void RightCircShift_Z_2d(MKL_Complex16* mat, int sx, int sy, int Nx, int Ny)
{
	int ix, iy, jx, jy, kx, ky, gcdx, gcdy;
	MKL_Complex16 temp;
	gcdx = gcd(sx, Nx);
	gcdy = gcd(sy, Ny);
	for (ix = 0; ix < gcdx; ix++)
	{
		for (iy = 0; iy < Ny; iy++)
		{
			temp = mat[iy * Nx + ix];
			jx = ix;
			while (1)
			{
				kx = (jx - sx + Nx) % Nx;
				if (kx == ix)
					break;
				mat[iy * Nx + jx] = mat[iy * Nx + kx];
				jx = kx;
			}
			mat[iy * Nx + jx] = temp;
		}
	}
	for (iy = 0; iy < gcdy; iy++)
	{
		for (ix = 0; ix < Nx; ix++)
		{
			temp = mat[iy * Nx + ix];
			jy = iy;
			while (1)
			{
				ky = (jy - sy + Ny) % Ny;
				if (ky == iy)
					break;
				mat[jy * Nx + ix] = mat[ky * Nx + ix];
				jy = ky;
			}
			mat[jy * Nx + ix] = temp;
		}
	}
}

// matrix circular shift (single-precision real vector, sx or sy, positive for left shift, negative for right shift):
void CircShift_S_2d(float* mat, int sx, int sy, int Nx, int Ny)
{
	int ix, iy, jx, jy, kx, ky, gcdx, gcdy, stepx, stepy;
	float temp;
	gcdx = (sx >= 0) ? gcd(sx, Nx) : gcd(-sx, Nx);
	gcdy = (sy >= 0) ? gcd(sy, Ny) : gcd(-sy, Ny);
	stepx = (sx >= 0) ? sx : (sx + Nx);
	stepy = (sy >= 0) ? sy : (sy + Ny);
	for (ix = 0; ix < gcdx; ix++)
	{
		for (iy = 0; iy < Ny; iy++)
		{
			temp = mat[iy * Nx + ix];
			jx = ix;
			while (1)
			{
				kx = (jx + stepx) % Nx;
				if (kx == ix)
					break;
				mat[iy * Nx + jx] = mat[iy * Nx + kx];
				jx = kx;
			}
			mat[iy * Nx + jx] = temp;
		}
	}
	for (iy = 0; iy < gcdy; iy++)
	{
		for (ix = 0; ix < Nx; ix++)
		{
			temp = mat[iy * Nx + ix];
			jy = iy;
			while (1)
			{
				ky = (jy + stepy) % Ny;
				if (ky == iy)
					break;
				mat[jy * Nx + ix] = mat[ky * Nx + ix];
				jy = ky;
			}
			mat[jy * Nx + ix] = temp;
		}
	}
}

// matrix circular shift (double-precision real vector, sx or sy, positive for left shift, negative for right shift):
void CircShift_D_2d(double* mat, int sx, int sy, int Nx, int Ny)
{
	int ix, iy, jx, jy, kx, ky, gcdx, gcdy, stepx, stepy;
	double temp;
	gcdx = (sx >= 0) ? gcd(sx, Nx) : gcd(-sx, Nx);
	gcdy = (sy >= 0) ? gcd(sy, Ny) : gcd(-sy, Ny);
	stepx = (sx >= 0) ? sx : (sx + Nx);
	stepy = (sy >= 0) ? sy : (sy + Ny);
	for (ix = 0; ix < gcdx; ix++)
	{
		for (iy = 0; iy < Ny; iy++)
		{
			temp = mat[iy * Nx + ix];
			jx = ix;
			while (1)
			{
				kx = (jx + stepx) % Nx;
				if (kx == ix)
					break;
				mat[iy * Nx + jx] = mat[iy * Nx + kx];
				jx = kx;
			}
			mat[iy * Nx + jx] = temp;
		}
	}
	for (iy = 0; iy < gcdy; iy++)
	{
		for (ix = 0; ix < Nx; ix++)
		{
			temp = mat[iy * Nx + ix];
			jy = iy;
			while (1)
			{
				ky = (jy + stepy) % Ny;
				if (ky == iy)
					break;
				mat[jy * Nx + ix] = mat[ky * Nx + ix];
				jy = ky;
			}
			mat[jy * Nx + ix] = temp;
		}
	}
}

// matrix circular shift (single-precision complex vector, sx or sy, positive for left shift, negative for right shift):
void CircShift_C_2d(MKL_Complex8* mat, int sx, int sy, int Nx, int Ny)
{
	int ix, iy, jx, jy, kx, ky, gcdx, gcdy, stepx, stepy;
	MKL_Complex8 temp;
	gcdx = (sx >= 0) ? gcd(sx, Nx) : gcd(-sx, Nx);
	gcdy = (sy >= 0) ? gcd(sy, Ny) : gcd(-sy, Ny);
	stepx = (sx >= 0) ? sx : (sx + Nx);
	stepy = (sy >= 0) ? sy : (sy + Ny);
	for (ix = 0; ix < gcdx; ix++)
	{
		for (iy = 0; iy < Ny; iy++)
		{
			temp = mat[iy * Nx + ix];
			jx = ix;
			while (1)
			{
				kx = (jx + stepx) % Nx;
				if (kx == ix)
					break;
				mat[iy * Nx + jx] = mat[iy * Nx + kx];
				jx = kx;
			}
			mat[iy * Nx + jx] = temp;
		}
	}
	for (iy = 0; iy < gcdy; iy++)
	{
		for (ix = 0; ix < Nx; ix++)
		{
			temp = mat[iy * Nx + ix];
			jy = iy;
			while (1)
			{
				ky = (jy + stepy) % Ny;
				if (ky == iy)
					break;
				mat[jy * Nx + ix] = mat[ky * Nx + ix];
				jy = ky;
			}
			mat[jy * Nx + ix] = temp;
		}
	}
}

// matrix circular shift (single-precision complex vector, sx or sy, positive for left shift, negative for right shift):
void CircShift_Z_2d(MKL_Complex16* mat, int sx, int sy, int Nx, int Ny)
{
	int ix, iy, jx, jy, kx, ky, gcdx, gcdy, stepx, stepy;
	MKL_Complex16 temp;
	gcdx = (sx >= 0) ? gcd(sx, Nx) : gcd(-sx, Nx);
	gcdy = (sy >= 0) ? gcd(sy, Ny) : gcd(-sy, Ny);
	stepx = (sx >= 0) ? sx : (sx + Nx);
	stepy = (sy >= 0) ? sy : (sy + Ny);
	for (ix = 0; ix < gcdx; ix++)
	{
		for (iy = 0; iy < Ny; iy++)
		{
			temp = mat[iy * Nx + ix];
			jx = ix;
			while (1)
			{
				kx = (jx + stepx) % Nx;
				if (kx == ix)
					break;
				mat[iy * Nx + jx] = mat[iy * Nx + kx];
				jx = kx;
			}
			mat[iy * Nx + jx] = temp;
		}
	}
	for (iy = 0; iy < gcdy; iy++)
	{
		for (ix = 0; ix < Nx; ix++)
		{
			temp = mat[iy * Nx + ix];
			jy = iy;
			while (1)
			{
				ky = (jy + stepy) % Ny;
				if (ky == iy)
					break;
				mat[jy * Nx + ix] = mat[ky * Nx + ix];
				jy = ky;
			}
			mat[jy * Nx + ix] = temp;
		}
	}
}

// in-place 2D fftshift (single-precision real vector)
void fftshift_S_2D_IP(float* x, int Nx, int Ny)
{
	RightCircShift_S_2d(x, Nx / 2, Ny / 2, Nx, Ny);
}

// in-place 2D fftshift (double-precision real vector)
void fftshift_D_2D_IP(double* x, int Nx, int Ny)
{
	RightCircShift_D_2d(x, Nx / 2, Ny / 2, Nx, Ny);
}

// in-place 2D fftshift (single-precision complex vector)
void fftshift_C_2D_IP(MKL_Complex8* x, int Nx, int Ny)
{
	RightCircShift_C_2d(x, Nx / 2, Ny / 2, Nx, Ny);
}

// in-place 2D fftshift (double-precision complex vector)
void fftshift_Z_2D_IP(MKL_Complex16* x, int Nx, int Ny)
{
	RightCircShift_Z_2d(x, Nx / 2, Ny / 2, Nx, Ny);
}

// in-place 2D ifftshift (single-precision real vector)
void ifftshift_S_2D_IP(float* x, int Nx, int Ny)
{
	LeftCircShift_S_2d(x, Nx / 2, Ny / 2, Nx, Ny);
}

// in-place 2D ifftshift (double-precision real vector)
void ifftshift_D_2D_IP(double* x, int Nx, int Ny)
{
	LeftCircShift_D_2d(x, Nx / 2, Ny / 2, Nx, Ny);
}

// in-place 2D ifftshift (single-precision complex vector)
void ifftshift_C_2D_IP(MKL_Complex8* x, int Nx, int Ny)
{
	LeftCircShift_C_2d(x, Nx / 2, Ny / 2, Nx, Ny);
}

// in-place 2D ifftshift (double-precision complex vector)
void ifftshift_Z_2D_IP(MKL_Complex16* x, int Nx, int Ny)
{
	LeftCircShift_Z_2d(x, Nx / 2, Ny / 2, Nx, Ny);
}