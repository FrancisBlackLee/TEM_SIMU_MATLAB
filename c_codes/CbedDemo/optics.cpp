/*
	optics.cpp is the open source part of vtemlab v0.0 engine,
	providing Fourier optics operations for vtemlab.

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

#include "mathconst.h"
#include "mode.h"
#include "optics.h"
#include "status.h"
#include "vtemlabmath.h"

// Construct a Fresnel propagation kernel for given parameters:
int FresPropKer(MKL_Complex16* propKer, double wavLen, double propDist, 
	double Lx, double Ly, int Nx, int Ny)
{
	if (propKer == NULL)
	{
		PrintErrorMsg("FresPropKer", VTEMLAB_ENULLPTR);
		return VTEMLAB_ENULLPTR;
	}

	double dx, dy, fx, fy, fxStart, fyStart, dfx, dfy, KerFac;
	dx = Lx / (double)Nx;
	dy = Ly / (double)Ny;
	fxStart = -1.0 / (2.0 * dx);
	fyStart = -1.0 / (2.0 * dy);
	dfx = 1.0 / Lx;
	dfy = 1.0 / Ly;

	int ix, iy;
	for (iy = 0; iy < Ny; iy++)
	{
		fy = fyStart + (double)iy * dfy;
		for (ix = 0; ix < Nx; ix++)
		{
			fx = fxStart + (double)ix * dfx;
			KerFac = -MathPI * wavLen * propDist * (fx * fx + fy * fy);
			propKer[iy * Nx + ix] = { cos(KerFac),sin(KerFac) };
		}
	}

	return VTEMLAB_SUCCESS;
}

// Fresnel propagation:
// FFTSHIFTED_VEC: fftshifted wave;
// UNFFTSHIFTED_VEC: unfftshifted wave;
// FFTSHIFTED_VEC: fftshifted PropKer;
// UNFFTSHIFTED_VEC: unfftshifted PropKer;
int FresProp(MKL_Complex16* wave, int wMode, MKL_Complex16* propKer, 
	int pkMode, int Nx, int Ny)
{
	if (wMode == UNFFTSHIFTED_VEC)
		fftshift_Z_2D_IP(wave, Nx, Ny);
	if (pkMode == UNFFTSHIFTED_VEC)
		fftshift_Z_2D_IP(propKer, Nx, Ny);

	// Constructing fft descriptor:
	DFTI_DESCRIPTOR_HANDLE fftDescHandle;
	MKL_LONG fftStat, matSize[2] = { (MKL_LONG)Ny,(MKL_LONG)Nx };
	fftStat = DftiCreateDescriptor(&fftDescHandle, DFTI_DOUBLE, DFTI_COMPLEX,
		2, matSize);
	fftStat = DftiCommitDescriptor(fftDescHandle);

	int vecLen = Ny * Nx;
	fftStat = DftiComputeForward(fftDescHandle, wave);
	vzMul(vecLen, propKer, wave, wave);
	fftStat = DftiComputeBackward(fftDescHandle, wave);
	cblas_zdscal(vecLen, 1 / (double)vecLen, wave, 1);
	fftStat = DftiFreeDescriptor(&fftDescHandle);

	if (wMode == UNFFTSHIFTED_VEC)
		ifftshift_Z_2D_IP(wave, Nx, Ny);
	if (pkMode == UNFFTSHIFTED_VEC)
		ifftshift_Z_2D_IP(propKer, Nx, Ny);

	return VTEMLAB_SUCCESS;
}

// Calculate the wave intensity
void GetWaveInten(double* waveI, MKL_Complex16* wave, int vecLen)
{
	vzAbs(vecLen, wave, waveI);
	vdSqr(vecLen, waveI, waveI);
}