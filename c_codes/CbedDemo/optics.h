/*
	optics.h is the open source part of vtemlab v0.0 engine,
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

#ifndef OPTICS_H
#define OPTICS_H

// mkl dependence
#include <mkl.h>

// Construct a Fresnel propagation kernel for given parameters:
int FresPropKer(MKL_Complex16* propKer, double wavLen, double propDist,
	double Lx, double Ly, int Nx, int Ny);

// Fresnel propagation:
// FFTSHIFTED_VEC: fftshifted wave;
// UNFFTSHIFTED_VEC: unfftshifted wave;
// FFTSHIFTED_VEC: fftshifted PropKer;
// UNFFTSHIFTED_VEC: unfftshifted PropKer;
int FresProp(MKL_Complex16* wave, int wMode, MKL_Complex16* propKer,
	int pkMode, int Nx, int Ny);

// Calculate the wave intensity
void GetWaveInten(double* waveI, MKL_Complex16* wave, int vecLen);

#endif // !OPTICS_H
