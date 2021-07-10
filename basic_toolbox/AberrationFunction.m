function [otf] = AberrationFunction(params, Lx, Ly, Nx, Ny)
%ABERRATIONFUNCTION.M calculates the transfer function of the objective
%lens.
%   params --TEM&STEM parameters;
%   Lx, Ly, Nx, Ny -- sampling parameters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2021  Francis Black Lee and Li Xian

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cs = params.Cs * 1e7;
df = params.df;
wavLen = HighEnergyWavLen_X(params.KeV);
fx = InitFreqAxis(Lx, Nx);
fy = InitFreqAxis(Ly, Ny);
[Fx, Fy] = meshgrid(fx, fy);
freqSqr = Fx.^2 + Fy.^2;
aperture = CircApert_X(Lx, Ly, Nx, Ny, wavLen, params.amax);
otfPhase = pi * wavLen * freqSqr .* (0.5 * Cs * wavLen^2 * freqSqr - df);
otf = exp(-1i * otfPhase) .* aperture;

end

