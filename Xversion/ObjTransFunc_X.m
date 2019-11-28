%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019  Francis Black Lee and Li Xian

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outllok.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OTF] = ObjTransFunc_X(params, Lx, Ly, Nx, Ny)
%ObjTransFunc_X.m generates the objective transfer function in reciprocal
%space, only including Cs3 and Cs5, for more aberrations please use
%MultiAberrPhaseError_X.m.
%   Lx, Ly, Nx, Ny -- sampling parameters, L denotes side length and N the
%       sampling number in real space;
%   params -- TEM & STEM parameters: KeV, df, Cs3 and Cs5, NA is not
%       required, for the numerical aperture is generated outside:
%       KeV: beam energy in KeV;
%       df: defocus in angstrom (a negative defocus is used to eliminate
%           the effect of the spherical aberration);
%       Cs3, Cs5: 3rd and 5th order spherical aberration in millimeter;
% Note: X denotes an experimental version!

KeV = params.KeV;
Cs3 = params.Cs3 * 1e7;
Cs5 = params.Cs5 * 1e7;
df = params.df;

WavLen = 12.3986 / sqrt((2 * 511.0 + KeV) * KeV);  %wavelength

dx = Lx / Nx;
dy = Ly / Ny;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[FX, FY] = meshgrid(fx, fy);
AngFreqSqu = (FX.^2 + FY.^2) * WavLen^2; % squared angular frequency

PhaseError = 2*pi/WavLen * (0.5 * df * AngFreqSqu + 0.25 * Cs3 * AngFreqSqu.^2 ...
                            + 1/6 * Cs5 * AngFreqSqu.^3);
OTF = exp(1i * PhaseError);

end

