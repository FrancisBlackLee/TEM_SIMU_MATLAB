%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019  Francis Black Lee and Li Xian

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
function [Aberr_TF] = AberrationFunction(Params,Lx,Ly,Nx,Ny)
%ABERRATIONFUNCTION.M calculates the transfer function of the objective
%lens.
%   Params --TEM&STEM parameters

dx = Lx / Nx;
dy = Ly / Ny;
Cs = Params.Cs * 1e7;
df = Params.df;
KeV = Params.KeV;
AngleMax = Params.amax * 0.001;
lambda = 12.3986 / sqrt((2 * 511.0 + KeV) * KeV);  %wavelength
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[Fx, Fy] = meshgrid(fx, fy);
FreqSquare = Fx.^2 + Fy.^2;
Aperture = ones(size(Fx));
Aperture(FreqSquare >= (sin(AngleMax) / lambda)^2) = 0;
Chi = pi * lambda * FreqSquare .* (0.5 * Cs * lambda^2 * FreqSquare - df * ones(size(Fx)));
Aberr_TF = exp(-1i * Chi) .* Aperture;

end

