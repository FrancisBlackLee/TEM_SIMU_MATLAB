function [probe] = ProbeCreate(params, xp, yp, Lx, Ly, Nx, Ny)
%PROBECREATE.M calculates the probe function with the input STEM
%parameters.
%   xp, yp --position of the center of the probe
%   Lx, Ly --side lengths
%   Nx, Ny --number of sampling

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

dx = Lx / Nx;
dy = Ly / Ny;
aberr_TF = AberrationFunction(params, Lx, Ly, Nx, Ny);
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[Fx, Fy] = meshgrid(fx, fy);
probe = ifft2(fftshift(aberr_TF .* exp(-1i * 2 * pi * (Fx * xp + Fy * yp)))) / (dx * dy);
probe = ifftshift(probe);
normCoeff = sqrt(sum(abs(probe.^2), 'all') * dx * dy);
probe = probe / normCoeff;

end

