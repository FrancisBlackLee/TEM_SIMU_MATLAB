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
function [ u2, Lx_2, Ly_2 ] = propFF_1( u1, Lx_1, Ly_1, waveLength, propDist )
% Annotations to be updated

[Ny, Nx] = size(u1);
dx_1 = Lx_1 / Nx;
dy_1 = Ly_1 / Ny;
k = 2 * pi / waveLength;

Lx_2 = waveLength * propDist / dx_1;
Ly_2 = waveLength * propDist / dy_1;
dx_2 = waveLength * propDist / Lx_1;
dy_2 = waveLength * propDist / Ly_1;
x_2 = -Lx_2 / 2 : dx_2 : Lx_2 / 2 - dx_2;
y_2 = -Ly_2 / 2 : dy_2 : Ly_2 / 2 - dy_2;
[xGrid_2, yGrid_2] = meshgrid(x_2, y_2);

kernel = 1 / (1i * waveLength * propDist) * ...
    exp(1i * k / (2 * propDist) * (xGrid_2.^2 + yGrid_2.^2));
u2 = kernel .* ifftshift(fft2(fftshift(u1))) * dx_1 * dx_2;

end

