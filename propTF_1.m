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
function [ u2 ] = propTF_1( u1, Lx, Ly , lambda, z )
%propagation - transfer function approach
%   assumes same x and y side lengths and uniform sampling
%   u1 - source plane field
%   Lx, Ly -- source and observation plane side length
%   lambda -- wavelength
%   z - propagation distance
%   u2 - observation plane field

[Ny, Nx] = size(u1);          %get input field array size
dx = Lx / Nx;
dy = Ly / Ny;%sample interval
k = 2 * pi / lambda;           %wavenumber

fx = -1 / (2*dx) : 1 / Lx : 1 / (2*dx) - 1 / Lx;
fy = -1 / (2*dy) : 1 / Ly : 1 / (2*dy) - 1 / Ly;%freq coords
[FX, FY] = meshgrid(fx, fy);

H = exp(-1i * pi * lambda * z * (FX.^2 + FY.^2));      %trans func
H = fftshift(H);                             %shift trans func
U1 = fft2(fftshift(u1));                     %shift, fft src field
U2 = H .* U1;                                  %multiply
u2 = ifftshift(ifft2(U2));                   %inv fft, center obs field

end

