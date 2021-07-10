function [ u2, Lx_2, Ly_2 ] = propFF_1( u1, Lx_1, Ly_1, waveLength, propDist )
% Fraunhofer diffraction.
% Input:
%   u1 -- incident wave function;
%   Lx_1, Ly_1 -- incident wave side length;
%   waveLength -- wavelength;
%   propDist -- propagation distance;
% Output:
%   u2 -- exit wave function;
%   Lx_2, Ly_2 -- exit wave side length;

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

