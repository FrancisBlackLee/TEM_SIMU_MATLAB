function [PropKernel] = FresnelPropKernel_X(Lx, Ly, Nx, Ny, WavLen, PropDist)
%FresnelPropKernel_X.m computes the Fresnel propagation kernel.
%   Lx, Ly -- sampling side length;
%   Nx, Ny -- sampling number;
%   WavLen -- wavelength;
%   PropDist -- propagation distance
% Note: X denotes an experimental version!

dx = Lx / Nx;
dy = Ly / Ny;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[FX, FY] = meshgrid(fx, fy);
PropKernel = exp(-1i * pi * WavLen * PropDist * (FX.^2 + FY.^2));

end

