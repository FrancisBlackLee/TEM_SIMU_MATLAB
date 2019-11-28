function [Aperture] = CircApert_X(Lx, Ly, Nx, Ny, WavLen, NA)
%CircApert_X.m generates a circular aperture in reciprocal space.
%   Lx, Ly, Nx, Ny -- sampling parameters, L denotes side length and N the
%       sampling number in real space;
%   WavLen -- wavelength;
%   NA -- numerical aperture in mrad;
% Note: X denotes an experimental version!

dx = Lx / Nx;
dy = Ly / Ny;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[FX, FY] = meshgrid(fx, fy);
FreqSqu = FX.^2 + FY.^2;

Aperture = (FreqSqu < (NA * 1e-3 / WavLen)^2);

end

