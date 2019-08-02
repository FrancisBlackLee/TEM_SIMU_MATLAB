function [limited] = BandwidthLimit(unlimited, Lx, Ly, Nx, Ny, proportion)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dx = Lx / Nx;
dy = Ly / Ny;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[Fx, Fy] = meshgrid(fx, fy);
FreqSquare = Fx.^2 + Fy.^2;
Aperture = ones(size(Fx));
R_apert = proportion * min(1 / (2 * dx), 1 / (2 * dy));
Aperture(FreqSquare > R_apert^2) = 0;
reciprocal = fftshift(Aperture) .* fft2(fftshift(unlimited));
limited = ifftshift(ifft2(reciprocal));

end

