function [probe] = GenerateProbe_X(OTF, xp, yp, Lx, Ly, Nx, Ny)
%GenerateProbe_X.m generates an electron probe.
%   OTF -- prepared objective transfer function in reciprocal space;
%   Lx, Ly, Nx, Ny -- sampling parameters, L denotes side length and N the
%       sampling number in real space;
%   xp, yp -- probe position in real space;
% Note: X denotes an experimental version!

dx = Lx / Nx;
dy = Ly / Ny;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[FX, FY] = meshgrid(fx, fy);
probe = ifftshift(ifft2(fftshift(OTF .* exp(-1i * 2 * pi * (FX * xp + FY * yp)))));

NormCoeff = sqrt(sum(sum(abs(probe.^2))) * dx * dy);
probe = probe / NormCoeff;

end

