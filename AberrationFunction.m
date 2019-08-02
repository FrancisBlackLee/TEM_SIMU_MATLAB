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

