function [Probe] = BesselBeam(Params, xp, yp, Lx, Ly, Nx, Ny)
%PROBECREATE.M calculates the probe function with the input STEM
%parameters.
%   xp, yp --position of the center of the probe
%   Lx, Ly --side lengths
%   Nx, Ny --number of sampling

dx = Lx / Nx;
dy = Ly / Ny;
aberr_TF = BesselAberration(Params, Lx, Ly, Nx, Ny);
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[Fx, Fy] = meshgrid(fx, fy);
Probe = ifft2(fftshift(aberr_TF .* exp(-1i * 2 * pi * (Fx * xp + Fy * yp)))) / (dx * dy);
Probe = ifftshift(Probe);
NormEffi = sqrt(sum(sum(abs(Probe.^2))) * dx * dy);
Probe = Probe / NormEffi;

end

