% PhaseErrorTest_0.m -- test MultiAberrPhaseError_X.m
clc;
close all;
clear all;
%% Parameter setting:
% unit of aberrations is angstrom
Aberr{1} = [0, 0];
Aberr{2} = [0, 1e4];
Aberr{3} = [0, 0, 0];
Aberr{4} = [0, 0, 0];
Aberr{5} = [0, 0, 0, 0];
Lx = 12;
Ly = 12;
Nx = 1024;
Ny = 1024;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;

KeV = 200;
WavLen = 12.3986 / sqrt((2 * 511.0 + KeV) * KeV);  %wavelength

PhaseError = MultiAberrPhaseError_X(Aberr, WavLen, Lx, Ly, Nx, Ny);

figure;
subplot(1, 2, 1);
imagesc(fx, fy, PhaseError);
colormap('gray'); axis square;
subplot(1, 2, 2);
plot(fx, PhaseError(Ny / 2 + 1, : ));

amax = 20e-3; % rad;
[FX, FY] = meshgrid(fx, fy);
FreqSquare = FX.^2 + FY.^2;
Aperture = ones(size(FX));
Aperture(FreqSquare >= (sin(amax) / WavLen)^2) = 0;
Aberr_TF = exp(1i * PhaseError) .* Aperture;

xp = 0; yp = 0;
Probe = ifft2(fftshift(Aberr_TF .* exp(-1i * 2 * pi * (FX * xp + FY * yp)))) / (dx * dy);
Probe = ifftshift(Probe);
NormEffi = sqrt(sum(sum(abs(Probe.^2))) * dx * dy);
Probe = Probe / NormEffi;

figure;
imagesc(x, y, abs(Probe.^2));
colormap('gray'); axis square;