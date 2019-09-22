% DPC_nvstgt_5.m studys the probe:
clc;
clear all;
close all;
%% Sampling settings:
Lx = 15;
Ly = Lx;
Nx = 512;
Ny = 512;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
fx = -1 / (2*dx) : 1 / Lx : 1 / (2*dx) - 1 / Lx;
fy = -1 / (2*dy) : 1 / Ly : 1 / (2*dy) - 1 / Ly;
[FX, FY] = meshgrid(fx, fy);
Fsqu = FX.^2 + FY.^2;
%% STEM settings:
Params.KeV = 300;
InterCoeff = InteractionCoefficient(Params.KeV);
WaveLength = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
WaveNumber = 2 * pi / WaveLength;     %wavenumber
Params.amax = 23;
Params.Cs = 0;
Params.df = 0;
%% Probe:
% Generate a probe centering at the origin:
Probe = ProbeCreate(Params, 0, 0, Lx, Ly, Nx, Ny);
FT_Probe = ifftshift(fft2(fftshift(Probe))) * dx * dy;
abs_FT_Probe = abs(FT_Probe);
% Show FT_Probe:
figure;
subplot(4, 2, 1);
imagesc(fx, fy, real(FT_Probe));
colormap('gray'); axis square;
title('real(FT{Probe})');
subplot(4, 2, 2);
plot(fx, real(FT_Probe(Ny / 2 + 1, : )));
title('real');
subplot(4, 2, 3);
imagesc(fx, fy, imag(FT_Probe));
colormap('gray'); axis square;
title('imag(FT{Probe})');
subplot(4, 2, 4);
plot(fx, imag(FT_Probe(Ny / 2 + 1, : )));
title('imag');
subplot(4, 2, 5);
imagesc(fx, fy, abs_FT_Probe);
colormap('gray'); axis square;
title('abs(FT{Probe})');
subplot(4, 2, 6);
plot(fx, abs_FT_Probe(Ny / 2 + 1, : ));
title('abs');
rszCoeff = 0.5;
rsz_abs_FT_Probe = imresize(abs_FT_Probe, rszCoeff);
rsz_FX = imresize(FX, rszCoeff);
rsz_FY = imresize(FY, rszCoeff);
[rszNx, rszNy] = size(rsz_FX);
% Show the resized abs_FT_Probe:
subplot(4, 2, 7);
imagesc(rsz_FX(1, : ), rsz_FY( : , 1), rsz_abs_FT_Probe);
colormap('gray'); axis square;
title('rsz_abs_FT_Probe');
subplot(4, 2, 8);
plot(rsz_FX(1, : ), rsz_abs_FT_Probe(rszNy / 2 + 1, : ));
title('rsz');
%% ProbeInten:
ProbeInten = abs(Probe.^2);
% Show ProbeInten:
figure;
subplot(1, 2, 1);
imagesc(x, y, ProbeInten);
colormap('gray'); axis square;
title('ProbeInten');
subplot(1, 2, 2);
plot(x, ProbeInten(Ny / 2 + 1, : ));

% Fourier transform the probe intensity:
FT_ProbeInten = ifftshift(fft2(fftshift(ProbeInten))) * dx * dy;
% Show FT_ProbeInten:
figure;
subplot(3, 2, 1);
imagesc(fx, fy, real(FT_ProbeInten));
colormap('gray'); axis square;
title('real(FT{ProbeInten})');
subplot(3, 2, 2);
plot(fx, real(FT_ProbeInten(Ny / 2 + 1, : )));
title('real');
subplot(3, 2, 3);
imagesc(fx, fy, imag(FT_ProbeInten));
colormap('gray'); axis square;
title('imag(FT{ProbeInten})');
subplot(3, 2, 4);
plot(fx, imag(FT_ProbeInten(Ny / 2 + 1, : )));
title('imag');
% Show the modulus of FT_ProbeInten:
subplot(3, 2, 5);
imagesc(fx, fy, abs(FT_ProbeInten));
colormap('gray'); axis square;
title('abs(FT{ProbeInten})');
subplot(3, 2, 6);
plot(fx, abs(FT_ProbeInten(Ny / 2 + 1, : )));
title('abs');