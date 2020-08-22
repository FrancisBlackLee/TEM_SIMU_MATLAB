% in_situ_gold_4nm.m
clc;
close all;
clear all;
%% 4nm silica in water
ProjPotDir = 'Tests\ProjPot';
Lx = 60.0;
Ly = Lx;
Nx = 1024;
Ny = 1024;
KeV = 300;
InciWave = ones(Ny, Nx);
SliceDist = load('Tests\SliceDist_4nmGold_in_water.txt');
ExitWave = multislice_X(InciWave, KeV, Lx, Ly, 'files', SliceDist, 1, ProjPotDir, '*.txt');

dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[Fx, Fy] = meshgrid(fx, fy);
%% CTEM settings:
Params.KeV = 400;
InterCoeff = InteractionCoefficient(Params.KeV);
WaveLength = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
WaveNumber = 2 * pi / WaveLength;     %wavenumber
Params.amax = 10.67;
Params.Cs = 0;
Params.df = 500;

ReciTransWave = fft2(fftshift(ExitWave));
ObjLens = fftshift(AberrationFunction(Params, Lx, Ly, Nx, Ny));
ReciTFWave = ReciTransWave.*ObjLens;
TFWave=ifftshift(ifft2(ReciTFWave));
Image = abs(TFWave.^2);
figure;
imagesc(x, y, Image);
colormap('gray'); axis square;