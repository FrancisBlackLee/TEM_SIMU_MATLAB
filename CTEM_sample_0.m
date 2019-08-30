% BF-CTEM sample: silicon [110]
clc;
close all;
clear all;
%% Lattice generation: silicon [110]
LattConst = [3.8396, 5.4300, 0]; % [a b]
LayerDist = [1.9198, 1.9198]; % distance between each slice
CellNum = [3, 2]; % expand the unit cell by Expan_Nx = 3 and Expan_Ny = 2, adaptive
DistError = 1e-2;
% Laters: Each column for an atom
LayerA = [14, 14; 0, 0.5; 0, 0.75];
LayerB = [14, 14; 0, 0.5; 0.25, 0.5];
%% basic settings
% sampling:
Lx = CellNum(1) * LattConst(1);
Ly = CellNum(2) * LattConst(2);
Nx = 512;
Ny = 512;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[Fx, Fy] = meshgrid(fx, fy);
% STEM settings:
Params.KeV = 400;
InterCoeff = InteractionCoefficient(Params.KeV);
WaveLength = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
WaveNumber = 2 * pi / WaveLength;     %wavenumber
Params.amax = 12;
Params.Cs = 1.3;
Params.df = 566;

%% Transmission functions
% Layer A:
Proj_PotA = MultiProjPot_conv_0(LayerA, CellNum, LattConst, Lx, Ly, Nx, Ny);
% Layer B:
Proj_PotB = MultiProjPot_conv_0(LayerB, CellNum, LattConst, Lx, Ly, Nx, Ny);
% test
figure;
imagesc(x, y, Proj_PotA);
colormap('gray');
figure;
imagesc(x, y, Proj_PotB);
colormap('gray');

TF_A = exp(1i * InterCoeff * Proj_PotA / 1000);
TF_B = exp(1i * InterCoeff * Proj_PotB / 1000);
TF_A = BandwidthLimit(TF_A, Lx, Ly, Nx, Ny, 0.67);
TF_B = BandwidthLimit(TF_B, Lx, Ly, Nx, Ny, 0.67);
TransFuncs(:, :, 1) = TF_A;
TransFuncs(:, :, 2) = TF_B;

%% BF-CTEM
IncidentWave = ones(Ny, Nx);
TransWave = multislice(IncidentWave, WaveLength, Lx, Ly, TransFuncs, LayerDist, 20);
ReciTransWave = fft2(fftshift(TransWave));
ObjLens = fftshift(AberrationFunction(Params, Lx, Ly, Nx, Ny));
ReciTFWave = ReciTransWave.*ObjLens;
TFWave=ifftshift(ifft2(ReciTFWave));
Image = abs(TFWave.^2);
figure;
imagesc(x, y, Image);
colormap('gray');