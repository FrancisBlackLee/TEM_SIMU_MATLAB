% DPC_nvstgt_0: silicon [110]
% Image CBED pattern
clc;
close all;
clear all;
%% Lattice generation: silicon [110]
Lattice_Const = [3.84, 5.43]; % [a b]
LayerDist = [1.9198, 1.9198]; % distance between each slice
M = 1;
Cell_Num = [3 * M, 2 * M]; % expand the unit cell by Expan_Nx = 3M and Expan_Ny = 2M, adaptive integer M
DistError = 1e-2;
% Laters: Each column for an atom
LayerA = [0, 0.5; 0, 0.75];
LayerB = [0, 0.5; 0.25, 0.5];
% Expansion
LayerA = SquareLattExpan(LayerA, Lattice_Const, Cell_Num);
LayerB = SquareLattExpan(LayerB, Lattice_Const, Cell_Num);
%% basic settings
% sampling:
Lx = Cell_Num(1) * Lattice_Const(1);
Ly = Cell_Num(2) * Lattice_Const(2);
Nx = 256;
Ny = 256;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[Fx, Fy] = meshgrid(fx, fy);
% STEM settings:
Params.KeV = 100;
InterCoeff = InteractionCoefficient(Params.KeV);
WaveLength = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
WaveNumber = 2 * pi / WaveLength;     %wavenumber
Params.amax = 11.43;
Params.Cs = 1.3;
Params.df = 850;
%% Transmission functions
% Layer A:
Proj_PotA = ProjectedPotential(Lx, Ly, Nx, Ny, 14 * ones(size(LayerA, 2), 1), LayerA(1, :), LayerA(2, :));
% Layer B:
Proj_PotB = ProjectedPotential(Lx, Ly, Nx, Ny, 14 * ones(size(LayerB, 2), 1), LayerB(1, :), LayerB(2, :));
% test
% figure;
% imagesc(x, y, Proj_PotA);
% colormap('gray');
% figure;
% imagesc(x, y, Proj_PotB);
% colormap('gray');

TF_A = exp(1i * InterCoeff * Proj_PotA / 1000);
TF_B = exp(1i * InterCoeff * Proj_PotB / 1000);
TF_A = BandwidthLimit(TF_A, Lx, Ly, Nx, Ny, 0.67);
TF_B = BandwidthLimit(TF_B, Lx, Ly, Nx, Ny, 0.67);
TransFuncs(:, :, 1) = TF_A;
TransFuncs(:, :, 2) = TF_B;
%% Scanning module
SaveDir = 'D:\Francis. B. Lee\cooperation\Group Cooperation\ZJR\WaveFuncs';
SaveSliceSeries = [1, 3, 7, 8, 9, 15, 19, 20, 21];

Probe = ProbeCreate(Params, 0, 0, Lx, Ly, Nx, Ny);
TransWave = multislice(Probe, WaveLength, Lx, Ly, TransFuncs, LayerDist, 10, SaveDir, SaveSliceSeries);
Trans_Wave_Far = ifftshift(fft2(fftshift(TransWave)) * dx * dy);
DetectInten = log(abs(Trans_Wave_Far.^2));
% DetectInten = abs(Trans_Wave_Far.^2);

% Show the detected image:
figure;
imagesc(x, y, DetectInten);
colormap('gray');
axis square;