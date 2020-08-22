%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2020  Francis Black Lee and Li Xian

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BF-CTEM sample: silicon [110]
clc;
close all;
clear all;
%% Lattice generation: silicon [110]
LattConst = [3.8396, 5.4300, 0]; % [a b]
LayerDist = [1.9198, 1.9198]; % distance between each slice
CellNum = 10 * [3, 2]; % expand the unit cell by Expan_Nx = 3 and Expan_Ny = 2, adaptive
DistError = 1e-2;
% Laters: Each column for an atom
LayerA = [14, 14; 0, 0.5; 0, 0.75];
LayerB = [14, 14; 0, 0.5; 0.25, 0.5];
%% basic settings
% sampling:
Lx = CellNum(1) * LattConst(1);
Ly = CellNum(2) * LattConst(2);
Nx = 1024;
Ny = 1024;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[Fx, Fy] = meshgrid(fx, fy);
% STEM settings:
Params.KeV = 200;
InterCoeff = InteractionCoefficient(Params.KeV);
WaveLength = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
WaveNumber = 2 * pi / WaveLength;     %wavenumber

%% Transmission functions
% Layer A:
Proj_PotA = MultiProjPot_conv_0(LayerA, CellNum, LattConst, Lx, Ly, Nx, Ny);
% Layer B:
Proj_PotB = MultiProjPot_conv_0(LayerB, CellNum, LattConst, Lx, Ly, Nx, Ny);
% % test
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

%% Multislice
IncidentWave = ones(Ny, Nx);
TransWave = multislice(IncidentWave, WaveLength, Lx, Ly, TransFuncs, LayerDist, 100);

%% Imaging
Params.amax = 20;
Params.Cs = 0;

Params.df = 1180;

ReciTransWave = fft2(fftshift(TransWave));
ObjLens = fftshift(AberrationFunction(Params, Lx, Ly, Nx, Ny));
ReciTFWave = ReciTransWave.*ObjLens;
TFWave=ifftshift(ifft2(ReciTFWave));
Image = abs(TFWave.^2);
figure;
imagesc(x, y, Image);
colormap('gray');
axis square;
title(['df = ', num2str(Params.df), 'Angs.']);

% s = 11;
% for i = 1 : 20
%     Params.df = i * 5 + s * 100;
% 
%     ReciTransWave = fft2(fftshift(TransWave));
%     ObjLens = fftshift(AberrationFunction(Params, Lx, Ly, Nx, Ny));
%     ReciTFWave = ReciTransWave.*ObjLens;
%     TFWave=ifftshift(ifft2(ReciTFWave));
%     Image = abs(TFWave.^2);
%     figure;
%     imagesc(x, y, Image);
%     colormap('gray');
%     axis square;
%     title(['df = ', num2str(Params.df), 'Angs.']);
% end