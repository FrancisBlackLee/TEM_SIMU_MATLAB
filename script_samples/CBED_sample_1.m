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
% CBED_sample_1.m: silicon [111]
% Image CBED pattern
clc;
close all;
clear all;
%% Lattice generation: silicon [111]
LattConst = [3.8396, 6.6504, 0]; % [a b]
LayerDist = [3.1350, 3.1350, 3.1350]; % distance between each slice
CellNum = [10, 6]; % expand the unit cell by Expan_Nx = 3M and Expan_Ny = 2M, adaptive integer M
DistError = 1e-2;
% Laters: Each column for an atom
LayerA = [14, 14, 14, 14; 0.75, 0.75, 0.25, 0.25; 0.0833, 0.4167, 0.5833, 0.9167];
LayerB = [14, 14, 14, 14; 0.75, 0.25, 0.25, 0.75; 0.0833, 0.25, 0.5833, 0.75];
LayerC = [14, 14, 14, 14; 0.25, 0.75, 0.75, 0.25; 0.25, 0.4167, 0.75, 0.9167];
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
Params.KeV = 100;
InterCoeff = InteractionCoefficient(Params.KeV);
WaveLength = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
WaveNumber = 2 * pi / WaveLength;     %wavenumber
Params.amax = 8;
Params.Cs = 0;
Params.df = 0;
%% Transmission functions
% Layer A:
Proj_PotA = MultiProjPot_conv_0(LayerA, CellNum, LattConst, Lx, Ly, Nx, Ny);
% Layer B:
Proj_PotB = MultiProjPot_conv_0(LayerB, CellNum, LattConst, Lx, Ly, Nx, Ny);
% Layer B:
Proj_PotC = MultiProjPot_conv_0(LayerC, CellNum, LattConst, Lx, Ly, Nx, Ny);

% % test
% figure;
% imagesc(x, y, Proj_PotA);
% colormap('gray');
% figure;
% imagesc(x, y, Proj_PotB);
% colormap('gray');
% figure;
% imagesc(x, y, Proj_PotC);
% colormap('gray');

TF_A = exp(1i * InterCoeff * Proj_PotA / 1000);
TF_B = exp(1i * InterCoeff * Proj_PotB / 1000);
TF_C = exp(1i * InterCoeff * Proj_PotC / 1000);
TF_A = BandwidthLimit(TF_A, Lx, Ly, Nx, Ny, 0.67);
TF_B = BandwidthLimit(TF_B, Lx, Ly, Nx, Ny, 0.67);
TF_C = BandwidthLimit(TF_C, Lx, Ly, Nx, Ny, 0.67);
TransFuncs(:, :, 1) = TF_A;
TransFuncs(:, :, 2) = TF_B;
TransFuncs(:, :, 3) = TF_C;
%% CBED module
StackNum = 50;
Probe = ProbeCreate(Params, 0, 0, Lx, Ly, Nx, Ny);
TransWave = multislice(Probe, WaveLength, Lx, Ly, TransFuncs, LayerDist, StackNum);
% Trans_Wave_Far = ifftshift(fft2(fftshift(TransWave)) * dx * dy);
Trans_Wave_Far = ifftshift(fft2(fftshift(TransWave)));
DetectInten = log(1 + 0.1 * abs(Trans_Wave_Far.^2));
% DetectInten = abs(Trans_Wave_Far.^2);

% Show the detected image:
figure;
subplot(1, 2, 1);
imagesc(fx, fy, DetectInten);
colormap('gray');
axis square;
subplot(1, 2, 2);
plot(fx, DetectInten(Ny / 2 + 1, : ));

thickness = sum(LayerDist(:)) * StackNum;