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
% ADF-STEM sample 2: Si <1 1 1>
clc;
close all;
clear all;
%% Lattice generation: silicon [111]
LattConst = [3.8396, 6.6504, 0]; % [a b]
SliceDist = [3.1350, 3.1350, 3.1350]; % distance between each slice
CellNum = [5, 3]; % expand the unit cell by Expan_Nx = 3M and Expan_Ny = 2M, adaptive integer M
DistError = 1e-2;
% Laters: Each column for an atom
LayerA = [14, 14, 14, 14; 0.75, 0.75, 0.25, 0.25; 0.0833, 0.4167, 0.5833, 0.9167];
LayerB = [14, 14, 14, 14; 0.75, 0.25, 0.25, 0.75; 0.0833, 0.25, 0.5833, 0.75];
LayerC = [14, 14, 14, 14; 0.25, 0.75, 0.75, 0.25; 0.25, 0.4167, 0.75, 0.9167];

%% sampling:
Lx = CellNum(1) * LattConst(1);
Ly = CellNum(2) * LattConst(2);
Nx = 512;
Ny = 512;
dx = Lx / Nx;
dy = Ly / Ny;
%% STEM settings:
Params.KeV = 300;
InterCoeff = InteractionCoefficient(Params.KeV);
WavLen = HighEnergyWavLen_X(Params.KeV);
Params.aperture = CircApert_X(Lx, Ly, Nx, Ny, WavLen, 18);
Params.Cs3 = 0;
Params.Cs5 = 0;
Params.df = 0;
Params.scanx = linspace(0, 3.8396, 39);
Params.scany = linspace(0, 6.6504, 67);

HighAngle = 200; % im mrad
LowAngle = 40; %in mrad

Params.detector = AnnularDetector_X(LowAngle, HighAngle, WavLen, Lx, Ly, Nx, Ny);
%% Transmission functions:
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

tic;
TF_A = exp(1i * InterCoeff * Proj_PotA / 1000);
TF_B = exp(1i * InterCoeff * Proj_PotB / 1000);
TF_C = exp(1i * InterCoeff * Proj_PotC / 1000);
TF_A = BandwidthLimit(TF_A, Lx, Ly, Nx, Ny, 0.67);
TF_B = BandwidthLimit(TF_B, Lx, Ly, Nx, Ny, 0.67);
TF_C = BandwidthLimit(TF_C, Lx, Ly, Nx, Ny, 0.67);
TransFuncs(:, :, 1) = TF_A;
TransFuncs(:, :, 2) = TF_B;
TransFuncs(:, :, 3) = TF_C;
%% Imaging section:
stemImg = ADF_STEM_X(Lx, Ly, Params, TransFuncs, SliceDist, 21, 0);
toc;

% repStemImg = repmat(stemImg, 6, 10);

figure;
imagesc(stemImg);
colormap('gray'); axis square;