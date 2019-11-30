%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019  Francis Black Lee and Li Xian

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
% ADF-STEM sample 1: GaAs <1 1 0>
clc;
close all;
clear all;
%% specimen preparation:
LattConst = [3.995, 5.65, 0]; % [a b]
SliceDist = [1.998, 1.998]; % distance between each slice
CellNum = [6, 4]; % expand the unit cell by Expan_Nx = 3 and Expan_Ny = 2, adaptive
% Laters: Each column for an atom
SliceA = [31,  31,  31,  31,  33;...
          1,   1,   1,   1,   1;...
          0,   1,   1,   0,   0.5;...
          0,   0,   1,   1,   0.75];
SliceB = [33,  33,  31;...
          1,   1,   1;...
          0,   1,   0.5;...
          0.25,0.25,0.5];

%% sampling:
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
[FX, FY] = meshgrid(fx, fy);
FreqSqu = FX.^2 + FY.^2;
%% STEM settings:
Params.KeV = 200;
Params.NA = 23;
Params.Cs3 = 0;
Params.Cs5 = 0;
Params.df = 0;
Params.scanx = linspace(-5, 0, 32);
Params.scany = linspace(3, 8, 32);

HighAngle = 200; % im mrad
LowAngle = 40; %in mrad

InterCoeff = InteractionCoefficient(Params.KeV);
WavLen = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength

Params.detector = ((FreqSqu < (HighAngle * 1e-3 / WavLen)^2) & (FreqSqu > (LowAngle * 1e-3 / WavLen)^2));
%% Transmission functions:
ProjPotA = MultiProjPot_conv_X(SliceA, CellNum, LattConst, Lx, Ly, Nx, Ny, 1e-5);
% figure;
% imagesc(x, y, ProjPotA);
% colormap('gray'); axis square;
% title('Proj.Pot. A');

ProjPotB = MultiProjPot_conv_X(SliceB, CellNum, LattConst, Lx, Ly, Nx, Ny, 1e-5);
% figure;
% imagesc(x, y, ProjPotB);
% colormap('gray'); axis square;
% title('Proj.Pot. B');

TF_A = exp(1i * InterCoeff * ProjPotA / 1000);
TF_B = exp(1i * InterCoeff * ProjPotB / 1000);
TF_A = BandwidthLimit(TF_A, Lx, Ly, Nx, Ny, 0.67);
TF_B = BandwidthLimit(TF_B, Lx, Ly, Nx, Ny, 0.67);
TransFuncs(:, :, 1) = TF_A;
TransFuncs(:, :, 2) = TF_B;
%% Imaging section:
stemImg = ADF_STEM_X(Lx, Ly, Params, TransFuncs, SliceDist, 20, 0);

figure;
imagesc(x, y, stemImg);
colormap('gray'); axis square;