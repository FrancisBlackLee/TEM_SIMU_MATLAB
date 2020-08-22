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
% ADF-STEM sample 3: C:mp-568286
clc;
close all;
clear all;
%% Lattice generation: silicon [111]
unitCell = load('E:\Group_Affairs\Huang_Guangyi\C Coordinates.txt');
atomNumInUnitCell = size(unitCell, 1);
unitCell = [unitCell(:, 1), ones(atomNumInUnitCell, 1), unitCell(:, 3),...
    unitCell(:, 4), unitCell(:, 5)]';

lattConstA = 4.2745;
lattConstB = 4.9347;
lattConstC = 8.0283;

[slice, sliceDist, extraSlice] = CrystalSlicing_X(unitCell, unitCell, 0.2, 1, 0, 0);

lattConst = [lattConstA, lattConstB];
sliceDist = sliceDist * lattConstC;

expanNum = [4, 4];

%% sampling:
Lx = expanNum(1) * lattConst(1);
Ly = expanNum(2) * lattConst(2);
Nx = 428;
Ny = 490;
dx = Lx / Nx;
dy = Ly / Ny;
%% STEM settings:
Params.KeV = 300;
InterCoeff = InteractionCoefficient(Params.KeV);
WavLen = HighEnergyWavLen_X(Params.KeV);
Params.aperture = CircApert_X(Lx, Ly, Nx, Ny, WavLen, 23);
Params.Cs3 = 0;
Params.Cs5 = 0;
Params.df = 0;
Params.scanx = linspace(0, lattConstA, 46);
Params.scany = linspace(0, lattConstB, 52);

HighAngle = 200; % im mrad
LowAngle = 40; %in mrad

Params.detector = AnularDetector_X(LowAngle, HighAngle, WavLen, Lx, Ly, Nx, Ny);
%% Transmission functions:
% Layer A:
ProjPotA = MultiProjPot_conv_X(slice{1}, expanNum, lattConst, Lx, Ly, Nx, Ny, 1.0e-5);
% Layer B:
ProjPotB = MultiProjPot_conv_X(slice{2}, expanNum, lattConst, Lx, Ly, Nx, Ny, 1.0e-5);

% test
figure;
imagesc(ProjPotA);
colormap('gray'); axis equal;
figure;
imagesc(ProjPotB);
colormap('gray'); axis equal;

tic;
TF_A = exp(1i * InterCoeff * ProjPotA / 1000);
TF_B = exp(1i * InterCoeff * ProjPotB / 1000);

TF_A = BandwidthLimit(TF_A, Lx, Ly, Nx, Ny, 0.67);
TF_B = BandwidthLimit(TF_B, Lx, Ly, Nx, Ny, 0.67);

TransFuncs(:, :, 1) = TF_A;
TransFuncs(:, :, 2) = TF_B;

%% Imaging section:
stemImg = ADF_STEM_X(Lx, Ly, Params, TransFuncs, sliceDist, 10, 0);
toc;

% repStemImg = repmat(stemImg, 6, 10);

figure;
imagesc(stemImg);
colormap('gray');