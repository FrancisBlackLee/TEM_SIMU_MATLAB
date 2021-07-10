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
% ADF-STEM sample 0: GaAs [1 1 0]
clc;
close all;
clear all;
%% specimen preparation:
lattConst = [3.995, 5.65, 0]; % [a b]
sliceDist = [1.998, 1.998]; % distance between each slice
expanNum = 4 * [3, 2]; % expand the unit cell by Expan_Nx = 3 and Expan_Ny = 2, adaptive
% Laters: Each column for an atom
sliceA = [31,  31,  31,  31,  33;...
          1,   1,   1,   1,   1;...
          0,   1,   1,   0,   0.5;...
          0,   0,   1,   1,   0.75];
sliceB = [33,  33,  31;...
          1,   1,   1;...
          0,   1,   0.5;...
          0.25,0.25,0.5];

%% sampling:
Lx = expanNum(1) * lattConst(1);
Ly = expanNum(2) * lattConst(2);
Nx = 512;
Ny = 512;
dx = Lx / Nx;
dy = Ly / Ny;
%% STEM settings:
params.KeV = 200;
interCoeff = InteractionCoefficient(params.KeV);
wavLen = HighEnergyWavLen_X(params.KeV);
params.aperture = CircApert_X(Lx, Ly, Nx, Ny, wavLen, 18);
params.Cs3 = 0;
params.Cs5 = 0;
params.df = 0;
params.scanx = linspace(0, 3.995, 41);
params.scany = linspace(0, 5.65, 57);

adfHighAngle = 200; % im mrad
adfLowAngle = 40; %in mrad

params.detector = AnnularDetector_X(adfLowAngle, adfHighAngle, wavLen, Lx, Ly, Nx, Ny);
%% Transmission functions:
projPotA = MultiProjPot_conv_X(sliceA, expanNum, lattConst, Lx, Ly, Nx, Ny, 1e-5);
% figure;
% imagesc(x, y, projPotA);
% colormap('gray'); axis square;
% title('Proj.Pot. A');

projPotB = MultiProjPot_conv_X(sliceB, expanNum, lattConst, Lx, Ly, Nx, Ny, 1e-5);
% figure;
% imagesc(x, y, projPotB);
% colormap('gray'); axis square;
% title('Proj.Pot. B');
tic;
tfA = exp(1i * interCoeff * projPotA / 1000);
tfB = exp(1i * interCoeff * projPotB / 1000);
tfA = BandwidthLimit(tfA, Lx, Ly, Nx, Ny, 0.67);
tfB = BandwidthLimit(tfB, Lx, Ly, Nx, Ny, 0.67);
transFuncs(:, :, 1) = tfA;
transFuncs(:, :, 2) = tfB;
%% Imaging section:
stackNum = 5;
stemImg = ADF_STEM_X(Lx, Ly, params, transFuncs, sliceDist, stackNum, 0);
toc;

repStemImg = repmat(stemImg, 2, 3);

figure;
imagesc(repStemImg);
colormap('gray'); axis square;