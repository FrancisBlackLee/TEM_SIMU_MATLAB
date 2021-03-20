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
% CBED_sample_3.m: GaAs [110]
% Image CBED pattern
clc;
clear;
close all;
%% Lattice generation: GaAs [110]
LattConst = [3.995, 5.65, 0]; % [a b]
sliceDist = [1.998, 1.998]; % distance between each slice
cellNum = [10, 7]; % expand the unit cell
DistError = 1e-2;
% Laters: Each column for an atom
sliceA = [31, 33; 0, 0.5; 0, 0.75];
sliceB = [33, 31; 0, 0.5; 0.25, 0.5];
%% basic settings
% sampling:
Lx = cellNum(1) * LattConst(1);
Ly = cellNum(2) * LattConst(2);
Nx = 512;
Ny = 512;

% STEM settings:
params.KeV = 100;
interCoeff = InteractionCoefficient(params.KeV);
wavLen = HighEnergyWavLen_X(params.KeV);
params.amax = 5;
params.Cs = 0;
params.df = 0;
%% Transmission functions
% Layer A:
Proj_PotA = MultiProjPot_conv_0(sliceA, cellNum, LattConst, Lx, Ly, Nx, Ny);
% Layer B:
Proj_PotB = MultiProjPot_conv_0(sliceB, cellNum, LattConst, Lx, Ly, Nx, Ny);

TF_A = exp(1i * interCoeff * Proj_PotA / 1000);
TF_B = exp(1i * interCoeff * Proj_PotB / 1000);
TF_A = BandwidthLimit(TF_A, Lx, Ly, Nx, Ny, 0.67);
TF_B = BandwidthLimit(TF_B, Lx, Ly, Nx, Ny, 0.67);
TransFuncs(:, :, 1) = TF_A;
TransFuncs(:, :, 2) = TF_B;
%% Scanning module
Probe = ProbeCreate(params, 0, 0, Lx, Ly, Nx, Ny);
TransWave = multislice(Probe, wavLen, Lx, Ly, TransFuncs, sliceDist, 100);
Trans_Wave_Far = ifftshift(fft2(fftshift(TransWave)));
DetectInten = log(1 + 0.1 * abs(Trans_Wave_Far.^2));

% Show the CBED pattern:
figure;
imagesc(DetectInten);
colormap('gray');
axis square;
axis off;