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
lattConst = [3.8396, 6.6504, 0]; % [a b]
sliceDist = [3.1350, 3.1350, 3.1350]; % distance between each slice
expanNum = [12, 7]; % expand the unit cell by Expan_Nx = 3M and Expan_Ny = 2M, adaptive integer M
% slices: Each column stands for an atom
sliceA = [14, 14, 14, 14; 0.75, 0.75, 0.25, 0.25; 0.0833, 0.4167, 0.5833, 0.9167];
sliceB = [14, 14, 14, 14; 0.75, 0.25, 0.25, 0.75; 0.0833, 0.25, 0.5833, 0.75];
sliceC = [14, 14, 14, 14; 0.25, 0.75, 0.75, 0.25; 0.25, 0.4167, 0.75, 0.9167];
%% basic settings
% sampling:
Lx = expanNum(1) * lattConst(1);
Ly = expanNum(2) * lattConst(2);
Nx = 512;
Ny = 512;
fx = InitFreqAxis(Lx, Nx);
fy = InitFreqAxis(Ly, Ny);
% STEM settings:
params.KeV = 100;
interCoeff = InteractionCoefficient(params.KeV);
wavLen = HighEnergyWavLen_X(params.KeV);
params.amax = 8;
params.Cs = 0;
params.df = 0;

%% Transmission functions
% slice A projected potential:
projPotA = MultiProjPot_conv_0(sliceA, expanNum, lattConst, Lx, Ly, Nx, Ny);
% slice B projected potential:
projPotB = MultiProjPot_conv_0(sliceB, expanNum, lattConst, Lx, Ly, Nx, Ny);
% slice B projected potential:
projPotC = MultiProjPot_conv_0(sliceC, expanNum, lattConst, Lx, Ly, Nx, Ny);

% % test
% figure;
% imagesc(x, y, projPotA);
% colormap('gray');
% figure;
% imagesc(x, y, projPotB);
% colormap('gray');
% figure;
% imagesc(x, y, projPotC);
% colormap('gray');

tfA = exp(1i * interCoeff * projPotA / 1000);
tfB = exp(1i * interCoeff * projPotB / 1000);
tfC = exp(1i * interCoeff * projPotC / 1000);
tfA = BandwidthLimit(tfA, Lx, Ly, Nx, Ny, 0.67);
tfB = BandwidthLimit(tfB, Lx, Ly, Nx, Ny, 0.67);
tfC = BandwidthLimit(tfC, Lx, Ly, Nx, Ny, 0.67);
transFuncs(:, :, 1) = tfA;
transFuncs(:, :, 2) = tfB;
transFuncs(:, :, 3) = tfC;
%% CBED module
stackNum = 21;
probe = ProbeCreate(params, 0, 0, Lx, Ly, Nx, Ny);
wave = multislice(probe, wavLen, Lx, Ly, transFuncs, sliceDist, stackNum);
wave = ifftshift(fft2(fftshift(wave)));

%% show CBED pattern
logCoeff = 300;
waveI = abs(wave.^2);
logWaveI = log(1 + logCoeff * waveI / max(waveI, [], 'all'));

figure;
imagesc(fx, fy, logWaveI);
colormap('gray');
axis square;
xlabel('$ f_x (\AA^{-1}) $', 'interpreter', 'latex');
ylabel('$ f_y (\AA^{-1}) $', 'interpreter', 'latex');