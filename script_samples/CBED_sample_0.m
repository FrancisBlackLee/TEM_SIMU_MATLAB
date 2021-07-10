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
% CBED_sample_0.m: silicon [110]
% Image CBED pattern
clc;
close all;
clear all;
%% Lattice generation: silicon [110]
lattConst = [3.84, 5.43, 0]; % [a b]
sliceDist = [1.9198, 1.9198]; % distance between each slice
expanNum = [24, 17]; % expand the unit cell by Expan_Nx = 3M and Expan_Ny = 2M, adaptive integer M
% slices: Each column stands for an atom
sliceA = [14, 14; 0, 0.5; 0, 0.75];
sliceB = [14, 14; 0, 0.5; 0.25, 0.5];
%% basic settings
% sampling:
Lx = expanNum(1) * lattConst(1);
Ly = expanNum(2) * lattConst(2);
Nx = 1024;
Ny = 1024;
fx = InitFreqAxis(Lx, Nx);
fy = InitFreqAxis(Ly, Ny);
% STEM settings:
params.KeV = 100;
interCoeff = InteractionCoefficient(params.KeV);
wavLen = HighEnergyWavLen_X(params.KeV);
params.amax = 6;
params.Cs = 0;
params.df = 0;
%% Transmission functions
% slice A projected potential:
projPotA = MultiProjPot_conv_0(sliceA, expanNum, lattConst, Lx, Ly, Nx, Ny);
% slice B projected potential:
projPotB = MultiProjPot_conv_0(sliceB, expanNum, lattConst, Lx, Ly, Nx, Ny);

% % test
% figure;
% imagesc(x, y, projPotA);
% colormap('gray');
% figure;
% imagesc(x, y, projPotB);
% colormap('gray');

tfA = exp(1i * interCoeff * projPotA / 1000);
tfB = exp(1i * interCoeff * projPotB / 1000);
tfA = BandwidthLimit(tfA, Lx, Ly, Nx, Ny, 0.67);
tfB = BandwidthLimit(tfB, Lx, Ly, Nx, Ny, 0.67);
transFuncs(:, :, 1) = tfA;
transFuncs(:, :, 2) = tfB;
%% diffraction
probe = ProbeCreate(params, 0, 0, Lx, Ly, Nx, Ny);
stackNum = 150;
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