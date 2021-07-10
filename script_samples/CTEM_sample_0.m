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
% BF-CTEM sample: GaAs [110]
clc;
close all;
clear all;
%% Lattice generation: GaAs [110]
lattConst = [3.995, 5.65, 0]; % [a b]
sliceDist = [1.998, 1.998]; % distance between each slice
expanNum = [6, 4]; % expand the unit cell by Expan_Nx = 3 and Expan_Ny = 2, adaptive
% slices: Each column stands for an atom
sliceA = [31, 33; 0, 0.5; 0, 0.75];
sliceB = [33, 31; 0, 0.5; 0.25, 0.5];
%% basic settings
% sampling:
Lx = expanNum(1) * lattConst(1);
Ly = expanNum(2) * lattConst(2);
Nx = 512;
Ny = 512;
dx = Lx / Nx;
dy = Ly / Ny;
x = InitAxis(Lx, Nx);
y = InitAxis(Ly, Ny);

% STEM settings:
params.KeV = 200;
interCoeff = InteractionCoefficient(params.KeV);
wavLen = HighEnergyWavLen_X(params.KeV);
params.amax = 10.67;
params.Cs = 1.3;
params.df = 700;

%% Transmission functions
% slice A projected potential:
projPotA = MultiProjPot_conv_0(sliceA, expanNum, lattConst, Lx, Ly, Nx, Ny);
% slice B projected potential:
projPotB = MultiProjPot_conv_0(sliceB, expanNum, lattConst, Lx, Ly, Nx, Ny);
% test
figure;
imagesc(x, y, projPotA);
colormap('gray');
figure;
imagesc(x, y, projPotB);
colormap('gray');

tfA = exp(1i * interCoeff * projPotA / 1000);
tfB = exp(1i * interCoeff * projPotB / 1000);
tfA = BandwidthLimit(tfA, Lx, Ly, Nx, Ny, 0.67);
tfB = BandwidthLimit(tfB, Lx, Ly, Nx, Ny, 0.67);
transFuncs(:, :, 1) = tfA;
transFuncs(:, :, 2) = tfB;

%% BF-CTEM
wave = ones(Ny, Nx);
stackNum = 50;
wave = multislice(wave, wavLen, Lx, Ly, transFuncs, sliceDist, stackNum);
wave = fft2(fftshift(wave));
otf = fftshift(AberrationFunction(params, Lx, Ly, Nx, Ny));
wave = wave .* otf;
wave = ifftshift(ifft2(wave));
waveI = abs(wave.^2);
figure;
imagesc(x, y, waveI);
colormap('gray');
axis square;
xlabel('$ x (\AA) $', 'interpreter', 'latex');
ylabel('$ y (\AA) $', 'interpreter', 'latex');