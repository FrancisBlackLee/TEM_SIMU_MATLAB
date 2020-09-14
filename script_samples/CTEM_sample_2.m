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
% BF-CTEM sample: Monolayer graphene
clc;
close all;
clear all;
%% Lattice generation: Monolayer graphene
LattConst = [2.46, 4.26, 0]; % [a b]
CellNum = [15, 9]; % expand the unit cell
DistError = 1e-2;
% Laters: Each column for an atom
slice = [6, 6, 6, 6; 0, 0.5, 0, 0.5; 0, 0.75, 0.25, 0.5];
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
Params.KeV = 200;
InterCoeff = InteractionCoefficient(Params.KeV);
WaveLength = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
WaveNumber = 2 * pi / WaveLength;     %wavenumber
Params.amax = 40;
Params.Cs = 0;
Params.df = 700;

%% Transmission functions
% Layer A:
Proj_PotA = MultiProjPot_conv_0(slice, CellNum, LattConst, Lx, Ly, Nx, Ny);
% % test
% figure;
% imagesc(x, y, Proj_PotA);
% colormap('gray');

TF_A = exp(1i * InterCoeff * Proj_PotA / 1000);

%% BF-CTEM
IncidentWave = ones(Ny, Nx);
TransWave = IncidentWave .* TF_A;
ReciTransWave = fft2(fftshift(TransWave));
ObjLens = fftshift(AberrationFunction(Params, Lx, Ly, Nx, Ny));
ReciTFWave = ReciTransWave.*ObjLens;
TFWave=ifftshift(ifft2(ReciTFWave));
Image = abs(TFWave.^2);
figure;
imagesc(x, y, Image);
colormap('gray');
axis square;