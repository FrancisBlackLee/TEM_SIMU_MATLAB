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
% CBED_sample_2.m: silicon [111], add thermal displacements
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
LayerA = [0.75, 0.75, 0.25, 0.25; 0.0833, 0.4167, 0.5833, 0.9167];
LayerB = [0.75, 0.25, 0.25, 0.75; 0.0833, 0.25, 0.5833, 0.75];
LayerC = [0.25, 0.75, 0.75, 0.25; 0.25, 0.4167, 0.75, 0.9167];

LayerA = SquareLattExpan(LayerA,LattConst, CellNum);
% figure;
% scatter(LayerA(1, :), LayerA(2, :));
LayerB = SquareLattExpan(LayerB,LattConst, CellNum);
% figure;
% scatter(LayerB(1, :), LayerB(2, :));
LayerC = SquareLattExpan(LayerC,LattConst, CellNum);
% figure;
% scatter(LayerC(1, :), LayerC(2, :));
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
StackNum = 10;
TransFuncs = zeros(Ny, Nx, 3*StackNum);
for i = 1 : StackNum
    TempLayer = LayerA + normrnd(0.0, 0.075, size(LayerA));
    ProjPot = ProjectedPotential(Lx, Ly, Nx, Ny, 14 * ones(1, size(TempLayer, 2)), TempLayer(1, :), TempLayer(2, :));
    tempTF = exp(1i * InterCoeff * ProjPot / 1.0e3);
    TransFuncs(:, :, 3*(i-1) + 1) = tempTF;
    
    TempLayer = LayerB + normrnd(0.0, 0.075, size(LayerB));
    ProjPot = ProjectedPotential(Lx, Ly, Nx, Ny, 14 * ones(1, size(TempLayer, 2)), TempLayer(1, :), TempLayer(2, :));
    tempTF = exp(1i * InterCoeff * ProjPot / 1.0e3);
    TransFuncs(:, :, 3*(i-1) + 2) = tempTF;
    
    TempLayer = LayerC + normrnd(0.0, 0.075, size(LayerC));
    ProjPot = ProjectedPotential(Lx, Ly, Nx, Ny, 14 * ones(1, size(TempLayer, 2)), TempLayer(1, :), TempLayer(2, :));
    tempTF = exp(1i * InterCoeff * ProjPot / 1.0e3);
    TransFuncs(:, :, 3*(i-1) + 3) = tempTF;
end

LayerDist = repmat(LayerDist, 1, StackNum);
%% CBED module
Probe = ProbeCreate(Params, 0, 0, Lx, Ly, Nx, Ny);
TransWave = multislice(Probe, WaveLength, Lx, Ly, TransFuncs, LayerDist, 1);
% Trans_Wave_Far = ifftshift(fft2(fftshift(TransWave)) * dx * dy);
Trans_Wave_Far = ifftshift(fft2(fftshift(TransWave)));
DetectInten = log(1 + 0.1 * abs(Trans_Wave_Far.^2));

% Show the detected image:
figure;
subplot(1, 2, 1);
imagesc(fx, fy, DetectInten);
colormap('gray');
axis square;
subplot(1, 2, 2);
plot(fx, DetectInten(Ny / 2 + 1, : ));

thickness = sum(LayerDist(:)) * StackNum;