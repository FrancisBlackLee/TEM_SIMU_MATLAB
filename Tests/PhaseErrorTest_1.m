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
% PhaseErrorTest_1.m -- test MultiAberrPhaseError_X.m
clc;
close all;
clear all;
%       Aberr{1} = [C10, C12];
%       Aberr{2} = [C21, C23];
%       Aberr{3} = [C30, C32, C34];
%       Aberr{4} = [C41, C43, C45];
%       Aberr{5} = [C50, C52, C54, C56];
%       C10  C1  Defocus
%       C12  A1  2-Fold astigmatism
%       C21  B2  Axial coma
%       C23  A2  3-Fold astigmatism
%       C30  C3  3rd order spherical aberration
%       C32  S3  Axial star aberration
%       C34  A3  4-Fold astigmatism
%       C41  B4  4th order axial coma
%       C43  D4  3-Lobe aberration
%       C45  A4  5-Fold astigmatism
%       C50  C5  5th order spherical aberration
%       C52  S5  5th order axial star
%       C54  R5  5th order rosette
%       C56  A5  6-Fold astigmatism
%% Parameter setting:
% unit of aberrations is angstrom
aberr{1} = [-2350, 3.865];
aberr{2} = [264.7, 364.6];
aberr{3} = [1742, 6092, 1.08e7];
aberr{4} = [0, 0, 25.12e7];
aberr{5} = [0, 0, 0, 0];

dx = 0.1407; % Angstrome, expand reciprocal space to about 70 mrad.
dy = dx;
Nx = 512;
Ny = 512;
Lx = Nx * dx;
Ly = Ny * dy;

x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;

KeV = 300;
WavLen = 12.3986 / sqrt((2 * 511.0 + KeV) * KeV);  %wavelength

PhaseError = MultiAberrPhaseError_X(aberr, WavLen, Lx, Ly, Nx, Ny);
WrapPhaseError = wrapTo2Pi(PhaseError);

imagesc(fx, fy, WrapPhaseError);
colormap('gray'); axis square;