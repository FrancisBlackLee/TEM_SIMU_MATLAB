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
% PhaseErrorTest_0.m -- test MultiAberrPhaseError_X.m
clc;
close all;
clear all;
%       Aberr{1} = [C10, C12];
%       Aberr{2} = [C21, C23];
%       Aberr{3} = [C30, C32, C34];
%       Aberr{4} = [C41, C43, C45];
%       Aberr{5} = [C50, C52, C54, C56];
%% Parameter setting:
% unit of aberrations is angstrom
Aberr{1} = [0, 0];
Aberr{2} = [0, 1e4];
Aberr{3} = [0, 0, 0];
Aberr{4} = [0, 0, 0];
Aberr{5} = [0, 0, 0, 0];
Lx = 12;
Ly = 12;
Nx = 1024;
Ny = 1024;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;

KeV = 300;
WavLen = 12.3986 / sqrt((2 * 511.0 + KeV) * KeV);  %wavelength

PhaseError = MultiAberrPhaseError_X(Aberr, WavLen, Lx, Ly, Nx, Ny);

figure;
imagesc(fx, fy, PhaseError);
colormap('gray'); axis square;

amax = 20e-3; % rad;
[FX, FY] = meshgrid(fx, fy);
FreqSquare = FX.^2 + FY.^2;
Aperture = ones(size(FX));
Aperture(FreqSquare >= (sin(amax) / WavLen)^2) = 0;
Aberr_TF = exp(1i * PhaseError) .* Aperture;

xp = 0; yp = 0;
Probe = ifft2(fftshift(Aberr_TF .* exp(-1i * 2 * pi * (FX * xp + FY * yp)))) / (dx * dy);
Probe = ifftshift(Probe);
NormCoeff = sqrt(sum(sum(abs(Probe.^2))) * dx * dy);
Probe = Probe / NormCoeff;

figure;
imagesc(x, y, abs(Probe.^2));
colormap('gray'); axis square;