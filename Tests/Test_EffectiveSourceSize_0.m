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
% Test_EffectiveSourceSize_0.m
clc;
close all;
clear all;
%% main
% Raw perfect data:
RawData = load('GaAs110adf_1.txt');

[Ny, Nx] = size(RawData);
Lx = 3 * 3.995;
Ly = 2 * 5.65;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
Rsq = X.^2 + Y.^2;

figure;
imagesc(x, y, RawData);
colormap('gray'); axis square;
title('GaAs<110> ADF-STEM image without ess');
xlabel('x / (Angs.)'); ylabel('y / (Angs.)');

ESS = 0.8;
sigma = ESS / 2.35482 / sqrt(dx*dy);
DataConv = imgaussfilt(RawData, sigma);

figure;
imagesc(x, y, DataConv);
colormap('gray'); axis square;
title('GaAs<110> ADF-STEM image without ess');
xlabel('x / (Angs.)'); ylabel('y / (Angs.)');

figure;
plot(x, RawData(Ny / 2 + 1, : ), x, DataConv(Ny / 2 + 1, : ));
legend('ideal', 'including ess');
xlabel('x / (Angs.)');