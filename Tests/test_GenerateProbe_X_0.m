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
%% test_GenerateProbe_X_0.m
clc;
close all;
clear all;
%% Sampling and STEM parameters:
Lx = 20;
Ly = 20;
Nx = 512;
Ny = 512;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;

params.KeV = 300;
params.Cs3 = 0;
params.Cs5 = 0;
params.df = 0;

% additional paramters for old interfaces
params.Cs = 0;
params.amax = 23;

NA = 23; % numerical aperture in mrad

WavLen = 12.3986 / sqrt((2 * 511.0 + params.KeV) * params.KeV);

OTF = CircApert_X(Lx, Ly, Nx, Ny, WavLen, NA) .* ObjTransFunc_X(params, Lx, Ly, Nx, Ny);

% draw probe profile:
NewProbeSample = GenerateProbe_X(OTF, 0, 0, Lx, Ly, Nx, Ny);
NewSampleInten = abs(NewProbeSample.^2);

OldProbeSample = ProbeCreate(params, 0, 0, Lx, Ly, Nx, Ny);
OldSampleInten = abs(OldProbeSample.^2);

figure;
plot(x, NewSampleInten(Ny / 2 + 1, : ), 'r-.', x, OldSampleInten(Ny / 2 + 1, : ) - 1, 'b-');
legend('New', 'Old - 1');

ScanCoordX = linspace(-7, 7, 32);
ScanCoordY = linspace(-8, 6, 32);
figure;
for iy = 1 : length(ScanCoordX)
    for ix = 1 : length(ScanCoordY)
        TempProbe = GenerateProbe_X(OTF, ScanCoordX(ix), ScanCoordY(iy), Lx, Ly, Nx, Ny);
        TempProbeInten = abs(TempProbe.^2);
        imagesc(x, y, TempProbeInten); colormap('gray'); axis square;
        drawnow;
    end
end