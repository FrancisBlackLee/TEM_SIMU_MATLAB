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
% ProbeTest_0.m -- test the normalization of the probe -- whether it should
% be the sum over the plane or the maximum of the signal
clc;
close all;
clear all;
%% Sampling settings:
Lx = 20;
Ly = 20;
Nx = 1024;
Ny = 1024;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;

Params.KeV = 200;
Params.Cs = 1.3;
Params.df = 700;
Params.amax = 10.37;

ProbeNS = ProbeCreate(Params, 0, 0, Lx, Ly, Nx, Ny);
% ProbeNSreal = real(ProbeNS);
% ProbeNSimag = imag(ProbeNS);
ProbeNSinten = abs(ProbeNS.^2);
figure;
subplot(1, 2, 1);
plot(x, ProbeNSreal(Ny / 2 + 1, : ), 'r-', x, ProbeNSimag(Ny / 2 + 1, : ), 'b--');
legend('real', 'imag');
subplot(1, 2, 2);
plot(x, ProbeNSinten(Ny / 2 + 1, : ));
