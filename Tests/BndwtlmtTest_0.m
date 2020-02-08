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
% BndwtlmtTest_0.m -- test the effect of bandwidth limit on projected
% potential and transmission function:
clc;
close all;
clear all;
%% Sampling settings and loading LPCMO<101> ProjPot:
LattConst = [7.7018, 7.6855];
CellNum = [4, 4];
Lx = CellNum(1) * LattConst(1);
Ly = CellNum(2) * LattConst(2);
Nx = 512;
Ny = 512;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;

Pot1 = load('E:\Group Communication\Wei_Cao\SliceProjPot_1.txt');
Pot2 = load('E:\Group Communication\Wei_Cao\SliceProjPot_2.txt');
Pot3 = load('E:\Group Communication\Wei_Cao\SliceProjPot_3.txt');
Pot4 = load('E:\Group Communication\Wei_Cao\SliceProjPot_4.txt');

Pot = Pot1 + Pot2 + Pot3 + Pot4;

% %% Apply bandwidth limit to Pot.'s:
% Pot1p = BandwidthLimit(Pot1, Lx, Ly, Nx, Ny, 1);
% Pot2p = BandwidthLimit(Pot2, Lx, Ly, Nx, Ny, 1);
% Pot3p = BandwidthLimit(Pot3, Lx, Ly, Nx, Ny, 1);
% Pot4p = BandwidthLimit(Pot4, Lx, Ly, Nx, Ny, 1);
% 
% Potp = Pot1p + Pot2p + Pot3p + Pot4p;
% 
% figure;
% subplot(1, 2, 1);
% imagesc(x, y, Pot);
% colormap('gray'); axis square;
% subplot(1, 2, 2);
% imagesc(x, y, Potp);
% colormap('gray'); axis square;
% 
% figure;
% plot(x, Pot(Ny / 2 + 1, : ), 'r-', x, Potp(Ny / 2 + 1, : ), 'b--');
% legend('unlimited', 'limited');

%% Apply bandwidth limit to TF's only:
TestDataDir = 'E:\Topics and Practice on Electron Microscopy\Projects\TEM_SIMU_MATLAB\Tests\BndwtlmtTestData\';

TF_1 = exp(1i * 0.6 * Pot1 / 1000);
TF_2 = exp(1i * 0.6 * Pot2 / 1000);
TF_3 = exp(1i * 0.6 * Pot3 / 1000);
TF_4 = exp(1i * 0.6 * Pot4 / 1000);

Phase1 = 0.6 * Pot1 / 1000; save(strcat(TestDataDir, 'Ori1.txt'), 'Phase1', '-ascii', '-double', '-tabs');
Phase2 = 0.6 * Pot2 / 1000; save(strcat(TestDataDir, 'Ori2.txt'), 'Phase2', '-ascii', '-double', '-tabs');
Phase3 = 0.6 * Pot3 / 1000; save(strcat(TestDataDir, 'Ori3.txt'), 'Phase3', '-ascii', '-double', '-tabs');
Phase4 = 0.6 * Pot4 / 1000; save(strcat(TestDataDir, 'Ori4.txt'), 'Phase4', '-ascii', '-double', '-tabs');

TF_1 = BandwidthLimit(TF_1, Lx, Ly, Nx, Ny, 2/3);
TF_2 = BandwidthLimit(TF_2, Lx, Ly, Nx, Ny, 2/3);
TF_3 = BandwidthLimit(TF_3, Lx, Ly, Nx, Ny, 2/3);
TF_4 = BandwidthLimit(TF_4, Lx, Ly, Nx, Ny, 2/3);

Phase1p = unwrap(angle(TF_1)); save(strcat(TestDataDir, 'BWL1.txt'), 'Phase1p', '-ascii', '-double', '-tabs');
Phase2p = unwrap(angle(TF_2)); save(strcat(TestDataDir, 'BWL2.txt'), 'Phase2p', '-ascii', '-double', '-tabs');
Phase3p = unwrap(angle(TF_3)); save(strcat(TestDataDir, 'BWL3.txt'), 'Phase3p', '-ascii', '-double', '-tabs');
Phase4p = unwrap(angle(TF_4)); save(strcat(TestDataDir, 'BWL4.txt'), 'Phase4p', '-ascii', '-double', '-tabs');