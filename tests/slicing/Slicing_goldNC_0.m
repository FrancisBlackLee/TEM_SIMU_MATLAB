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
%% Slicing_goldNC_0.m 2020/04/05
% 4nm gold in water
clc;
close all;
clear all;
%% main:
latt = load('GoldClusterCoord_4nm.txt');
latt = [latt(:, 1), ones(size(latt(:, 1))), latt(:, 6), latt(:, 7), latt(:, 8)]';

% figure;
% scatter3(latt(3, : ), latt(4, : ), latt(5, : ));

% move the nano cluster to the center of vacuum box
latt(3:4, :) = latt(3:4, :) + 30.0;
latt(5, : ) = latt(5, : ) + 50.0;
% figure;
% scatter3(latt(3, : ), latt(4, : ), latt(5, : ));


water = load('water_20nm.txt');
water = water';
waterAN = size(water, 2);
water = [water(1, :); ones(1, waterAN); water(2:4, :)];
% figure;
% scatter3(water(3, :), water(4, :), water(5, :));

water(:, (water(3, :)>60) | (water(4, :)>60)) = [];
water(:, water(5, :)>100) = [];
water(:, sqrt((water(3,:)-30).^2 + (water(4,:)-30).^2 + (water(5,:)-50).^2) < 20) = [];
% figure;
% scatter3(water(3, :), water(4, :), water(5, :));

sample = [latt, water];

[slice, SliceDist, ExtraSlice] = CrystalSlicing_X(sample, sample, 2.0, 100.0, 1, 0);

Lx = 60.0;
Ly = Lx;
Nx = 1024;
Ny = 1024;
CellNum = [1, 1];
LattConst = [Lx, Ly];
destDir = 'in_situ_gold_in_water';
mkDirStat = mkdir(destDir);
save(fullfile(destDir, 'SliceDist_4nmGold_in_water.txt'),...
    'SliceDist',  '-ascii', '-double', '-tabs');

wbHandle = waitbar(0, 'Generating projected potential...');
sliceNum = length(SliceDist);
for sliceIdx = 1 : sliceNum
    filename = strcat('ppj', num2str(sliceIdx), '.bin');
    filename = fullfile(destDir, filename);
    tempSlice = slice{1, sliceIdx};
    tempSlice(3, : ) = tempSlice(3, : ) / Lx;
    tempSlice(4, : ) = tempSlice(4, : ) / Ly;
    ProjPot = MultiProjPot_conv_X(tempSlice, CellNum, LattConst, Lx, Ly, Nx, Ny, 1.0e-5);
    WriteBinaryFile(filename, ProjPot);
    waitbar(sliceIdx / sliceNum, wbHandle, 'Generating projected potential...');
end

close(wbHandle);