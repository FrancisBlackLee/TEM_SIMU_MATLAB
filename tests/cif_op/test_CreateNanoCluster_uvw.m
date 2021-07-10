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
% test_CreateNanoCluster_uvw.m
clc;
clear;
close all;
%% filenames:
filename_1 = 'tests\cif_op\AlGaAs2_mp-1228891_computed.cif';
filename_2 = 'tests\cif_op\AlGaAs2_mp-1228891_conventional_standard.cif';
filename_3 = 'tests\cif_op\AlGaAs2_mp-1228891_primitive.cif';
filename_4 = 'tests\cif_op\AlGaAs2_mp-1228891_symmetrized.cif';

filename_5 = 'tests\cif_op\Si_mp-149_computed.cif';
filename_6 = 'tests\cif_op\Si_mp-149_conventional_standard.cif';
filename_7 = 'tests\cif_op\Si_mp-149_primitive.cif';
filename_8 = 'tests\cif_op\Si_mp-149_symmetrized.cif';

%% main:
crysInfo = LoadCif(filename_4);

[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
initAtomSiteMat = ExtractAtomSiteFromCrysInfo(crysInfo);
fullAtomSiteMat = AddEquivAtomSites(initAtomSiteMat);

radius = 20;
uvw = [1, 1, 1];

atomCoordMat = CreateNanoCluster_uvw(fullAtomSiteMat, cellLengths, cellAngles, uvw, radius);

figure;
scatter(atomCoordMat(3, :), atomCoordMat(4, :), 10, 'filled');
% scatter3(atomCoordMat(3, :), atomCoordMat(4, :), atomCoordMat(5, :), 'filled');
axis equal;