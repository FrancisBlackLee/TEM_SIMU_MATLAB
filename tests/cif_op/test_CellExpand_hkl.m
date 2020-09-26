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
% test_CellExpand_hkl.m
clc;
clear;
close all;
%% filenames:
filename_1 = ['E:\practice\crystallography\cif\cif\@cif\',...
    'AlGaAs2_mp-1228891_computed.cif'];

filename_2 = ['E:\practice\crystallography\cif\cif\@cif\',...
    'AlGaAs2_mp-1228891_conventional_standard.cif'];

filename_3 = ['E:\practice\crystallography\cif\cif\@cif\',...
    'AlGaAs2_mp-1228891_primitive.cif'];

filename_4 = ['E:\practice\crystallography\cif\cif\@cif\',...
    'AlGaAs2_mp-1228891_symmetrized.cif'];

filename_5 = ['E:\practice\crystallography\cif\cif\@cif\',...
    'Langbeinite (ortho).cif'];

%% main:
crysInfo = LoadCif(filename_2);

[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
initAtomSiteMat = ExtractAtomSiteFromCrysInfo(crysInfo);
fullAtomSiteMat = AddEquivAtomSites(initAtomSiteMat);

sideLength = 100;
hkl = [0, 2, 1];

atomCoordMat = CellExpand_hkl(fullAtomSiteMat, cellLengths, cellAngles, hkl, sideLength);

deleteIndices = find((atomCoordMat(3, :) > sideLength / 2) |...
    (atomCoordMat(3, :) < -sideLength / 2) |...
    (atomCoordMat(4, :) > sideLength / 2) |...
    (atomCoordMat(4, :) < -sideLength / 2));
atomCoordMat(:, deleteIndices) = [];

figure;
scatter(atomCoordMat(3, :), atomCoordMat(4, :), 'filled');
% scatter3(atomCoordMat(3, :), atomCoordMat(4, :), atomCoordMat(5, :), 'filled');
axis equal;