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
% test_CellExpand_uvw.m
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
crysInfo = LoadCif(filename_8);

[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
initAtomSiteMat = ExtractAtomSiteFromCrysInfo(crysInfo);
fullAtomSiteMat = AddEquivAtomSites(initAtomSiteMat);

sideLength = 50;
uvw = [1, 1, 0];

convMat = ConversionMatrix_uvw(cellLengths, cellAngles, uvw);
viewDirection = uvw(1) * convMat(:, 1) +...
        uvw(2) * convMat(:, 2) +...
        uvw(3) * convMat(:, 3);

atomCoordMat = CellExpand_uvw(fullAtomSiteMat, cellLengths, cellAngles, uvw, sideLength);

deleteIndices = find((atomCoordMat(3, :) > sideLength / 2) |...
    (atomCoordMat(3, :) < -sideLength / 2) |...
    (atomCoordMat(4, :) > sideLength / 2) |...
    (atomCoordMat(4, :) < -sideLength / 2));
atomCoordMat(:, deleteIndices) = [];

maxSliceSpacing = 2.0;
zMax = norm(viewDirection);
lattMode = 0;
plotYN = 0;
[slice, sliceDist, extraSlice] = CrystalSlicing_X(atomCoordMat, atomCoordMat,...
    maxSliceSpacing, zMax, lattMode, plotYN);

figure;
scatter(atomCoordMat(3, :), atomCoordMat(4, :), 10, 'filled');
% scatter3(atomCoordMat(3, :), atomCoordMat(4, :), atomCoordMat(5, :), 'filled');
axis equal;