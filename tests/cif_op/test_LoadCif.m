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
% test_LoadCif.m
clc;
clear;
close all;
%% main:
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

crysInfo = LoadCif(filename_4);

[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
atomSiteMat = ExtractAtomSiteFromCrysInfo(crysInfo);
fullAtomSiteMat = AddEquivAtomSites(atomSiteMat);

hkl = [-3, 2, 3];
convMat = ConversionMatrix_hkl(cellLengths, cellAngles, hkl);
atomCartCoord = convMat * fullAtomSiteMat(3 : 5, :);

figure;
% scatter(atomCartCoord(1, :), atomCartCoord(2, :), 'filled');
scatter3(atomCartCoord(1, :), atomCartCoord(2, :), atomCartCoord(3, :), 'filled');
axis equal;
title(['hkl = ', num2str(hkl)]);