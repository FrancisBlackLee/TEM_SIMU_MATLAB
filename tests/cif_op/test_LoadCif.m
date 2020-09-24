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