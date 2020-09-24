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