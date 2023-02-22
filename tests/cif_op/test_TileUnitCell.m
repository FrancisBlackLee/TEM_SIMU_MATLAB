% test_TileUnitCell.m
clc;
clear;
close all;
%% main:
crysInfo = LoadCif('Si_mp-149_symmetrized.cif');
[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
initAtomSites = ExtractAtomSiteFromCrysInfo(crysInfo);

tiles = [2, 3, 4];
tiledAtomSites = TileUnitCell(initAtomSites, tiles, 1.0e-6);

convMat = ConversionMatrix(cellLengths, cellAngles);

figure;
PlotUnitCell2D(convMat, tiledAtomSites);

figure;
PlotUnitCell3D(convMat, tiledAtomSites);