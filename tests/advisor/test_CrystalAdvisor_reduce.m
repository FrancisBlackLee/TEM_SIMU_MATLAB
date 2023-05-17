% test_CrystalAdvisor_reduce.m
clc;
clear;
close all;
%% unit cell:
cifFilename = 'mp-20351_InP.cif';

crysInfo = LoadCif(cifFilename);

[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
atomSites = ExtractAtomSiteFromCrysInfo(crysInfo);
fullAtomSites = AddEquivAtomSites(atomSites);

uvw = [1, 1, 2];
convMat = ConversionMatrix_uvw(cellLengths, cellAngles, uvw);

%% plot the initial unit cell:
figure;
PlotUnitCell2D(convMat, fullAtomSites);
title('Initial unit cell (2D)');
figure;
PlotUnitCell3D(convMat, fullAtomSites);
title('Initial unit cell (3D)');

%% calculate the super cell and plot it:
[fracCoords, convMat, cutoff] = CrystalAdvisor(cifFilename, uvw, true);
fullFracCoords = AddEquivAtomSites(fracCoords);

figure;
PlotUnitCell2D(convMat, fullFracCoords);
title('Transformed unit cell (2D)');
figure;
PlotUnitCell3D(convMat, fullFracCoords);
title('Transformed unit cell (3D)');
