% test_ReduceUnitCell.m
clc;
clear;
close all;
%% init transform
cifFilename_1 = 'tests/advisor/Si_mp-149_symmetrized.cif';
cifFilename_2 = 'tests/advisor/mp-20351_InP.cif';
uvw = [1, 1, 0];
[fracCoords, convMat, cutoff] = CrystalAdvisor(cifFilename_1, uvw);
fullFracCoords = AddEquivAtomSites(fracCoords);

figure;
PlotUnitCell3D(convMat, fullFracCoords);
title('Transformed unit cell (3D)');

%% reduce unit cell
[reducedUnitCell, newConvMat] = ReduceUnitCell(fullFracCoords, convMat, 16);

figure;
PlotUnitCell3D(newConvMat, reducedUnitCell);
title('Reduced unit cell (3D)');