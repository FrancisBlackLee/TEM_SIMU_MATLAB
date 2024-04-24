% test_WriteEjkXyz.m
clc;
clear;
close all;
%% main:
crysInfo = LoadCif('tests/cif_op/Si_mp-149_symmetrized.cif');
[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
convMat = ConversionMatrix(cellLengths, cellAngles);
initAtomSites = ExtractAtomSiteFromCrysInfo(crysInfo);
atomSites = RemoveSymmetricAtoms(initAtomSites);
atomSites(3 : 5, :) = convMat * atomSites(3 : 5, :);
wobbles = 0.08 * ones(1, size(atomSites, 2));
WriteEjkXyz('tests/cif_op/test_si100.xyz', cellLengths, atomSites, wobbles);