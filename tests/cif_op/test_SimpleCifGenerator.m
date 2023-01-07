% test_SimpleCifGenerator.m
clc;
clear;
close all;
%% Read data from Ref CIF:
refDir = 'tests\cif_op';
filename1 = 'AlGaAs2_mp-1228891_computed.cif';
filename2 = 'AlGaAs2_mp-1228891_conventional_standard.cif';
filename3 = 'AlGaAs2_mp-1228891_primitive.cif';
filename4 = 'AlGaAs2_mp-1228891_symmetrized.cif';

filename5 = 'Si_mp-149_computed.cif';
filename6 = 'Si_mp-149_conventional_standard.cif';
filename7 = 'Si_mp-149_primitive.cif';
filename8 = 'Si_mp-149_symmetrized.cif';
refFilename = filename8;

crysInfo = LoadCif(fullfile(refDir, refFilename));

[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
fracCoords = ExtractAtomSiteFromCrysInfo(crysInfo);
fullFracCoords = AddEquivAtomSites(fracCoords);

uvw = [1, 1, 2];
convMat = ConversionMatrix_uvw(cellLengths, cellAngles, uvw);

%% Write simple CIF using convMat and fracCoords:
testFilename = ['test_simple_', refFilename];
SimpleCifGenerator(fullfile(refDir, testFilename), convMat, fullFracCoords);
