% test_PlotUnitCell.m
clc;
clear;
close all;
%% Load CIF:
filename_1 = 'tests\cif_op\AlGaAs2_mp-1228891_computed.cif';
filename_2 = 'tests\cif_op\AlGaAs2_mp-1228891_conventional_standard.cif';
filename_3 = 'tests\cif_op\AlGaAs2_mp-1228891_primitive.cif';
filename_4 = 'tests\cif_op\AlGaAs2_mp-1228891_symmetrized.cif';

filename_5 = 'tests\cif_op\Si_mp-149_computed.cif';
filename_6 = 'tests\cif_op\Si_mp-149_conventional_standard.cif';
filename_7 = 'tests\cif_op\Si_mp-149_primitive.cif';
filename_8 = 'tests\cif_op\Si_mp-149_symmetrized.cif';

crysInfo = LoadCif(filename_1);

[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
atomSites = ExtractAtomSiteFromCrysInfo(crysInfo);
initAtomNum = size(atomSites, 2);

uvw = [2, 1, 1];
convMat = ConversionMatrix_uvw(cellLengths, cellAngles, uvw);

%% plot the initial unit cell:
figure;
PlotUnitCell2D(convMat, atomSites);
title('Unit cell (2D)');
figure;
PlotUnitCell3D(convMat, atomSites);
title('Unit cell (3D)');

%% Add translated unit cell
cellNum = 4;
expandAtomSites = repmat(atomSites, 1, cellNum);

for cellIdx = 1 : cellNum
    startIdx = (cellIdx - 1) * initAtomNum + 1;
    tmpRange = startIdx : startIdx + initAtomNum - 1;
    expandAtomSites(3, tmpRange) = expandAtomSites(3, tmpRange) + cellIdx - 2;
    expandAtomSites(4, tmpRange) = expandAtomSites(4, tmpRange) + cellIdx - 2;
    expandAtomSites(5, tmpRange) = expandAtomSites(5, tmpRange) + cellIdx - 2;
end

%% plot the expanded cells:
figure;
PlotUnitCell2D(convMat, expandAtomSites);
title('Expanded cell (2D)');
figure;
PlotUnitCell3D(convMat, expandAtomSites);
title('Expanded cell (3D)');