% test_WriteVtlCryst.m
clc;
clear;
close all;
%% main:
crysInfo = LoadCif('tests/cif_op/AlGaAs2_mp-1228891_symmetrized.cif');
[unitCell, convMat] = CrysInfoToUnitCell(crysInfo);

figure;
PlotUnitCell3D(convMat, unitCell);

bases = ConvMatToBases(convMat);
vbases.a = [1, 0, 0];
vbases.b = [0, 1, 0];
vbases.c = [0, 0, 1];
nAtom = size(unitCell, 2);
anisoU = zeros(6, nAtom);
anisoU(1 : 2, :) = 0.005;
anisoU(3, :) = 0.008;
anisoU(6, :) = 0.0025;

WriteVtlCryst('tests/vtlfiles/AlGaAs2_mp-1228891_symmetrized.cryst', bases, ...
    unitCell, vbases, anisoU);