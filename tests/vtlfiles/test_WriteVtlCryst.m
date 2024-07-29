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
nAtom = size(unitCell, 2);
wobbles = rand(3, nAtom);

WriteVtlCryst('tests/vtlfiles/AlGaAs2_mp-1228891_symmetrized.cryst', bases, ...
    unitCell, wobbles);