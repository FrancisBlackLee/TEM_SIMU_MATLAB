% test_ReadVtlCryst.m
clc;
clear;
close all;
%% main:
[bases, unitCell, vbases, anisoU] = ReadVtlCryst('tests/vtlfiles/AlGaAs2_mp-1228891_symmetrized.cryst');
convMat = BasesToConvMat(bases);

figure;
PlotUnitCell3D(convMat, unitCell);