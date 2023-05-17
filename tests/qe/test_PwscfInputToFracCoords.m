% test_PwscfInputToFracCoords.m
clc;
clear;
close all;
%% main:
pwscf = ReadPwscfInput('tests/qe/hBN.1_scf.in');
[fracCoords, a1, a2, a3, ~, ~, ~] = PwscfInputToFracCoords(pwscf);

bases.a = a1;
bases.b = a2;
bases.c = a3;
convMat = BasesToConvMat(bases);

figure;
PlotUnitCell3D(convMat, fracCoords);