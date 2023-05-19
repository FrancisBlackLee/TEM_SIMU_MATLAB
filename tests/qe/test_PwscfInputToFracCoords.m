% test_PwscfInputToFracCoords.m
clc;
clear;
close all;
%% main:
pwscf = ReadPwscfInput('tests/qe/hBN.1_scf.in');
[fracCoords, a1, a2, a3, ~, ~, ~] = PwscfInputToFracCoords(pwscf);

convMat = BasesToConvMat(a1, a2, a3);

figure;
PlotUnitCell3D(convMat, fracCoords);