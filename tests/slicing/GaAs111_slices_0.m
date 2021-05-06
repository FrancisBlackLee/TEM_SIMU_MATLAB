% GaAs111_slices_0.m
clc;
close all;
clear all;
%% Load Crystal
latt = load('E:\projects\vTEMLAB_v0_testdata\AlGaAs110_TDS\GaAs Coordinates_0.txt');
atomNum = size(latt, 1);
latt = [latt(:, 1), ones(atomNum, 1), latt(:, 6), latt(:, 7), latt(:, 8)]';

zMax = max(latt(5, :)) - min(latt(5, :));

plotColor = ones(1, atomNum);
plotColor(latt(1, :) == 33) = 2;

[slice, SliceDist, ExtraSlice] = CrystalSlicing_X(latt, latt, 1.0, zMax, 1, 1, plotColor);