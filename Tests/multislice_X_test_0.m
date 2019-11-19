% multislice_X_test_0.m -- multislice_X.m loading file test:
clc;
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ProjPotDir = 'E:\practice\wasted\Posi_Pot_Img_1';
Lx = 30.8242;
Ly = 30.7420;
Nx = 512;
Ny = 512;
KeV = 300;
InciWave = ones(Ny, Nx);
SliceDist = [1.6410, 2.1053, 1.7477, 2.2120];
ExitWave = multislice_X(InciWave, KeV, Lx, Ly, 'files', SliceDist, 1, ProjPotDir, '*.txt');