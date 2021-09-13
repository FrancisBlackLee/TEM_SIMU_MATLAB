% SCEM_X_test_2_test_data.m
% Using SrTiO3 100
clc;
clear;
close all;
%% specimen preparation:
lattConst = [3.9051, 3.9051, 0]; % [a b]
sliceDists = [1.9525, 1.9525]; % distance between each slice
expanNum = 15 * [1, 1];
% Laters: Each column for an atom
sliceA = [38,   8;...
          1,    1;...
          0,    0.5;...
          0,    0.5];
sliceB = [22,   8,      8;...
          1,    1,      1;...
          0.5,  0.5,    0;...
          0.5,  0,      0.5];

%% sampling:
Lx = expanNum(1) * lattConst(1);
Ly = expanNum(2) * lattConst(2);
Nx = 1024;
Ny = 1024;
dx = Lx / Nx;
dy = Ly / Ny;

%% STEM settings:
params.KeV = 300;
interCoeff = InteractionCoefficient(params.KeV);
wavLen = HighEnergyWavLen_X(params.KeV);

% annular upper aperture
upperInnerAngle = 21;
upperOuterAngle = 24.5;
params.upperAperture = AnnularAperture_X(Lx, Ly, Nx, Ny, wavLen,...
    upperInnerAngle, upperOuterAngle);

% circular lower aperture
lowerOuterAngle = 20.5;
params.lowerAperture = CircApert_X(Lx, Ly, Nx, Ny, wavLen, lowerOuterAngle);

params.pinholeRadii = linspace(dx, Lx / 2, 5);

% Initialize aberrations
params.Cs3 = 0;
params.Cs5 = 0;
params.df = 0;
params.dfSeries = -50 : 20 : 50;
params.scanx = linspace(0, 3.9051, 10);
scanNx = length(params.scanx);
params.scany = linspace(0, 3.9051, 10);
scanNy = length(params.scany);

%% Transmission functions:
stackNum = 20;
projPotA = MultiProjPot_conv_X(sliceA, expanNum, lattConst, Lx, Ly, Nx, Ny, 1e-5);
projPotB = MultiProjPot_conv_X(sliceB, expanNum, lattConst, Lx, Ly, Nx, Ny, 1e-5);

tfA = exp(1i * interCoeff * projPotA / 1000);
tfB = exp(1i * interCoeff * projPotB / 1000);
tfA = BandwidthLimit(tfA, Lx, Ly, Nx, Ny, 0.67);
tfB = BandwidthLimit(tfB, Lx, Ly, Nx, Ny, 0.67);
transFuncs(:, :, 1) = tfA;
transFuncs(:, :, 2) = tfB;

%% Generate referrence data:
destDir = 'E:\practice\TEM_SIMU_MATLAB_testdata\scem_test_2';
SCEM_X(Lx, Ly, params, transFuncs, sliceDists, stackNum, destDir, 'reduced');