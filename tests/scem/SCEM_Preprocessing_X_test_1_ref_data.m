% SCEM_Preprocessing_X_test_1_ref_data.m
% Using SrTiO3 100
clc;
clear;
close all;
%% specimen preparation:
lattConst = [3.9051, 3.9051, 0]; % [a b]
sliceDist = [1.9525, 1.9525]; % distance between each slice
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
innerAngle = 21;
outerAngle = 24.5;
params.aperture = AnnularAperture_X(Lx, Ly, Nx, Ny, wavLen, innerAngle, outerAngle);
% Initialize aberrations
params.Cs3 = 0;
params.Cs5 = 0;
params.df = 0;
params.dfSeries = -100 : 20 : 100;
dfNum = length(params.dfSeries);
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
destDir = 'E:\practice\TEM_SIMU_MATLAB_testdata\scem_ref';
taskNum = dfNum * scanNy;
wbHandle = waitbar(0, 'scanning...');
for dfIdx = 1 : dfNum
    dfDir = sprintf('df=%.4fAngs', params.dfSeries(dfIdx));
    dfDir = fullfile(destDir, dfDir);
    status = mkdir(dfDir);
    if status == 0
        return;
    end
    
    params.df = params.dfSeries(dfIdx);
    otf = params.aperture .* ObjTransFunc_X(params, Lx, Ly, Nx, Ny);
    
    for yIdx = 1 : scanNy
        doneRatio = ((dfIdx - 1) * scanNy + yIdx - 1) / taskNum;
        wbMessage = sprintf('df: %d / %d, line: %d / %d completed',...
            dfIdx - 1, dfNum, yIdx - 1, scanNy);
        waitbar(doneRatio, wbHandle, wbMessage);
        for xIdx = 1 : scanNx
            wave = GenerateProbe_X(otf, params.scanx(xIdx),...
                params.scany(yIdx), Lx, Ly, Nx, Ny);
            wave = multislice_X(wave, params.KeV, Lx, Ly, transFuncs,...
                sliceDist, stackNum);
            
            filename = sprintf('ref_y%d_x%d.bin', yIdx, xIdx);
            filename = fullfile(dfDir, filename);
            WriteComplexBinaryFile(filename, wave, 'column');
        end
    end
end

delete(wbHandle);