% SCEM_Preprocessing_X_test_1_ref_test_cmp.m
clc;
clear;
close all;
%% main
refDir = 'E:\practice\TEM_SIMU_MATLAB_testdata\scem_ref_1';
testDir = 'E:\practice\TEM_SIMU_MATLAB_testdata\scem_test_1';

Nx = 1024;
Ny = 1024;
params.dfSeries = -100 : 20 : 100;
dfNum = length(params.dfSeries);

params.scanx = linspace(0, 3.9051, 10);
scanNx = length(params.scanx);
params.scany = linspace(0, 3.9051, 10);
scanNy = length(params.scany);

acceptError = 1.0e-8;
for dfIdx = 1 : dfNum
    dfDir = sprintf('df=%.4fAngs', params.dfSeries(dfIdx));
    refDfDir = fullfile(refDir, dfDir);
    testDfDir = fullfile(testDir, dfDir);
    
    for yIdx = 1 : scanNy
        for xIdx = 1 : scanNx
            refFilename = sprintf('ref_y%d_x%d.bin', yIdx, xIdx);
            refFilename = fullfile(refDfDir, refFilename);
            refData = ReadComplexBinaryFile(refFilename, [Ny, Nx], 'column');
            
            testFilename = sprintf('wave_y%d_x%d.bin', yIdx, xIdx);
            testFilename = fullfile(testDfDir, testFilename);
            testData = ReadComplexBinaryFile(testFilename, [Ny, Nx], 'column');
            
            testRefDiff = testData - refData;
            sigma = sqrt(mean(abs(testRefDiff.^2), 'all'));
            if sigma > acceptError
                disp('Unaccepted loss of accuracy');
                return;
            end
        end
    end
end