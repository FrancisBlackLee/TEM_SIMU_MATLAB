% SCEM_X_test_2_ref_test_cmp.m
clc;
clear;
close all;
%% main
refDir = 'E:\practice\TEM_SIMU_MATLAB_testdata\scem_ref_2';
testDir = 'E:\practice\TEM_SIMU_MATLAB_testdata\scem_test_2';

Nx = 1024;
Ny = 1024;
params.dfSeries = -50 : 20 : 50;
dfNum = length(params.dfSeries);

acceptError = 1.0e-8;
for dfIdx = 1 : dfNum
    dfDir = sprintf('df=%.4fAngs', params.dfSeries(dfIdx));
    refDfDir = fullfile(refDir, dfDir);
    testDfDir = fullfile(testDir, dfDir);
    
    refDataFilename = 'scem_images.mat';
    refDataFilename = fullfile(refDfDir, refDataFilename);
    refData = load(refDataFilename);
    
    testDataFilename = 'scem_images.mat';
    testDataFilename = fullfile(testDfDir, testDataFilename);
    testData = load(testDataFilename);
    
    testRefDiff = testData.scemImg - refData.scemImg;
    sigma = sqrt(mean(testRefDiff.^2, 'all'));
    if sigma > acceptError
        disp('Unaccepted loss of accuracy');
        return;
    end
end