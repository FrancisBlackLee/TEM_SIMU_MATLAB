% test_MonoProjPot_conv_vs_imtl_benchmark_1.m
clc;
clear;
close all;
%% basic settings:
Lx = 20;
Ly = Lx;
Nx = 1024;
Ny = 1024;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;

atomType = 12;
cellNum = [1, 1];
lattConst = [Lx, Ly];
testIteNum = 10;
%% benchmark:
testAtomNum = (1 : 9) * 1e2;
testNum = length(testAtomNum);
convTimeArray = zeros(1, testNum);
imtlTimeArray = zeros(1, testNum);

wbHandle = waitbar(0, 'computing...');
for testIdx = 1 : testNum
    atomNum = testAtomNum(testIdx);
    eleProp = ones(1, atomNum);
    xyCoords = zeros(2, atomNum);
    xyCoords(1, :) = Lx * rand(1, atomNum) - Lx / 2;
    scaledCoords = xyCoords / Lx;
    
    % conv method:
    tic;
    for iteIdx = 1 : testIteNum
        projPot = MonoProjPot_conv_X(atomType, eleProp, scaledCoords,...
            cellNum, lattConst, Lx, Ly, Nx, Ny);
    end
    convTimeArray(testIdx) = toc / testIteNum;
    
    % imtl method:
    tic;
    for iteIdx = 1 : testIteNum
        projPot = MonoProjPot_imtl_X(atomType, eleProp, xyCoords,...
            Lx, Ly, Nx, Ny);
    end
    imtlTimeArray(testIdx) = toc / testIteNum;
    
    waitbar(testIdx / testNum, wbHandle, 'computing...');
end
close(wbHandle);

figure;
plot(testAtomNum, convTimeArray, testAtomNum, imtlTimeArray);
legend('conv', 'imtl');