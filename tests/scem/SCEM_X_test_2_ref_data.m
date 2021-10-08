% SCEM_X_test_2_ref_data.m
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

% generate mesh for real-space coordinates for determining pinholes:
xAxis = InitAxis(Lx, Nx);
yAxis = InitAxis(Ly, Ny);
[xMesh, yMesh] = meshgrid(xAxis, yAxis);

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
pinholeNum = length(params.pinholeRadii);

% Initialize aberrations
params.Cs3 = 0;
params.Cs5 = 0;
params.df = 0;
params.dfSeries = -50 : 20 : 50;
dfNum = length(params.dfSeries);
params.scanx = linspace(0, 3.9051, 10);
scanNx = length(params.scanx);
params.scany = linspace(0, 3.9051, 10);
scanNy = length(params.scany);

%% Transmission functions:
stackNum = 20;
specimenThickness = stackNum * sum(sliceDist);
projPotA = MultiProjPot_conv_X(sliceA, expanNum, lattConst, Lx, Ly, Nx, Ny, 1e-5);
projPotB = MultiProjPot_conv_X(sliceB, expanNum, lattConst, Lx, Ly, Nx, Ny, 1e-5);

tfA = exp(1i * interCoeff * projPotA / 1000);
tfB = exp(1i * interCoeff * projPotB / 1000);
tfA = BandwidthLimit(tfA, Lx, Ly, Nx, Ny, 0.67);
tfB = BandwidthLimit(tfB, Lx, Ly, Nx, Ny, 0.67);
transFuncs(:, :, 1) = tfA;
transFuncs(:, :, 2) = tfB;

%% Generate referrence data:
destDir = 'E:\practice\TEM_SIMU_MATLAB_testdata\scem_ref_2';
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
    secondPropDist = -params.df - specimenThickness;
    secondPropKernel = FresnelPropKernel_X(Lx, Ly, Nx, Ny, wavLen, secondPropDist);
    otf = params.upperAperture .* ObjTransFunc_X(params, Lx, Ly, Nx, Ny);
    
    scemImg = zeros(scanNy, scanNx, pinholeNum);
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
            wave = ifftshift(fft2(fftshift(wave)));
            wave = secondPropKernel .* wave;
            wave = params.lowerAperture .* wave;
            wave = ifftshift(ifft2(fftshift(wave)));
            
            waveI = abs(wave.^2);
            
            % calculate relative position of the pinhole center to the
            % probe center: relative x
            relativeXMesh = xMesh - params.scanx(xIdx);
            % alter to satisfy periodic boundary condition
            alterIndices = find(relativeXMesh < -Lx / 2.0);
            relativeXMesh(alterIndices) = relativeXMesh(alterIndices) + Lx;

            alterIndices = find(relativeXMesh > Lx / 2.0);
            relativeXMesh(alterIndices) = relativeXMesh(alterIndices) - Lx;

            % relative y
            relativeYMesh = yMesh - params.scany(yIdx);
            % alter to satisfy periodic boundary condition
            alterIndices = find(relativeYMesh < -Ly / 2.0);
            relativeYMesh(alterIndices) = relativeYMesh(alterIndices) + Ly;

            alterIndices = find(relativeYMesh > Ly / 2.0);
            relativeYMesh(alterIndices) = relativeYMesh(alterIndices) - Ly;

            % relative distance
            rMesh = sqrt(relativeXMesh.^2 + relativeYMesh.^2);
            
            for pinholeIdx = 1 : pinholeNum
                pinhole = (rMesh < params.pinholeRadii(pinholeIdx));
                scemImg(yIdx, xIdx, pinholeIdx) = sum(waveI .* pinhole, 'all');
            end
        end
    end
    
    filename = 'scem_images.mat';
    filename = fullfile(dfDir, filename);
    save(filename, 'scemImg');
end

delete(wbHandle);