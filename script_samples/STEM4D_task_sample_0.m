% STEM4D_task_sample_0.m
clc;
clear;
close all;
%% specimen preparation:
lattConsts = [3.9051, 3.9051, 3.9051]; % [a b]
sliceDists = [1.9525, 1.9525]; % distance between each slice
expanNum = [4, 4];
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
lx = expanNum(1) * lattConsts(1);
ly = expanNum(2) * lattConsts(2);
nx = 512;
ny = 512;
dx = lx / nx;
dy = ly / ny;
%% STEM settings:
outDir = './script_samples/stem4d_out';
err = mkdir(outDir);
task = STEM4DTask(outDir);
task.keV = 300;
task.lx = lx;
task.ly = ly;
wavLen = HighEnergyWavLen_X(task.keV);
task.aperture = CircApert_X(lx, ly, nx, ny, wavLen, 21.4);
scanNx = 40;
scanNy = 40;
task.scanx = linspace(0, 3.9051, scanNx);
task.scany = linspace(0, 3.9051, scanNy);

task.nBin = 2;
task.useFullAberration = true;
task.useGPU = true;

%% specimen:
interCoeff = InteractionCoefficient(task.keV);
projPotA = MultiProjPot_conv_X(sliceA, expanNum, lattConsts, lx, ly, nx, ny, 1e-5);
% figure;
% imagesc(projPotA);
% colormap('gray'); axis square;
% title('Proj.Pot. A');

projPotB = MultiProjPot_conv_X(sliceB, expanNum, lattConsts, lx, ly, nx, ny, 1e-5);
% figure;
% imagesc(projPotB);
% colormap('gray'); axis square;
% title('Proj.Pot. B');
tic;
tfA = exp(1i * interCoeff * projPotA / 1000);
tfB = exp(1i * interCoeff * projPotB / 1000);
transFuncs(:, :, 1) = tfA;
transFuncs(:, :, 2) = tfB;

task.transFuncs = gpuArray(single(transFuncs));

task.depths = 50 : 50 : 200;
task.sliceDists = sliceDists;
task.nStack = round(max(task.depths) / lattConsts(3));

%% execute task:
task.ExecuteTask();

%% check results:
binNx = nx / task.nBin;
binNy = ny / task.nBin;
detector = AnnularDetector_X(76, 200, wavLen, lx, ly, nx, ny, task.nBin);

t = 200;
pacbed = zeros(binNy, binNx);
haadf = zeros(scanNy, scanNx);
for scanIy = 1 : scanNy
    for scanIx = 1 : scanNx
        filename = ['scan_x_', num2str(scanIx - 1), '_y_',...
            num2str(scanIy - 1), '_z=', num2str(t, '%.6f'), 'A.bin'];
        cbed = ReadBinaryFile(fullfile(outDir, filename), [binNy, binNx], 'row', 'float');
        imagesc(cbed);
        colormap('gray');
        axis square;
        drawnow;

        haadf(scanIy, scanIx) = VtemlabPacbedQstem(cbed, detector);
        pacbed = pacbed + cbed;
    end
end

%% plot
figure;
imagesc(haadf);
colormap('gray');
axis equal;
axis tight;

figure;
logPacbed = log(1 + 1000 * pacbed / max(pacbed, [], 'all'));
imagesc(logPacbed);
colormap('gray');
axis square;