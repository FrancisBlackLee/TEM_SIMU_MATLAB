% prototype_prism_post.m
clc;
clear;
close all;
%% sampling:
lattConst = [3.84, 5.43, 0]; % [a b]
expanNum = [4, 3];
lx = expanNum(1) * lattConst(1);
ly = expanNum(2) * lattConst(2);
nx = 512;
ny = 512;

keV = 300;
wavLen = HighEnergyWavLen_X(keV);

fx = InitFreqAxis(lx, nx);
fy = InitFreqAxis(ly, ny);

%% load data:
load("tests/prism/prism_waves.mat");

%% scanning:
scanNx = round(lattConst(1) / 0.1) + 1;
scanNy = round(lattConst(2) / 0.1) + 1;
scanx = linspace(0, lattConst(1), 39);
scany = linspace(0, lattConst(2), 55);
adf = zeros(scanNy, scanNx);
nWave = size(prismWaves, 3);
detector = AnnularDetector_X(76, 200, wavLen, lx,  ly, nx, ny);
wb = waitbar(0, 'starting...');
for scanIy = 1 : scanNy
    yp = scany(scanIy);
    for scanIx = 1 : scanNx
        xp = scanx(scanIx);
        wave = zeros(ny, nx);
        for iWave = 1 : nWave
            wave = wave + prismWaves(:, :, iWave) *...
                exp(-1i * 2 * pi * (fx(ifx(iWave)) * xp + fy(ify(iWave)) * yp));
        end
        waveI = abs(wave.^2);
        adf(scanIy, scanIx) = VtemlabPacbedQstem(waveI, detector);
    end

    waitbar(scanIy / scanNy, wb, [num2str(scanIy), '/', num2str(scanNy)]);
end
close(wb);

%% show adf:
figure;
imagesc(scanx, scany, adf);
colormap('gray');
axis equal;
axis tight;