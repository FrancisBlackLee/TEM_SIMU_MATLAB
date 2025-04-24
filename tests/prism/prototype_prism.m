% prototype_prism.m
clc;
clear;
close all;
%% si [110]
lattConst = [3.84, 5.43, 0]; % [a b]
sliceDist = [1.9198, 1.9198]; % distance between each slice
expanNum = [4, 3]; % expand the unit cell by Expan_Nx = 3M and Expan_Ny = 2M, adaptive integer M
% slices: Each column stands for an atom
sliceA = [14, 14; 0, 0.5; 0, 0.75];
sliceB = [14, 14; 0, 0.5; 0.25, 0.5];

%% sampling and keV
lx = expanNum(1) * lattConst(1);
ly = expanNum(2) * lattConst(2);
nx = 512;
ny = 512;

keV = 300;
wavLen = HighEnergyWavLen_X(keV);
interCoeff = InteractionCoefficient(keV);

%% Transmission functions
% slice A projected potential:
projPotA = MultiProjPot_conv_0(sliceA, expanNum, lattConst, lx, ly, nx, ny);
% slice B projected potential:
projPotB = MultiProjPot_conv_0(sliceB, expanNum, lattConst, lx, ly, nx, ny);

% % test
% figure;
% imagesc(x, y, projPotA);
% colormap('gray');
% figure;
% imagesc(x, y, projPotB);
% colormap('gray');

tfA = exp(1i * interCoeff * projPotA / 1000);
tfB = exp(1i * interCoeff * projPotB / 1000);
tfA = BandwidthLimit(tfA, lx, ly, nx, ny, 2/3);
tfB = BandwidthLimit(tfB, lx, ly, nx, ny, 2/3);
transFuncs(:, :, 1) = tfA;
transFuncs(:, :, 2) = tfB;

%% create a set of incident beams
amax = 15;
aperture = CircApert_X(lx, ly, nx, ny, wavLen, amax);
fi = find(aperture == 1);
[ify, ifx] = ind2sub([ny, nx], fi);
nWave = length(fi);
prismWaves = 1 + 1i * ones(ny, nx, nWave);
wb = waitbar(0, 'starting...');
for iWave = 1 : nWave
    wave = zeros(ny, nx);
    wave(ify(iWave), ifx(iWave)) = 1;
    wave = ifftshift(ifft2(fftshift(wave)));
    wave = multislice(wave, wavLen, lx, ly, transFuncs, sliceDist, 10);
    wave = ifftshift(fft2(fftshift(wave)));
    prismWaves(:, :, iWave) = wave;
    waitbar(iWave / nWave, wb, [num2str(iWave), '/', num2str(nWave)]);
end
close(wb);

save("tests/prism/prism_waves.mat", "prismWaves", "ifx", "ify", "-v7.3");