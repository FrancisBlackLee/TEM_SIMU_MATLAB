% prototype_vortex_beam_binary_mask.m
clc;
clear;
close all;
%% sampling:
keV = 300;
wavLen = HighEnergyWavLen_X(keV);
lx = 50e4;
ly = 50e4;
nx = 1024;
ny = 1024;
x = InitAxis(lx, nx);
y = InitAxis(ly, ny);
[X, Y] = meshgrid(x, y);
R = sqrt(X.^2 + Y.^2);
kx = InitFreqAxis(lx, nx);
ky = InitFreqAxis(ly, ny);
[KX, KY] = meshgrid(kx, ky);
azi = atan2(Y, X);

%% binary mask
l = 1;
rMax = 5e4;
kPerp = 4 / rMax;
waveT = exp(1i * (l * azi - pi / 2));
figure;
imagesc(angle(waveT));
axis square;
colorbar;

waveRef = exp(1i * 2 * pi * kPerp * X);
holoI = 0.25 * abs((waveT + waveRef).^2);

figure;
imagesc(holoI);
axis square;
colorbar;

mask = (holoI > 0.5) & (R < rMax);

figure;
imshow(mask);

%% generate a probe:
probe = ifftshift(fft2(fftshift(mask)));
probeI = abs(probe.^2);
probePhase = angle(probe);

figure;
subplot(1, 2, 1);
imagesc(probeI);
axis square;

subplot(1, 2, 2);
imagesc(probePhase);
axis square;