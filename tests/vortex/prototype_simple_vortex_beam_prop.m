% prototype_simple_vortex_beam_prop.m
clc;
clear;
close all;
%% sampling:
lx = 64 * 4;
ly = 64 * 4;
nx = 4096;
ny = 4096;

%% vortex beam:
params = InitProbeParams_X();
params.KeV = 200;
params.Lx = lx;
params.Ly = ly;
params.Nx = nx;
params.Ny = ny;


wavLen = HighEnergyWavLen_X(params.KeV);
x = InitAxis(lx, nx);
y = InitAxis(ly, ny);
kx = InitFreqAxis(lx, nx);
ky = InitFreqAxis(ly, ny);
[KX, KY] = meshgrid(kx, ky);
l = 5;
vortexMask = exp(1i * l * atan2(KY, KX));

params.aperture = CircApert_X(lx, ly, nx, ny, wavLen, 1.5) .* vortexMask;
figure;
imagesc(angle(params.aperture));
axis square;

probe = GenerateProbe_X(params);
probeI = abs(probe.^2);
probePhase = angle(probe);

figure;
subplot(1, 2, 1);
imagesc(x, y, probeI);
axis square;

subplot(1, 2, 2);
imagesc(x, y, probePhase);
axis square;

%% propagation:
figure('units','normalized','outerposition',[0 0 1 1]);
for z = 1000 : 1000 : 10000
    wave = propTF_1(probe, lx, ly, wavLen, z);
    waveI = abs(wave.^2);
    wavePhase = angle(wave);
    subplot(1, 2, 1);
    imagesc(x, y, waveI);
    axis square;
    title(['z = ', num2str(z), ' \AA'], 'Interpreter', 'latex');
    subplot(1, 2, 2);
    imagesc(x, y, wavePhase);
    axis square;
    drawnow;
    pause(0.1);
end