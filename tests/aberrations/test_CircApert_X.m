% test_CircApert_X.m
clc;
clear;
close all;
%% sampling:
Lx = 64;
Ly = 64;
Nx = 2048;
Ny = 2048;

x = InitAxis(Lx, Nx);
y = InitAxis(Ly, Ny);

%% STEM:
wavLen = HighEnergyWavLen_X(300);
pol = 2;
azi = 0;
numApert = 21.4;
aperture = CircApert_X(Lx, Ly, Nx, Ny, wavLen, numApert, pol, azi);

figure;
imagesc(x, y, aperture);
axis square;
colormap('gray');
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');

probe = ifftshift(fft2(fftshift(aperture)));
probeI = abs(probe.^2);
figure;
imagesc(x, y, probeI);
axis square;
colormap('gray');
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');

z = 1000;
farProbe = propTF_1(probe, Lx, Ly, wavLen, z, 2/3);
farProbeI = abs(farProbe.^2);
figure;
imagesc(x, y, farProbeI);
axis square;
colormap('gray');
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');