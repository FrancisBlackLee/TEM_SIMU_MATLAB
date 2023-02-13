% test_tilt_fresnel_propagation.m
clc;
clear;
close all;
%% main:
Lx = 32;
Ly = 32;
Nx = 1024;
Ny = 1024;
x = InitAxis(Lx, Nx);
y = InitAxis(Ly, Ny);

aberrations = InitObjectiveLensAberrations_X();
wavLen = HighEnergyWavLen_X(300);
otfPhase = -AberrationPhaseShift_X(aberrations, wavLen, Lx, Ly, Nx, Ny);
otf = exp(1i * otfPhase) .* CircApert_X(Lx, Ly, Nx, Ny, wavLen, 21.4);
probe = GenerateProbe_X(otf, 0.0, 0.0, Lx, Ly, Nx, Ny);

dist = 1000;
tiltX = 1;
tiltY = 1;
propKer = FresnelPropKernel_X(Lx, Ly, Nx, Ny, wavLen, dist, tiltX, tiltY);
wave = ifftshift(FresnelProp_X(fftshift(probe), fftshift(propKer)));
waveI = abs(wave.^2);

figure;
imagesc(x, y, waveI);
colormap('gray');
axis square;
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');