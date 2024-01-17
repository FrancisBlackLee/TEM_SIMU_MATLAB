% prototype_fourier_domain_rotation.m
clc;
clear;
close all;
%% wave properties:
keV = 300;
wavLen = HighEnergyWavLen_X(keV);

%% sampling on source plane:
srcLx = 10;
srcLy = 10;
srcNx = 512;
srcNy = 512;
srcX = InitAxis(srcLx, srcNx);
srcY = InitAxis(srcLy, srcNy);
srcFx = InitFreqAxis(srcLx, srcNx);
srcFy = InitFreqAxis(srcLy, srcNy);
[srcFxMesh, srcFyMesh] = meshgrid(srcFx, srcFy);
srcFzMesh = sqrt(wavLen^-2 - srcFxMesh.^2 - srcFyMesh.^2);

%% wave on source plane:
r_srcWave = ones(srcNy, srcNx);
k_srcWave = ifftshift(fft2(fftshift(r_srcWave)));

%% sampling on reference plane:
phi = pi / 50;
refLx = srcLx;
refLy = srcLy;
refNx = 512;
refNy = 512;
refX = InitAxis(refLx, refNx);
refY = InitAxis(refLy, refNy);
refFx = InitFreqAxis(refLx, refNx);
refFy = InitFreqAxis(refLy, refNy);
[refFxMesh, refFyMesh] = meshgrid(refFx, refFy);
refFzMesh = sqrt(wavLen^-2 - refFxMesh.^2 - refFyMesh.^2);

% a1  a2  a3
% a4  a5  a6
% a7  a8  a9
invRotMat = [cos(phi), 0, sin(phi);...
             0,          1, 0;...
            -sin(phi), 0, cos(phi)];

%% wave on reference plane:
alpha = invRotMat(1, 1) * refFxMesh + invRotMat(1, 2) * refFyMesh +...
    invRotMat(1, 3) * refFzMesh;
beta = invRotMat(2, 1) * refFxMesh + invRotMat(2, 2) * refFyMesh +...
    invRotMat(2, 3) * refFzMesh;

thr = 5e-2;
k_refWave = zeros(refNy, refNx);
k_refWave(sqrt(alpha.^2 + beta.^2) <= thr) = 1;

k_refWaveI = abs(k_refWave.^2);
figure;
imagesc(refFx, refFy, k_refWaveI);
colormap('gray');
axis equal;
axis tight;
xlabel('fx (1 / \AA)', 'Interpreter', 'latex');
ylabel('fy (1 / \AA)', 'Interpreter', 'latex');
title('Intensity (k-space)');

r_refWave = ifftshift(ifft2(fftshift(k_refWave)));
r_refWaveI = abs(r_refWave.^2);
r_refWavePhase = angle(r_refWave);

figure;
subplot(1, 2, 1);
imagesc(refX, refY, r_refWaveI);
colormap('gray');
axis equal;
axis tight;
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title('Intensity');

subplot(1, 2, 2);
imagesc(refX, refY, r_refWavePhase);
colormap('gray');
axis equal;
axis tight;
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title('Phase');