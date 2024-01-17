% prototype_focused_probe_propagation_between_tilted_planes.m
clc;
clear;
close all;
%% wave properties:
keV = 300;
wavLen = HighEnergyWavLen_X(keV);
aberrs = InitObjectiveLensAberrations_X();
aberrs.C1 = 0;
c2Angle = 20;

%% sampling on reference plane:
phi = pi / 180;
lx = 15;
ly = 15;
nx = 1024;
ny = 1024;
x = InitAxis(lx, nx);
y = InitAxis(ly, ny);
fx = InitFreqAxis(lx, nx);
fy = InitFreqAxis(ly, ny);
[fxMesh, fyMesh] = meshgrid(fx, fy);
fzMesh = sqrt(wavLen^-2 - fxMesh.^2 - fyMesh.^2);

% a1  a2  a3
% a4  a5  a6
% a7  a8  a9
invRotMat = [cos(phi), 0, sin(phi);...
             0,          1, 0;...
            -sin(phi), 0, cos(phi)];

%% wave on reference plane:
alpha = invRotMat(1, 1) * fxMesh + invRotMat(1, 2) * fyMesh +...
    invRotMat(1, 3) * fzMesh;
beta = invRotMat(2, 1) * fxMesh + invRotMat(2, 2) * fyMesh +...
    invRotMat(2, 3) * fzMesh;

aperture = MeshedCircApert(alpha, beta, wavLen, c2Angle);
probe = MeshedProbe(aberrs, wavLen, aperture, 0, 0, alpha, beta);

figure;
subplot(1, 2, 1);
imagesc(x, y, abs(probe.^2));
colormap('gray');
axis equal;
axis tight;
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title('Intensity (ref plane, r-space)');

subplot(1, 2, 2);
imagesc(x, y, angle(probe));
colormap('gray');
axis equal;
axis tight;
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title('Phase (ref plane, r-space)');

%% propagation:
propDist = 100;
wave = MeshedFresnelProp(probe, fxMesh, fyMesh, wavLen, propDist, -phi * 1e3, 0);
waveI = abs(wave.^2);

figure;
imagesc(x, y, waveI);
colormap('gray');
axis equal;
axis tight;
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title(['Intensity z = ', num2str(propDist), ' \AA'], 'Interpreter', 'latex');