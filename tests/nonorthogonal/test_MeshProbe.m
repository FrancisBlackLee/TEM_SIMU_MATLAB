% test_MeshProbe.m
clc;
clear;
close all;
%% common settings:
aberrs = InitObjectiveLensAberrations_X();
aberrs.C1 = 100;
keV = 300;
wavLen = HighEnergyWavLen_X(keV);
numApert = 30;

%% nonorthogonal mesh:
a = 2;
a1 = a * [1/2, sqrt(3)/2, 0.0];
a2 = a * [1/2, -sqrt(3)/2, 0.0];
na1 = 20;
na2 = 20;
n1 = 1024;
n2 = 1024;
[xMesh, yMesh, fxMesh, fyMesh] = InitNonorthoMesh2D(a1, a2, na1, na2, n1, n2);

aperture = MeshedCircApert(fxMesh, fyMesh, wavLen, numApert);

figure;
mesh(fxMesh, fyMesh, aperture);
xlabel('$f_x$ (1 / \AA)', 'Interpreter', 'latex');
ylabel('$f_y$ (1 / \AA)', 'Interpreter', 'latex');
title('Nonorthogonal mesh: aperture');
axis equal;
axis tight;

otfPhase = aperture .* MeshedAberrPhase(aberrs, wavLen, fxMesh, fyMesh);
figure;
mesh(fxMesh, fyMesh, otfPhase);
xlabel('$f_x$ (1 / \AA)', 'Interpreter', 'latex');
ylabel('$f_y$ (1 / \AA)', 'Interpreter', 'latex');
title('Nonorthogonal mesh: OTF phase');
axis equal;
axis tight;

probe = MeshedProbe(aberrs, wavLen, aperture, 0.0, 0.0, fxMesh, fyMesh);
probeI = abs(probe.^2);
figure;
mesh(xMesh, yMesh, probeI);
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title('Nonorthogonal Mesh: probe intensity');
axis equal;
axis tight;

%% orthogonal mesh
lx = 40;
ly = 40;
nx = 1024;
ny = 1024;
x = InitAxis(lx, nx);
y = InitAxis(ly, ny);
fx = InitFreqAxis(lx, nx);
fy = InitFreqAxis(ly, ny);

aperture = CircApert_X(lx, ly, nx, ny, wavLen, numApert);
figure;
mesh(fx, fy, aperture);
xlabel('$f_x$ (1 / \AA)', 'Interpreter', 'latex');
ylabel('$f_y$ (1 / \AA)', 'Interpreter', 'latex');
title('Orthogonal mesh: aperture');
axis equal;
axis tight;

otfPhase = aperture .* AberrationPhaseShift_X(aberrs, wavLen, lx, ly, nx, ny);
figure;
mesh(fx, fy, otfPhase);
xlabel('$f_x$ (1 / \AA)', 'Interpreter', 'latex');
ylabel('$f_y$ (1 / \AA)', 'Interpreter', 'latex');
title('Orthogonal mesh: OTF phase');
axis equal;
axis tight;

params = InitProbeParams_X();
params.KeV = 300;
params.type = 'full';
params.aberration = aberrs;
params.aperture = aperture;
params.Lx = lx;
params.Ly = ly;
params.Nx = nx;
params.Ny = ny;

probe = GenerateProbe_X(params);
probeI = abs(probe.^2);
figure;
mesh(x, y, probeI);
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title('Orthogonal Mesh: probe intensity');
axis equal;
axis tight;