% prototype_focused_probe_on_tilted_plane.m
clc;
clear;
close all;
%% wave properties:
keV = 300;
wavLen = HighEnergyWavLen_X(keV);
aberrs = InitObjectiveLensAberrations_X();
aberrs.C1 = 0;
c2Angle = 21.4;

%% sampling on source plane:
srcLx = 15;
srcLy = 15;
srcNx = 1024;
srcNy = 1024;
srcX = InitAxis(srcLx, srcNx);
srcY = InitAxis(srcLy, srcNy);
srcFx = InitFreqAxis(srcLx, srcNx);
srcFy = InitFreqAxis(srcLy, srcNy);
[srcFxMesh, srcFyMesh] = meshgrid(srcFx, srcFy);
srcFzMesh = sqrt(wavLen^-2 - srcFxMesh.^2 - srcFyMesh.^2);

%% wave on source plane:
srcAperture = MeshedCircApert(srcFxMesh, srcFyMesh, wavLen, c2Angle);
srcProbe = MeshedProbe(aberrs, wavLen, srcAperture, 0, 0, srcFxMesh, srcFyMesh);
figure;
subplot(1, 2, 1);
imagesc(srcX, srcY, abs(srcProbe.^2));
colormap('gray');
axis equal;
axis tight;
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title('Intensity (src plane, r-space)');

subplot(1, 2, 2);
imagesc(srcX, srcY, angle(srcProbe));
colormap('gray');
axis equal;
axis tight;
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title('Phase (src plane, r-space)');

%% sampling on reference plane:
phi = pi / 5;
refLx = srcLx;
refLy = srcLy;
refNx = srcNx;
refNy = srcNy;
refX = InitAxis(refLx, refNx);
refY = InitAxis(refLy, refNy);
refFx = InitFreqAxis(refLx, refNx);
refFy = InitFreqAxis(refLy, refNy);
[refFxMesh, refFyMesh] = meshgrid(refFx, refFy);
refFrMesh = sqrt(refFxMesh.^2 + refFyMesh.^2);
refFzMesh = sqrt(wavLen^-2 - (sin(refFrMesh * wavLen) / wavLen).^2);

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

refAperture = MeshedCircApert(alpha, beta, wavLen, c2Angle);
refProbe = MeshedProbe(aberrs, wavLen, refAperture, 2, 0, alpha, beta);

figure;
subplot(1, 2, 1);
imagesc(refX, refY, abs(refProbe.^2));
colormap('gray');
axis equal;
axis tight;
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title('Intensity (ref plane, r-space)');

subplot(1, 2, 2);
imagesc(refX, refY, angle(refProbe));
colormap('gray');
axis equal;
axis tight;
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title('Phase (ref plane, r-space)');