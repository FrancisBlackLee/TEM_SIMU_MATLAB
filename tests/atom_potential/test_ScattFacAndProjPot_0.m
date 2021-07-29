% test_ScattFacAndProjPot_0.m
clc;
clear;
close all;
%% main:
L = 36;
N = 512;
d = L / N;

atomType = 30;
refProjPot = ProjectedPotential(L, L, N, N, atomType, 0, 0);

f = InitFreqAxis(L, N);
[fMeshX, fMeshY] = meshgrid(f);
fMeshR = sqrt(fMeshX.^2 + fMeshY.^2);

a = 0.529; % Bohr radius in angstrom
e = 14.4; % elemental charge in volt - angstrom
scaleCoeff = 2 * pi * e * a;
scattFacMesh = ScatteringFactor(atomType, fMeshR);
testProjPot = real(ifftshift(ifft2(fftshift(scattFacMesh)))) / (d^2) * scaleCoeff;

x = InitAxis(L, N);
figure;
plot(x, refProjPot(N / 2 + 1, :), x, testProjPot(N / 2 + 1, :));
xlabel('$ x(\AA) $', 'Interpreter', 'latex');
ylabel('$ V_z (V \cdot \AA) $', 'Interpreter', 'latex');
legend('ref', 'test');