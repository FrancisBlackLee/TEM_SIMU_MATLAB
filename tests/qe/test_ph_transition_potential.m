% test_ph_transition_potential.m
clc;
clear;
close all;
%% load pwscfin
pwscf = ReadPwscfInput('tests/qe/si_cubic.1_scf.in');
[unitCell, a1, a2, a3, b1, b2, b3] = PwscfInputToFracCoords(pwscf);
pwscf = ExpandPwscfAtomMass(pwscf);
nks = [5, 5, 5];
nCell = nks(1) * nks(2) * nks(3);
superCell = TileAsymUnitCell(unitCell, nks);
convMat = BasesToConvMat(a1, a2, a3);
fracs = superCell(3 : 5, :) - nks' / 2;
carts = convMat * fracs;

figure;
PlotUnitCell3D(convMat, superCell);

%% matdyn.modes
[qs, bands, eigenVecs] = ReadMatdynModes('tests/qe/matdyn5x5x5.modes');
thr = 0.2;
nPh = 1;

%% sampling:
lx = norm(a1) * nks(1);
ly = norm(a2) * nks(2);
nx = 512;
ny = 512;
x = InitAxis(lx, nx);
y = InitAxis(ly, ny);
fx = InitFreqAxis(lx, nx);
fy = InitFreqAxis(ly, ny);

%% transition function:
iq = 26;
iBand = 24;
eLoss = 4.136 * abs(bands(iq, iBand));
h = PhTransitionPotential(superCell(1, :), carts, fracs, qs, bands, eigenVecs, ...
    iq, iBand, lx, ly, nx, ny, 300, nPh, thr);

%% plot transition probability
rH = ifftshift(fft2(fftshift(h)));
rHSqr = abs(rH.^2);
blurredRhSqr = imgaussfilt(rHSqr, 1);
figure;
subplot(1, 2, 1);
imagesc(x, y, blurredRhSqr);
colormap('jet');
axis square;
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title('$|H(\mathbf{r})|^2$', 'Interpreter', 'latex');

hSqr = abs(h.^2);
blurredHSqr = imgaussfilt(hSqr, 1);
subplot(1, 2, 2);
imagesc(fx, fy, blurredHSqr);
colormap('jet');
axis square;
xlabel('$f_x$ (1 / \AA)', 'Interpreter', 'latex');
ylabel('$f_y$ (1 / \AA)', 'Interpreter', 'latex');
title('$|H(\mathbf{q})|^2$', 'Interpreter', 'latex');