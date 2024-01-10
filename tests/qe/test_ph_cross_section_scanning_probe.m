% test_ph_cross_section_scanning_probe.m
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
carts = convMat * superCell(3 : 5, :);
fracs = superCell(3 : 5, :);

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

%% scanning and cross section:
params = InitProbeParams_X();
params.KeV = 300;
wavLen = HighEnergyWavLen_X(params.KeV);
params.aperture = CircApert_X(lx, ly, nx, ny, wavLen, 15);
nProbe = 5;
pxs = linspace(0, norm(a1) / 2, nProbe);
params.yp = 0;
params.Lx = lx;
params.Ly = ly;
params.Nx = nx;
params.Ny = ny;
detector = CircApert_X(lx, ly, nx, ny, wavLen, 5);
figure;
for iProbe = 1 : nProbe
    params.xp = pxs(iProbe);
    probe = GenerateProbe_X(params);
    [phCs, eLosses] = PhCrossSection(probe, superCell(1, :), carts, fracs, qs, ...
        bands, eigenVecs, lx, ly, params.KeV, detector, nPh, thr);
    filename = ['eloss_data_probe_', num2str(iProbe), '.mat'];
    save(fullfile('tests/qe', filename), 'eLosses', 'phCs');
end