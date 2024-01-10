% test_ph_cross_section_plane_wave.m
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

%% cross section:
keV = 300;
wavLen0 = HighEnergyWavLen_X(keV);
wave = ones(ny, nx);
detector = CircApert_X(lx, ly, nx, ny, wavLen0, 5);
[phCs, eLosses] = PhCrossSection(wave, superCell(1, :), carts, fracs, qs, ...
    bands, eigenVecs, lx, ly, keV, detector, nPh, thr);

%% sum the results:
dmeV = 0.1;
meV = 0 : dmeV : 80;
nmeV = length(meV);
phSpectra = zeros(1, nmeV);
for imeV = 1 : nmeV
    phSpectra(imeV) = sum(phCs(eLosses > meV(imeV) - 0.5 * dmeV & ...
        eLosses < meV(imeV) + 0.5 * dmeV), 'all');
end

filtPhSpectra = gaussfilt(meV, phSpectra, 1);

figure;
plot(meV, filtPhSpectra);
xlabel('energy loss (meV)');
ylabel('cross section (a.u.)');