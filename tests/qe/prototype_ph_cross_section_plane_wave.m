% prototype_ph_cross_section.m
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
fx = InitFreqAxis(lx, nx);
fy = InitFreqAxis(ly, ny);

%% cross section:
keV = 300;
r_wave = ones(ny, nx);
[nq, nBand] = size(bands);
eLosses = 4.136 * abs(bands);
wavLen0 = HighEnergyWavLen_X(keV);
k0 = 2 * pi / wavLen0;
phSpectra = zeros(nq, nBand);

detector = CircApert_X(lx, ly, nx, ny, wavLen0, 2);

nTasks = nq * nBand;
wb = waitbar(0, 'init...');
for iq = 1 : nq
    for iBand = 1 : nBand
        k_h = PhTransitionPotential(superCell(1, :), carts, fracs, qs, bands, eigenVecs, ...
            iq, iBand, lx, ly, nx, ny, keV, nPh, thr);
        r_h = ifftshift(ifft2(fftshift(k_h)));
        r_phWave = r_wave .* r_h;
        k_phWave = ifftshift(fft2(fftshift(r_phWave)));
        k_phWaveI = abs(k_phWave.^2);

        wavLenN = HighEnergyWavLen_X(300 - 1e-6 * eLosses(iq, iBand));
        kN = 2 * pi / wavLenN;

        phSpectra(iq, iBand) = kN / k0 * sum(k_phWaveI .* detector, 'all');

        doneTasks = (iq - 1) * nBand + iBand;
        waitbar(doneTasks / nTasks, wb, [num2str(doneTasks), ' / ', num2str(nTasks)]);
    end
end
close(wb);

%% sum the results:
dmeV = 0.1;
meV = 0 : dmeV : 80;
nmeV = length(meV);
phSpectraProf = zeros(1, nmeV);
for imeV = 1 : nmeV
    phSpectraProf(imeV) = sum(phSpectra(eLosses > meV(imeV) - 0.5 * dmeV & ...
        eLosses < meV(imeV) + 0.5 * dmeV), 'all');
end

filtPhSpectraProf = gaussfilt(meV, phSpectraProf, 1);

figure;
plot(meV, filtPhSpectraProf);
xlabel('energy loss (meV)');
ylabel('cross section (a.u.)');