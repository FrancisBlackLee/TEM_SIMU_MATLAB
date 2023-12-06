% prototype_ph_transition_potential.m
clc;
clear;
close all;
%% load pwscfin
pwscf = ReadPwscfInput('tests/qe/si_cubic.1_scf.in');
[unitCell, a1, a2, a3, b1, b2, b3] = PwscfInputToFracCoords(pwscf);
pwscf = ExpandPwscfAtomMass(pwscf);
nks = [8, 8, 8];
nCell = nks(1) * nks(2) * nks(3);
superCell = TileAsymUnitCell(unitCell, nks);
convMat = BasesToConvMat(a1, a2, a3);
carts = convMat * superCell(3 : 5, :);
fracs = superCell(3 : 5, :);

figure;
PlotUnitCell3D(convMat, superCell);

%% matdyn.modes
[qs, bands, eigenVecs] = ReadMatdynModes('tests/qe/matdyn.modes');
thr = 0.2;
nj = 1;

%% sampling:
lx = norm(a1) * nks(1);
ly = norm(a2) * nks(2);
nx = 512;
ny = 512;
fx = InitFreqAxis(lx, nx);
fy = InitFreqAxis(ly, ny);
[fxMesh, fyMesh] = meshgrid(fx, fy);
frMesh = sqrt(fxMesh.^2 + fyMesh.^2);

%% transition function:
iq = 1;
q = qs(iq, :);
nAtom = size(superCell, 2);
nBand = size(bands, 2);
h = zeros(ny, nx);
hbar = 1.054571817e-34;
nUnitCellAtom = size(unitCell, 2);
for iAtom = 1 : nAtom
    h1 = exp(-2 * pi * 1i * (fxMesh * carts(1, iAtom) + ...
        fyMesh * carts(2, iAtom))) .* ...
    ScatteringFactor(superCell(1, iAtom), frMesh);

    h2 = ones(ny, nx);
    iUnitCellAtom = mod(iAtom - 1, nUnitCellAtom) + 1;
    for iBand = 1 : nBand
        if abs(bands(iq, iBand)) > thr
            omega = 2 * pi * bands(iq, iBand) * 1e12;
            dwf = pi^2 * hbar / omega;
            epsilon = 1 / sqrt(nCell) * real(eigenVecs(iUnitCellAtom, :, iBand, iq) * ...
                exp(2 * pi * 1i * dot(q, fracs(:, iAtom)')));
            qEpsilon = fxMesh * epsilon(1) + fyMesh * epsilon(2);
            h2 = h2 * (-1i * sqrt(2 * dwf) * qEpsilon).^nj / factorial(nj) .* ...
                exp(-dwf * qEpsilon.^2);
        end
    end
    h = h + h1 .* h2;
end

%% plot transition probability
hSqr = abs(h.^2);
figure;
imagesc(fx, fy, hSqr);
