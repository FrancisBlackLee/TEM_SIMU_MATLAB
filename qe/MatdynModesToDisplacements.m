function [displacements] = MatdynModesToDisplacements(superCell, a1, a2, a3,...
    b1, b2, b3, qs, bands, eigenVecs, mass, T, nConfig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

thr = 1e-4;

superCellAtomNum = size(superCell, 2);
unitCellAtomNum = length(mass);

displacements = zeros(3, superCellAtomNum, nConfig);

hbar = 1.054571817e-34;
m = mass * 1.66053906660e-27;
kb = 1.380649e-23;

nq = size(qs, 1);
nBand = 3 * unitCellAtomNum;

a1 = reshape(a1, 1, []);
a2 = reshape(a2, 1, []);
a3 = reshape(a3, 1, []);

b1 = reshape(b1, 1, []);
b2 = reshape(b2, 1, []);
b3 = reshape(b3, 1, []);

convMat = [a1; a2; a3]';

xyz = convMat * superCell(3 : 5, :);

wb = waitbar(0, 'wait...');

for iConfig = 1 : nConfig
    for iAtom = 1 : superCellAtomNum
        u = zeros(1, 3);
        for iq = 1 : nq
            if qs(iq, 1) >= 0 % only count positive-x space
                for iBand = 1 : nBand
                    q = qs(iq, :) * [b1; b2; b3];
                    omega = 2 * pi * bands(iq, iBand) * 1e12;
                    if abs(bands(iq, iBand)) > thr
                        iUnitCellAtom = mod(iAtom - 1, unitCellAtomNum) + 1;
                        sigma = sqrt(hbar / (2 * omega * m(iUnitCellAtom) *...
                            tanh(hbar * omega / (2 * kb * T))));
                        rnd1 = normrnd(0, sigma);
                        rnd2 = normrnd(0, sigma);
                        u = u + sqrt(2) * real(eigenVecs(iUnitCellAtom, :, iBand, iq) *...
                            exp(2 * pi * 1i * dot(q, xyz(:, iAtom)'))) * rnd1 -...
                            sqrt(2) * imag(eigenVecs(iUnitCellAtom, :, iBand, iq) *...
                            exp(2 * pi * 1i * dot(q, xyz(:, iAtom)'))) * rnd2;
                    end
                end
            end
        end
        displacements(:, iAtom, iConfig) = u';
    end

    wbMsg = [num2str(iConfig), '/', num2str(nConfig), ' done'];
    waitbar(iConfig / nConfig, wb, wbMsg);
end

close(wb);

nCell = superCellAtomNum / unitCellAtomNum;
displacements = displacements / sqrt(nCell) * 1e10; % angstrom

end