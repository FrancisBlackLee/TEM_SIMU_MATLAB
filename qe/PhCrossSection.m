function [phCs, eLosses] = PhCrossSection(wave, aTypes, aCarts, aFracs, qs, ...
    bands, eigenVecs, lx, ly, keV, detector, nPh, thr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[ny, nx] = size(wave);

[nq, nBand] = size(bands);
eLosses = 4.136 * abs(bands);
wavLen0 = HighEnergyWavLen_X(keV);
k0 = 2 * pi / wavLen0;

phCs = zeros(nq, nBand);
wb = waitbar(0, 'init...');
shiftedWave = fftshift(wave);
shiftedDetector = fftshift(detector);
for iq = 1 : nq
    parfor iBand = 1 : nBand
        k_h = PhTransitionPotential(aTypes, aCarts, aFracs, qs, bands, eigenVecs, ...
            iq, iBand, lx, ly, nx, ny, keV, nPh, thr);
        r_h = ifft2(fftshift(k_h));
        r_phWave = shiftedWave .* r_h;
        k_phWave = fft2(r_phWave);
        k_phWaveI = abs(k_phWave.^2);

        wavLenN = HighEnergyWavLen_X(keV - 1e-6 * eLosses(iq, iBand));
        kN = 2 * pi / wavLenN;

        phCs(iq, iBand) = kN / k0 * sum(k_phWaveI .* shiftedDetector, 'all');
    end
    waitbar(iq / nq, wb, [num2str(iq), ' / ', num2str(nq)]);
end
close(wb);

end