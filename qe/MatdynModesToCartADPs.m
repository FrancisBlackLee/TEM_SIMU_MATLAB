function [cartAdps] = MatdynModesToCartADPs(qs, bands, eigenVecs, mass, T, thr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 5
    thr = 1e-1;
end

nAtom = length(mass);
cartAdps = zeros(3, 3, nAtom);

hbar = 1.054571817e-34;
m = mass * 1.66053906660e-27;
kb = 1.380649e-23;

nq = size(qs, 1);
nBand = 3 * nAtom;

for iAtom = 1 : nAtom
    for alpha = 1 : 3
        for beta = 1 : alpha
            for iq = 1 : nq
                for iBand = 1 : nBand
                    omega = 2 * pi * bands(iq, iBand) * 1e12;
                    if abs(bands(iq, iBand)) > thr
                        cartAdps(alpha, beta, iAtom) = cartAdps(alpha, beta, iAtom) +...
                            coth(hbar * omega / (2 * kb * T)) *...
                            eigenVecs(iAtom, alpha, iBand, iq) *...
                            conj(eigenVecs(iAtom, beta, iBand, iq)) /omega / m(iAtom);
                    end
                end
            end

            cartAdps(beta, alpha, iAtom) = cartAdps(alpha, beta, iAtom);
        end
    end
end

cartAdps = real(cartAdps * hbar / (2 * nq) * 1e20);

end