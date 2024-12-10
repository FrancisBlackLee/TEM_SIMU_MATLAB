function [dwfs] = MatdynModesToDebyeWallerFactor_X(qs, bands, eigenVecs, mass, convMat, T, thr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 6
    thr = 1e-1;
end

refDwfs = MatdynModesToDebyeWallerFactor(qs, bands, eigenVecs, mass, T, thr);

nAtom = length(mass);
dwfs = zeros(3, 3, nAtom);

hbar = 1.054571817e-34;
m = mass * 1.66053906660e-27;
kb = 1.380649e-23;

nq = size(qs, 1);
nBand = 3 * nAtom;

% invConvMat = inv(convMat);

% lattConsts = [norm(convMat(:, 1)), norm(convMat(:, 2)), norm(convMat(:, 3))];

for iAtom = 1 : nAtom
    for alpha = 1 : 3
        for beta = 1 : alpha
            for iq = 1 : nq
                for iBand = 1 : nBand
                    tmpVec = eigenVecs(iAtom, :, iBand, iq)';
                    tmpVec = convMat \ tmpVec;
                    omega = 2 * pi * bands(iq, iBand) * 1e12;
                    if abs(bands(iq, iBand)) > thr
                        dwfs(alpha, beta, iAtom) = dwfs(alpha, beta, iAtom) +...
                            coth(hbar * omega / (2 * kb * T)) *...
                            real(tmpVec(alpha) *...
                            conj(tmpVec(beta))) / omega / m(iAtom);
                    end
                end
            end

            dwfs(beta, alpha, iAtom) = dwfs(alpha, beta, iAtom);
        end
    end
end

dwfs = 8 * pi^2 * dwfs * hbar / (2 * nq) * 1e20;
for iAtom = 1 : nAtom
    scalars = zeros(1, 3);
    for i = 1 : 3
        scalars(i) = sqrt(refDwfs(i, i, iAtom) / dwfs(i, i, iAtom));
    end

    for alpha = 1 : 3
        for beta = 1 : 3
            dwfs(alpha, beta, iAtom) = scalars(alpha) * scalars(beta) * dwfs(alpha, beta, iAtom);
        end
    end
end

end