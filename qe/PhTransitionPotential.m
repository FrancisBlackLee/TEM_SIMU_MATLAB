function [h] = PhTransitionPotential(aTypes, aCarts, aFracs, qs, bands, ...
    eigenVecs, iq, iBand, lx, ly, nx, ny, keV, nPh, thr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

h = zeros(ny, nx);
if abs(bands(iq, iBand)) > thr
    % wavLen = HighEnergyWavLen_X(keV);
    % k0 = 2 * pi / wavLen;
    % mesh:
    fx = InitFreqAxis(lx, nx);
    fy = InitFreqAxis(ly, ny);
    [fxMesh, fyMesh] = meshgrid(fx, fy);
    frMesh = sqrt(fxMesh.^2 + fyMesh.^2);
    % fzMesh = sqrt(k0^2 - fxMesh.^2 - fyMesh.^2);
    
    % transition potential
    q = qs(iq, :);
    nAtom = length(aTypes);
    nBand = size(bands, 2);
    nUnitCellAtom = nBand / 3;
    nCell = nAtom / nUnitCellAtom;
    
    hbar = 1.054571817e-34;
    
    omega = 2 * pi * bands(iq, iBand) * 1e12;
    dwf = pi^2 * hbar / omega;
    
    for iAtom = 1 : nAtom
        h1 = exp(-2 * pi * 1i * (fxMesh * aCarts(1, iAtom) + ...
            fyMesh * aCarts(2, iAtom))) .* ...
            ScatteringFactor(aTypes(iAtom), frMesh);
    
        h2 = ones(ny, nx);
        iUnitCellAtom = mod(iAtom - 1, nUnitCellAtom) + 1;
        epsilon = sqrt(2 / nCell) * real(eigenVecs(iUnitCellAtom, :, iBand, iq) * ...
            exp(2 * pi * 1i * dot(q, aFracs(:, iAtom)')));
        qEpsilon = fxMesh * epsilon(1) + fyMesh * epsilon(2);
        for iPh = 1 : nPh
            h2 = h2 .* (-1i * sqrt(2 * dwf) * qEpsilon).^nPh / factorial(nPh) .* ...
                exp(-dwf * qEpsilon.^2);
        end
        h = h + h1 .* h2;
    end
end

end