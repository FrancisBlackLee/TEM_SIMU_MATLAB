function probe = MeshedProbe(aberrs, wavLen, aperture, xp, yp, ...
    fxMesh, fyMesh)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

otfPhase = MeshedAberrPhase(aberrs, wavLen, fxMesh, fyMesh);
otf = aperture .* exp(-1i * otfPhase);
probe = ifftshift(ifft2(fftshift(otf .* exp(-1i * 2 * pi * (fxMesh * xp + fyMesh * yp)))));
normCoeff = sqrt(sum(abs(probe.^2), 'all'));
probe = probe / normCoeff;

end