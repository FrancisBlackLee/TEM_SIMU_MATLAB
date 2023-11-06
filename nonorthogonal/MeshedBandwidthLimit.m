function limited = MeshedBandwidthLimit(unlimited, fxMesh, fyMesh, bwl)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

frMesh = sqrt(fxMesh.^2 + fyMesh.^2);
tmp = fftshift(frMesh <= bwl) .* fft2(fftshift(unlimited));
limited = ifftshift(ifft2(tmp));

end