function [u2] = propTF_2(u1, kxMesh, kyMesh, wavLen, z)
%propagation - transfer function approach
%   assumes same x and y side lengths and uniform sampling
%   u1 - source plane field
%   kxMesh, kyMesh -- k-space meshes
%   lambda -- wavelength
%   z - propagation distance
%   u2 - observation plane field

H = exp(-1i * pi * wavLen * z * (kxMesh.^2 + kyMesh.^2));      %trans func
H = fftshift(H);                             %shift trans func
U1 = fft2(fftshift(u1));                     %shift, fft src field
U2 = H .* U1;                                  %multiply
u2 = ifftshift(ifft2(U2));                   %inv fft, center obs field

end