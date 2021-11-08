function [ u2 ] = propTF_1( u1, Lx, Ly, lambda, z, bwlProp )
%propagation - transfer function approach
%   assumes same x and y side lengths and uniform sampling
%   u1 - source plane field
%   Lx, Ly -- source and observation plane side length
%   lambda -- wavelength
%   z - propagation distance
%   u2 - observation plane field

[Ny, Nx] = size(u1);          %get input field array size
dx = Lx / Nx;
dy = Ly / Ny;%sample interval
k = 2 * pi / lambda;           %wavenumber

fx = -1 / (2*dx) : 1 / Lx : 1 / (2*dx) - 1 / Lx;
fy = -1 / (2*dy) : 1 / Ly : 1 / (2*dy) - 1 / Ly;%freq coords
[FX, FY] = meshgrid(fx, fy);

H = exp(-1i * pi * lambda * z * (FX.^2 + FY.^2));      %trans func
H = fftshift(H);                             %shift trans func
U1 = fft2(fftshift(u1));                     %shift, fft src field
U2 = H .* U1;                                  %multiply
if nargin == 6
    FR = fftshift(sqrt(FX.^2 + FY.^2));
    U2 = U2 .* (FR <= bwlProp / (2 * max(dx, dy)));
end
u2 = ifftshift(ifft2(U2));                   %inv fft, center obs field

end

