function [shifted_u2] = FresnelProp_X(shifted_u1, shifted_PropKernel)
%FresnelProp_X.m computes Fresnel propagation with input wave and prop
%kernel.
%   shifted_u1 -- fftshifted incident wave;
%   shifted_PropKernel -- fftshifted propagation kernel in Fourier space,
%       corresponding to the wave propagation theory, not restricted to
%       Fresnel propagation.
%   shifted_u2 -- fftshifted exit wave.
% Note: X denotes an experimental version! All the functions are kept
% fftshifted deliberately to reduce unnecessary computation, for it is not
% an illustration version!

shifted_u2 = ifft2(shifted_PropKernel .* fft2(shifted_u1));

end

