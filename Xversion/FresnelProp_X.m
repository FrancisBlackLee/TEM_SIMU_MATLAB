%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019  Francis Black Lee and Li Xian

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outllok.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

