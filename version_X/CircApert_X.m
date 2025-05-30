function [aperture] = CircApert_X(Lx, Ly, Nx, Ny, wavLen, numApert, pol, azi)
%CircApert_X.m generates a circular aperture in reciprocal space.
%   Lx, Ly, Nx, Ny -- sampling parameters, L denotes side length and N the
%       sampling number in real space;
%   wavLen -- wavelength;
%   numApert -- numerical aperture in mrad;
% Note: X denotes an experimental version!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2023  Francis Black Lee (Li Xian)

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fx = InitFreqAxis(Lx, Nx);
fy = InitFreqAxis(Ly, Ny);
[FX, FY] = meshgrid(fx, fy);
if nargin == 6
    freqSqr = FX.^2 + FY.^2;
    aperture = (freqSqr < (numApert * 1.0e-3 / wavLen)^2);
else
    tiltFreq = pol * 1.0e-3 / wavLen;
    tiltFreqX = tiltFreq * cosd(azi);
    tiltFreqY = tiltFreq * sind(azi);
    freqSqr = (FX - tiltFreqX).^2 + (FY - tiltFreqY).^2;
    aperture = (freqSqr < (numApert * 1.0e-3 / wavLen)^2);
end

end



