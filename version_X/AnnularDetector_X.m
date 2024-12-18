function [detector] = AnnularDetector_X(lowAngle, highAngle, wavLen, Lx, Ly, Nx, Ny, nBin)
%AnularDetector_X.m generates the detector for ADF-STEM mode.
%   lowAngle, highAngle -- describe the shape of the detector in mrad;
%   wavLen -- wavelength of the electron beam;
%   Lx, Ly, Nx, Ny -- sampling parameters;
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

if nargin == 7
    nBin = 1;
end

fx = InitFreqAxis(Lx, Nx, nBin);
fy = InitFreqAxis(Ly, Ny, nBin);
[FX, FY] = meshgrid(fx, fy);
freqSqr = FX.^2 + FY.^2;
detector = ((freqSqr < (highAngle * 1e-3 / wavLen)^2) &...
    (freqSqr > (lowAngle * 1e-3 / wavLen)^2));

end



