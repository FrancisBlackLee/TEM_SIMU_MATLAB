function [detectors] = Dpc4Detectors_X(lowAngle, highAngle, rotAngle,...
    wavLen, Lx, Ly, Nx, Ny, nBin)
%Dpc4Detectors_X generates 4-piece segmented DPC detectors
%   lowAngle, highAngle -- describe the shape of the detector in mrad;
%   rotAngle -- rotation angle of the detectors;
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

if nargin == 8
    nBin = 1;
end

fx = InitFreqAxis(Lx, Nx, nBin);
fy = InitFreqAxis(Ly, Ny, nBin);
[FX, FY] = meshgrid(fx, fy);
freqSqr = FX.^2 + FY.^2;

lowFreqSqr = (lowAngle * 1e-3 / wavLen)^2;
highFreqSqr = (highAngle * 1e-3 / wavLen)^2;

mask = (freqSqr > lowFreqSqr) & (freqSqr < highFreqSqr);

detectors = zeros(round(Ny / nBin), round(Nx / nBin), 4);
% quadrant 1:
detectors(:, :, 1) = (mask & (FX > 0) & (FY > 0));
detectors(:, :, 1) = imrotate(detectors(:, :, 1), rotAngle, 'crop');
% quadrant 2:
detectors(:, :, 2) = (mask & (FX < 0) & (FY > 0));
detectors(:, :, 2) = imrotate(detectors(:, :, 2), rotAngle, 'crop');
% quadrant 3:
detectors(:, :, 3) = (mask & (FX < 0) & (FY < 0));
detectors(:, :, 3) = imrotate(detectors(:, :, 3), rotAngle, 'crop');
% quadrant 4:
detectors(:, :, 4) = (mask & (FX > 0) & (FY < 0));
detectors(:, :, 4) = imrotate(detectors(:, :, 4), rotAngle, 'crop');

end


