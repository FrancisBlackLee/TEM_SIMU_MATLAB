function [projPot] = MonoProjPot_conv_0(atomType, fracCoord, expanNum,...
    lattConst, Lx, Ly, Nx, Ny, method)
%MonoProjPot_conv_0.m calculates the projected potential for one type of
%atom, more comprehensive function to be released soon.
%   fracCoord -- fractional xy coordinates, syntax: [fracX1,..., fracXN;
%       fracY1, fracYN];
%   expanNum -- number of unit cells to be expanned, sytax: [numX,numY];
%   lattConst -- lattice constants, syntax: [a, b];
%   Lx, Ly, Nx, Ny -- sampling parameters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2021  Francis Black Lee and Li Xian

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
    method = 'sf';
end

atomNum = size(fracCoord, 2);
% Sampling:
dx = Lx / Nx;
dy = Ly / Ny;
fx = InitFreqAxis(Lx, Nx);
fy = InitFreqAxis(Ly, Ny);
[FX, FY] = meshgrid(fx, fy);

if strcmp(method, 'rs')
    projPotFFT = fft2(fftshift(ProjectedPotential(Lx, Ly, Nx, Ny, atomType, 0, 0)));
elseif strcmp(method, 'sf')
    FR = sqrt(FX.^2 + FY.^2);

    a = 0.529; % Bohr radius in angstrom
    e = 14.4; % elemental charge in volt - angstrom
    scaleCoeff = 2 * pi * e * a;

    projPotFFT = fftshift(scaleCoeff * ScatteringFactor(atomType, FR));
else
    errID = 'myComponent:inputError';
    msgtext = 'Invaid parameter of method!';
    ME = MException(errID, msgtext);
    throw(ME);
end

% Build the convolution kernel:
kernel = 0;
scaleXshift = expanNum(1) / 2;
scaleYshift = expanNum(2) / 2;
for atomIdx = 1 : atomNum
    for xEpanIdx = 0 : expanNum(1) - 1
        coordX = lattConst(1) * (fracCoord(1, atomIdx) + xEpanIdx - scaleXshift);
        for yExpanIdx = 0 : expanNum(2) - 1
            coordY = lattConst(2) * (fracCoord(2, atomIdx) + yExpanIdx - scaleYshift);
            kernel = kernel + exp(-1i * 2 * pi *(FX * coordX + FY * coordY));
        end
    end
end
kernel = fftshift(kernel);
projPot = real(ifftshift(ifft2(kernel .* projPotFFT)));

if strcmp(method, 'sf')
    projPot = projPot / (dx * dy);
end

end

