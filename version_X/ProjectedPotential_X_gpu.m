function projPot = ProjectedPotential_X_gpu(Lx, Ly, Nx, Ny, crysMat, deltaSq)
%ProjectedPotential.m calculates the projected potential of a series of
%atoms on a slice.
%   Lx, Ly, Nx, Ny -- sampling parameters;
%   AtomCoordType -- [T1, ..., TN; P1, ..., PN; x1, ..., xN; y1, ..., yN; z1, ..., zN];
%       where T denotes the atomic types represented by the atomic numbers
%       ranging from 1 to 103; P denotes the proportion of the elemnetal
%       concentration; third to fifth row are the atomic coordinates,
%       whether fractional or orthogonal;
%   deltaSq -- criterion to remove the singularity induced by coord 0.
% X denotes an experimental version.

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

atomNum = size(crysMat, 2);

dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);

a = 0.529; % Bohr radius in angstrom
e = 14.4; % elemental charge in volt - angstrom

scattParam = load('Scattering_Factors.txt');

projPot = zeros(size(X), "single", "gpuArray");
for i = 1 : atomNum
    atomX = crysMat(3, i);
    atomY = crysMat(4, i);
    atomType = crysMat(1, i);
    eleProp = crysMat(2, i);
    radiusSqr = gpuArray(single((X - atomX).^2 + (Y - atomY).^2));
    radiusSqr(radiusSqr < deltaSq) = deltaSq;
    startIndex = 3 * (atomType - 1) + 1;
    A = [scattParam(startIndex, 1), scattParam(startIndex, 3), scattParam(startIndex + 1, 1)];
    B = [scattParam(startIndex, 2), scattParam(startIndex, 4), scattParam(startIndex + 1, 2)];
    C = [scattParam(startIndex + 1, 3), scattParam(startIndex + 2, 1), scattParam(startIndex + 2, 3)];
    D = [scattParam(startIndex + 1, 4), scattParam(startIndex + 2, 2), scattParam(startIndex + 2, 4)];
    for j = 1:3
        projPot = projPot + eleProp * (4 * pi^2 * A(j) * besselk(0, 2 * pi * sqrt(radiusSqr) * sqrt(B(j)))...
                   + 2 * pi^2 * C(j) / D(j) * exp(-pi^2 * (radiusSqr) / D(j)));
    end
end
projPot = a * e * projPot;

end




