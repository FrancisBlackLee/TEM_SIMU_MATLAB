function projPot = ProjectedPotential(Lx, Ly, Nx, Ny, atomType, atomX, atomY)
%ProjectedPotential.m calculates the projected potential of a series of
%atoms on a slice.
%   Lx, Ly, Nx, Ny -- sampling parameters;
%   AtomZ -- atomic numbers of the input atom series;
%   AtomX, AtomY -- atomic coordinates corresponding to the atomic series;

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

atomNum = length(atomType);

dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
deltaSqr = (0.5 * min(dx, dy))^2;

a = 0.529; % Bohr radius in angstrom
e = 14.4; % elemental charge in volt - angstrom

scattParam = load('Scattering_Factors.txt');

projPot = zeros(size(X));
for i = 1 : atomNum
    radiusSqr = (X - atomX(i)).^2 + (Y - atomY(i)).^2;
    radiusSqr(radiusSqr < deltaSqr) = deltaSqr;
    startIndex = 3 * (atomType(i) - 1) + 1;
    A = [scattParam(startIndex, 1), scattParam(startIndex, 3), scattParam(startIndex + 1, 1)];
    B = [scattParam(startIndex, 2), scattParam(startIndex, 4), scattParam(startIndex + 1, 2)];
    C = [scattParam(startIndex + 1, 3), scattParam(startIndex + 2, 1), scattParam(startIndex + 2, 3)];
    D = [scattParam(startIndex + 1, 4), scattParam(startIndex + 2, 2), scattParam(startIndex + 2, 4)];
    for j = 1:3
        projPot = projPot + 4 * pi^2 * A(j) * besselk(0, 2 * pi * sqrt(radiusSqr) * sqrt(B(j)))...
                   + 2 * pi^2 * C(j) / D(j) * exp(-pi^2 * (radiusSqr) / D(j));
    end
end
projPot = a * e * projPot;

end

