function [projPot] = RadialProjectedPotential(atomType, radius)
%RadialProjectedPotential.m calculates the radial projected atomic potential.
%   atomType -- atomic type, Z;
%   radius -- polar radial coordinates;
%   projPot -- radial projected potential;

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

a = 0.529; % Bohr radius in angstrom
e = 14.4; % elemental charge in volt - angstrom
scattParam = load('Scattering_Factors.txt');
startIndex = 3 * (atomType - 1) + 1;
A = [scattParam(startIndex, 1), scattParam(startIndex, 3), scattParam(startIndex + 1, 1)];
B = [scattParam(startIndex, 2), scattParam(startIndex, 4), scattParam(startIndex + 1, 2)];
C = [scattParam(startIndex + 1, 3), scattParam(startIndex + 2, 1), scattParam(startIndex + 2, 3)];
D = [scattParam(startIndex + 1, 4), scattParam(startIndex + 2, 2), scattParam(startIndex + 2, 4)];

projPot = zeros(size(radius));
for i = 1:3
    projPot = projPot + 4 * pi^2 * A(i) * besselk(0, 2 * pi * radius * sqrt(B(i)))...
               + 2 * pi^2 * C(i) / D(i) * exp(-pi^2 * radius.^2 / D(i));
end

projPot = a * e * projPot;

end

