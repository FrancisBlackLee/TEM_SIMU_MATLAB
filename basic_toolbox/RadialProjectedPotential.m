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
function [projPot] = RadialProjectedPotential(atomType, rho)
%RadialProjectedPotential.m calculates the radial projected atomic potential.
%   atomType -- atomic type, Z;
%   rho -- radial coordinates;
%   projPot -- radial projected potential;

a = 0.529; % Bohr radius in angstrom
e = 14.4; % elemental charge in volt - angstrom
scattFac = load('Scattering_Factors.txt');
startIndex = 3 * (atomType - 1) + 1;
A = [scattFac(startIndex, 1), scattFac(startIndex, 3), scattFac(startIndex + 1, 1)];
B = [scattFac(startIndex, 2), scattFac(startIndex, 4), scattFac(startIndex + 1, 2)];
C = [scattFac(startIndex + 1, 3), scattFac(startIndex + 2, 1), scattFac(startIndex + 2, 3)];
D = [scattFac(startIndex + 1, 4), scattFac(startIndex + 2, 2), scattFac(startIndex + 2, 4)];

projPot = zeros(size(rho));
for i = 1:3
    projPot = projPot + 4 * pi^2 * A(i) * besselk(0, 2 * pi * rho * sqrt(B(i)))...
               + 2 * pi^2 * C(i) / D(i) * exp(-pi^2 * rho.^2 / D(i));
end

projPot = a * e * projPot;

end

