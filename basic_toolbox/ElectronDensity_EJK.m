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
function [rho] = ElectronDensity_EJK(atomType, r)
%ElectronDensity_EJK.m calculates the electron density using the
%parameterization given by E. J. Kirkland.
%   atomType -- atomic type, Z;
%   r -- radial coordinates;
%   rho -- electron density;

a = 0.529; % Bohr radius in angstrom
scattFac = load('Scattering_Factors.txt');
startIndex = 3 * (atomType - 1) + 1;
A = [scattFac(startIndex, 1), scattFac(startIndex, 3), scattFac(startIndex + 1, 1)];
B = [scattFac(startIndex, 2), scattFac(startIndex, 4), scattFac(startIndex + 1, 2)];
C = [scattFac(startIndex + 1, 3), scattFac(startIndex + 2, 1), scattFac(startIndex + 2, 3)];
D = [scattFac(startIndex + 1, 4), scattFac(startIndex + 2, 2), scattFac(startIndex + 2, 4)];

rho = zeros(size(r));
% note: the delta function part is excluded:
for i = 1 : 3
    rho = rho + 2 * pi^3 * a * A(i) * B(i) * exp(-2 * pi * sqrt(B(i)) * r) ./ r +...
        pi^5.5 * a * C(i) / (D(i))^3.5 * (2 * r.^2 - 3 * D(i) / pi^2) .*...
        exp(-pi^2 * r.^2 / D(i));
end

end

