function [electronDensity] = ElectronDensity_EJK(atomType, r)
%ElectronDensity_EJK.m calculates the electron density using the
%parameterization given by E. J. Kirkland.
%   atomType -- atomic type, Z;
%   r -- radial coordinates;
%   rho -- electron density;

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
scattParam = load('Scattering_Factors.txt');
startIndex = 3 * (atomType - 1) + 1;
A = [scattParam(startIndex, 1), scattParam(startIndex, 3), scattParam(startIndex + 1, 1)];
B = [scattParam(startIndex, 2), scattParam(startIndex, 4), scattParam(startIndex + 1, 2)];
C = [scattParam(startIndex + 1, 3), scattParam(startIndex + 2, 1), scattParam(startIndex + 2, 3)];
D = [scattParam(startIndex + 1, 4), scattParam(startIndex + 2, 2), scattParam(startIndex + 2, 4)];

electronDensity = zeros(size(r));
% note: the delta function part is excluded:
for i = 1 : 3
    electronDensity = electronDensity + 2 * pi^3 * a * A(i) * B(i) * exp(-2 * pi * sqrt(B(i)) * r) ./ r +...
        pi^5.5 * a * C(i) / (D(i))^3.5 * (2 * r.^2 - 3 * D(i) / pi^2) .*...
        exp(-pi^2 * r.^2 / D(i));
end

end

