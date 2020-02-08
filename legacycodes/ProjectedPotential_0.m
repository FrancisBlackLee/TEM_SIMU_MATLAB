%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2020  Francis Black Lee and Li Xian

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
function Proj_Pot = ProjectedPotential_0(Lx, Ly, Nx, Ny, AtomCoordType)
%ProjectedPotential.m calculates the projected potential of a series of
%atoms on a slice.
%   Lx, Ly, Nx, Ny -- sampling parameters;
%   AtomCoordType -- [x1, ..., xN; y1, ..., yN; z1, ..., zN; Type1, ..., TypeN];

AtomNum = size(AtomCoordType, 2);

dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
% delta = 0.1 * dx;
deltaSq = 0.5 * (dx * dx + dy * dy);

a = 0.529; % Bohr radius in angstrom
e = 14.4; % elemental charge in volt - angstrom

FileName = mfilename('fullpath');
FileName = strcat(FileName, '.m');
[filepath, name, ext] = fileparts(FileName);
Pot_txt_name = fullfile(filepath, 'Scattering_Factors.txt');
Scatt_Fac = load(Pot_txt_name);

Proj_Pot = zeros(size(X));
for i = 1 : AtomNum
    AtomX = AtomCoordType(1, i);
    AtomY = AtomCoordType(2, i);
    AtomType = AtomCoordType(4, i);
    RHOsq = (X - AtomX).^2 + (Y - AtomY).^2;
    RHOsq(RHOsq < deltaSq) = deltaSq;
    StartIndex = 3 * (AtomType - 1) + 1;
    A = [Scatt_Fac(StartIndex, 1), Scatt_Fac(StartIndex, 3), Scatt_Fac(StartIndex + 1, 1)];
    B = [Scatt_Fac(StartIndex, 2), Scatt_Fac(StartIndex, 4), Scatt_Fac(StartIndex + 1, 2)];
    C = [Scatt_Fac(StartIndex + 1, 3), Scatt_Fac(StartIndex + 2, 1), Scatt_Fac(StartIndex + 2, 3)];
    D = [Scatt_Fac(StartIndex + 1, 4), Scatt_Fac(StartIndex + 2, 2), Scatt_Fac(StartIndex + 2, 4)];
    for j = 1:3
        Proj_Pot = Proj_Pot + 4 * pi^2 * A(j) * besselk(0, 2 * pi * sqrt(RHOsq) * sqrt(B(j)))...
                   + 2 * pi^2 * C(j) / D(j) * exp(-pi^2 * (RHOsq) / D(j));
    end
end
Proj_Pot = a * e * Proj_Pot;

end

