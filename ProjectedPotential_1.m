%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019  Francis Black Lee and Li Xian

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outllok.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Proj_Pot = ProjectedPotential_1(Lx, Ly, Nx, Ny, CrysMat)
%ProjectedPotential.m calculates the projected potential of a series of
%atoms on a slice.
%   Lx, Ly, Nx, Ny -- sampling parameters;
%   AtomCoordType -- [T1, ..., TN; P1, ..., PN; x1, ..., xN; y1, ..., yN; z1, ..., zN];
%       where T denotes the atomic types represented by the atomic numbers
%       ranging from 1 to 103; P denotes the proportion of the elemnetal
%       concentration; third to fifth row are the atomic coordinates,
%       whether fractional or orthogonal;

AtomNum = size(CrysMat, 2);

dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
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
    AtomX = CrysMat(3, i);
    AtomY = CrysMat(4, i);
    AtomType = CrysMat(1, i);
    EleProp = CrysMat(2, i);
    RHOsq = (X - AtomX).^2 + (Y - AtomY).^2;
    RHOsq(RHOsq < deltaSq) = deltaSq;
    StartIndex = 3 * (AtomType - 1) + 1;
    A = [Scatt_Fac(StartIndex, 1), Scatt_Fac(StartIndex, 3), Scatt_Fac(StartIndex + 1, 1)];
    B = [Scatt_Fac(StartIndex, 2), Scatt_Fac(StartIndex, 4), Scatt_Fac(StartIndex + 1, 2)];
    C = [Scatt_Fac(StartIndex + 1, 3), Scatt_Fac(StartIndex + 2, 1), Scatt_Fac(StartIndex + 2, 3)];
    D = [Scatt_Fac(StartIndex + 1, 4), Scatt_Fac(StartIndex + 2, 2), Scatt_Fac(StartIndex + 2, 4)];
    for j = 1:3
        Proj_Pot = Proj_Pot + EleProp * (4 * pi^2 * A(j) * besselk(0, 2 * pi * sqrt(RHOsq) * sqrt(B(j)))...
                   + 2 * pi^2 * C(j) / D(j) * exp(-pi^2 * (RHOsq) / D(j)));
    end
end
Proj_Pot = a * e * Proj_Pot;

end

