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
function [sf] = ScatteringFactor(AtomType, q)
%Calcuate scattering factors with respect to element type and scattering 
%angle.
%   AtomType: type of the atom (atomic number);
%   q: scattering angle;

Scatt_Fac = load('Scattering_Factors.txt');

StartIndex = 3 * (AtomType - 1) + 1;
A = [Scatt_Fac(StartIndex, 1), Scatt_Fac(StartIndex, 3), Scatt_Fac(StartIndex + 1, 1)];
B = [Scatt_Fac(StartIndex, 2), Scatt_Fac(StartIndex, 4), Scatt_Fac(StartIndex + 1, 2)];
C = [Scatt_Fac(StartIndex + 1, 3), Scatt_Fac(StartIndex + 2, 1), Scatt_Fac(StartIndex + 2, 3)];
D = [Scatt_Fac(StartIndex + 1, 4), Scatt_Fac(StartIndex + 2, 2), Scatt_Fac(StartIndex + 2, 4)];

sf = 0;
for i = 1 : 3
    sf = sf + A(i)./(q.^2 + B(i)) + C(i)*exp(-D(i)*q.^2);
end
end

