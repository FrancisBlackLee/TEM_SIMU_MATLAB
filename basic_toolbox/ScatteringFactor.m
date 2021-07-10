function [sf] = ScatteringFactor(atomType, q)
%Calcuate scattering factors with respect to element type and scattering 
%angle.
%   AtomType: type of the atom (atomic number);
%   q: scattering angle;

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

scattParam = load('Scattering_Factors.txt');

startIndex = 3 * (atomType - 1) + 1;
paramA = [scattParam(startIndex, 1),...
    scattParam(startIndex, 3),...
    scattParam(startIndex + 1, 1)];
paramB = [scattParam(startIndex, 2),...
    scattParam(startIndex, 4),...
    scattParam(startIndex + 1, 2)];
paramC = [scattParam(startIndex + 1, 3),...
    scattParam(startIndex + 2, 1),...
    scattParam(startIndex + 2, 3)];
paramD = [scattParam(startIndex + 1, 4),...
    scattParam(startIndex + 2, 2),...
    scattParam(startIndex + 2, 4)];

sf = 0;
for i = 1 : 3
    sf = sf + paramA(i) ./ (q.^2 + paramB(i)) + paramC(i) *...
        exp(-paramD(i) * q.^2);
end
end

