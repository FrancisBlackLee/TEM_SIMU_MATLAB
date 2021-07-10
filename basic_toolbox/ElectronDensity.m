function [electronDensity] = ElectronDensity(atomType, rCoords)
%ElectronDensity.m calculates the radial atomic electron density.
%   atomType -- atomic type, Z;
%   rCoords -- radial coordinates;
%   electronDensity -- electron density;

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

a0 = 0.529; % Bohr radius in angstrom
scattParam = load('VanDyckScattParam.txt');
paramA = scattParam(2 * atomType - 1, :);
paramB = scattParam(2 * atomType, :);

electronDensity = zeros(size(rCoords));
for i = 1 : 5
    electronDensity = electronDensity + paramA(i) / paramB(i)^2.5 *...
        exp(-2 * pi * rCoords / sqrt(paramB(i)));
end
electronDensity = 2 * pi^4 * a0 * electronDensity;

end

