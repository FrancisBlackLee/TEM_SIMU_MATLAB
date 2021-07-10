function [projPot] = MonoProjPot_imtl_X(atomType, eleProp, xyCoords,...
    Lx, Ly, Nx, Ny)
%MonoProjPot_imtl_X.m calculates the projected potential for one type of
%atom based on functions ProjectedPotential and imtranslate.
%   atomType -- atomic type;
%   eleProp -- elemental proportion;
%   xyCoords -- x and y Cartesian coordinates, format: [x1, ..., xN; y1,
%       ..., yN];
%   Lx, Ly -- sampling side lengths;
%   Nx, Ny -- sampling number;
% Note: X denotes an experimental version!

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

atomNum = size(xyCoords, 2);
singleProjPot = ProjectedPotential(Lx, Ly, Nx, Ny, atomType, 0, 0);
% % move the center of the specimen to the origin:
% tmpXyCoords = xyCoords - mean(xyCoords, 2);
% iterate over all atoms:
dx = Lx / Nx;
dy = Ly / Ny;
unitTransVec = [dx; dy];
projPot = zeros(Ny, Nx);
for atomIdx = 1 : atomNum
    transVec = xyCoords(:, atomIdx) ./ unitTransVec;
    tmpProjPot = imtranslate(singleProjPot, transVec', 'nearest');
    projPot = projPot + eleProp(atomIdx) * tmpProjPot;
end

end

