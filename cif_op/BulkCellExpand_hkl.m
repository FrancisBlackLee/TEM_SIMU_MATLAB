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
function [atomCoordMat] = BulkCellExpand_hkl(atomSiteMat, cellLengths,...
    cellAngles, hkl, width, thickness)
%CellExpand_hkl() rotates the unit cell to the given orientation, duplicate
%the unit cell and rehapes the cell.
% Input:
%   atomSiteMat -- atomic site matrix;
%   cellLengths -- element 1 for cell length a, 2 for cell length b and 3
%       for cell length c;
%   cellAngles -- element 1 for cell angle alpha (between bases a and c);
%       2 for cell angle beta (between bases b and c) and
%       3 for cell angle gamma (between bases a and b);
%   hkl -- Miller indices;
%   sideLength -- the size of the square sample you want;
% Output:
%   atomCoordMat -- atomic coordinate matrix (cartesian);

radius = sqrt(2 * width^2 + thickness^2) / 2;
atomCoordMat = CreateNanoCluster_hkl(atomSiteMat, cellLengths, cellAngles, hkl, radius);

tolerance = 1e-8;
atomCoordMat(:, (atomCoordMat(5, :) > thickness / 2 + tolerance) |...
    (atomCoordMat(5, :) < -thickness / 2 - tolerance)) = [];

end

