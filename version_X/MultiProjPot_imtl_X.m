function [projPot] = MultiProjPot_imtl_X(atomTypeCoords, Lx, Ly, Nx, Ny)
%MultiProjPot_imtl_X.m calculates the projected potential of multiple types
%of atoms, this function iteratively calls MonoProjPot_imtl_X.m to
%calculate the component for each atomic type.
%   atomTypeCoord -- matrix that contains the specimen information, where
%       the first row denotes the atomic types, the second row denotes the
%       elemental proportion, and the third to the fifth row denotes the
%       atomic cartesian coordinates. Syntax: [type; prop; X; Y; Z];
%       Since all the atoms included in this matrix are all on one slice,
%       the fifth row is not required strictly;
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

% sort atomTypeCoord in an ascending order by atomic type:
[sortedType, typeOrder] = sort(atomTypeCoords(1, : ));
atomTypeCoords = atomTypeCoords( : , typeOrder);
projPot = zeros(Ny, Nx);
if length(sortedType) > 1
    lastType = 1;
    for atomIdx = 1 : length(sortedType)
        if sortedType(atomIdx) ~= sortedType(lastType)
            projPot = projPot + MonoProjPot_imtl_X(sortedType(lastType),...
                atomTypeCoords(2, lastType : atomIdx - 1),...
                atomTypeCoords(3 : 4, lastType : atomIdx - 1),...
                Lx, Ly, Nx, Ny);
            lastType = atomIdx;
        end
    end
    projPot = projPot + MonoProjPot_imtl_X(sortedType(lastType),...
        atomTypeCoords(2, lastType : atomIdx),...
        atomTypeCoords(3 : 4, lastType : atomIdx),...
        Lx, Ly, Nx, Ny);
else
    projPot = MonoProjPot_imtl_X(sortedType, atomTypeCoords(2, :),...
        atomTypeCoords(3 : 4, :), Lx, Ly, Nx, Ny);
end

end

