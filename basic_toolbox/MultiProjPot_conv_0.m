function [projPot] = MultiProjPot_conv_0(fracTypeCoord, expanNum, lattConst,...
    Lx, Ly, Nx, Ny)
%MultiProjPot_conv_0.m calculates the projected potential of one slice
%using convolution if this slice contains more than one type of atoms.
%   fracTypeCoord -- matrix that contains the lattice information, where
%       the first row denotes the atomic types, the second to the forth row
%       denote the fractional atomic coordinates, since the involving atoms
%       are all on the same slice, the forth row is not required;
%   expanNum -- expansion of the unit cell, syntax: [numX, numY];
%   lattConst -- lattice constants, syntax: [a, b];
%   Lx, Ly, Nx, Ny -- sampling parameters;

%   P.S.: I have come up with another way to further speed this process
%       when massive atoms are to be calculated, anyway, there are still
%       some bugs. When the bugs are fixed, the released version will be
%       named as MonoProjPot_conv_1.m and provide the same interface as
%       MonoProjPot_conv_0.m does, so it can be called in this function by
%       directly change the version number.

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

% sort ScaleTypeCoord in an ascending order by atomic type:
[sortedType, typeOrder] = sort(fracTypeCoord(1, : ));
fracTypeCoord = fracTypeCoord( : , typeOrder);
projPot = 0;
if length(sortedType) > 1
    lastType = 1;
    for atomIdx = 1 : length(sortedType)
        if sortedType(atomIdx)~=sortedType(lastType)
            projPot = projPot + MonoProjPot_conv_0(sortedType(lastType), ...
                fracTypeCoord(2 : 3 , lastType : atomIdx-1), expanNum,...
                lattConst, Lx, Ly, Nx, Ny);
            lastType = atomIdx;
        end
    end
    projPot = projPot + MonoProjPot_conv_0(sortedType(lastType), ...
            fracTypeCoord(2 : 3 , lastType : atomIdx), expanNum, lattConst,...
            Lx, Ly, Nx, Ny);
else
    projPot = MonoProjPot_conv_0(sortedType, fracTypeCoord(2 : 3, : ),...
        expanNum, lattConst, Lx, Ly, Nx, Ny);
end

end

