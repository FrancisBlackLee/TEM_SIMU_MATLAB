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
function [ProjPot] = MultiProjPot_conv_0(ScaleTypeCoord, CellNum, LattConst, Lx, Ly, Nx, Ny)
%MultiProjPot_conv_0.m calculates the projected potential of one slice
%using convolution if this slice contains more than one type of atoms.
%   ScaleTypeCoord -- matrix that contains the lattice information, where
%       the first row denotes the atomic types, the second to the forth row
%       denote the scaled atomic coordinates, since the involving atoms are 
%       all on one slice, the forth row may not be input;
%   CellNum -- expansion of the unit cell, syntax: [CellNumX, CellNumY];
%   LattConst -- lattice constants, syntax: [a, b];
%   Lx, Ly, Nx, Ny -- sampling parameters;
%   P.S.: I have come up with another way to further speed this process
%       when massive atoms are to be calculated, anyway, there are still
%       some bugs. When the bugs are fixed, the released version will be
%       named as MonoProjPot_conv_1.m and provide the same interface as
%       MonoProjPot_conv_0.m does, so it can be called in this function by
%       directly change the version number.

% sort ScaleTypeCoord in an ascending order by atomic type:
[SortedType, TypeOrder] = sort(ScaleTypeCoord(1, : ));
ScaleTypeCoord = ScaleTypeCoord( : , TypeOrder);
ProjPot = 0;
if length(SortedType) > 1
    LastType = 1;
    for Atom_Idx = 1 : length(SortedType)
        if SortedType(Atom_Idx)~=SortedType(LastType)
            ProjPot = ProjPot + MonoProjPot_conv_0(SortedType(LastType), ...
                ScaleTypeCoord(2 : 3 , LastType : Atom_Idx-1), CellNum, LattConst, Lx, Ly, Nx, Ny);
            LastType = Atom_Idx;
        end
    end
    ProjPot = ProjPot + MonoProjPot_conv_0(SortedType(LastType), ...
            ScaleTypeCoord(2 : 3 , LastType : Atom_Idx), CellNum, LattConst, Lx, Ly, Nx, Ny);
else
    ProjPot = MonoProjPot_conv_0(SortedType, ScaleTypeCoord(2 : 3, : ), CellNum, LattConst, Lx, Ly, Nx, Ny);
end

end

