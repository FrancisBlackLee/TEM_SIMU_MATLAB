function [name] = CrystalClassifier(cellLengths, cellAngles)
%CrystalClassifier.m classify the crystal systems based on the cell
%constants instead of the '_symmetry_cell_setting',
%'_space_group_crystal_system', '_symmetry_Int_Tables_number',
%'_space_group_IT_number', '_symmetry_space_group_name_H-M' and
%'_space_group_name_H-M_alt' that may or may not be contained or wrongly
%presented in CIF. Note that for trigonal crystal system, the primitive
%cell in hexagonal axes is classified as hexagonal, and the primitive cell
%in rhombohedral axes is classified as rhombohedral.
% Input:
%   cellLengths -- element 1 for cell length a, 2 for cell length b and 3
%       for cell length c;
%   cellAngles -- element 1 for cell angle alpha (between bases b and c);
%       2 for cell angle beta (between bases a and c) and
%       3 for cell angle gamma (between bases a and b);
% Output:
%   name -- name of the crystal system.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2023  Francis Black Lee (Li Xian)

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

a = cellLengths(1);
b = cellLengths(2);
c = cellLengths(3);

alpha = cellAngles(1);
beta = cellAngles(2);
gamma = cellAngles(3);

if (a == b) && (b == c) && (alpha == 90) && (beta == 90) && (gamma == 90)
    name = 'cubic';
elseif (a == b) && (b == c) && (alpha == beta) && (beta == gamma)
    name = 'rhombohedral';
elseif (a == b) && (alpha == beta) && (alpha == 90) && (gamma == 120)
    name = 'hexagonal';
elseif (a == b) && (alpha == 90) && (beta == 90) && (gamma == 90)
    name = 'tetragonal';
elseif (alpha == 90) && (beta == 90) && (gamma == 90)
    name = 'orthorhombic';
elseif ((alpha == 90) && (beta == 90)) || ((alpha == 90) && (gamma == 90))
    name = 'monoclinic';
else
    name = 'triclinic';
end

end


