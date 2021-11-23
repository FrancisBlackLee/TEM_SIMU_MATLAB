function [newTypeCoords] = RemoveSymmetricAtoms(oldTypeCoords, tolerance, varargin)
%RemoveSymmetricAtoms.m removes the extra symmetric atoms within a unit
%cell (2D/3D).
%   oldTypeCoord -- atomic type-coordinate matrix (default: fractional):
%       format: [T; P; X; Y; Z], note that T denotes atomic types 
%       represented by their atomic numbers, P denotes the atomic
%       concentration, X, Y, Z denote the atomic coordinates.
%   tolerance -- tolerance for distinguishing two different atoms;
%   varargin -- (coordType, lattConsts), syntax:
%       no input
%       ('frac') -- fractional coordinates:
%       (lattConsts) -- coordType = 'cart', and the lattice constants;
%       ('cart', lattConsts) -- cartesian coordinates, and the lattice
%           constants.

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

newTypeCoords = oldTypeCoords;
coordType = 'frac'; % default
switch nargin
    case 1
        tolerance = 1.0e-8;
    case 2
        % do nothing
    case 3
        if isnumeric(varargin{1})
            if length(varargin{1}) == 3
                coordType = 'cart';
                lattConsts = varargin{1};
                newTypeCoords(3, :) = newTypeCoords(3, :) / varargin{1}(1);
                newTypeCoords(4, :) = newTypeCoords(4, :) / varargin{1}(2);
                newTypeCoords(5, :) = newTypeCoords(5, :) / varargin{1}(3);
            else
                error('Invalid input of lattice constants');
            end
        elseif ischar(varargin{1})
            if ~strcmp(varargin{1}, 'frac')
                error('Invalid input of coordinate type');
            end
        else
            error('Invalid input');
        end
    case 4
        if ischar(varargin{1}) && isnumeric(varargin{2})
            if strcmp(varargin{1}, 'cart') && (length(varargin{2}) == 3)
                coordType = 'cart';
                lattConsts = varargin{2};
                newTypeCoords(3, :) = newTypeCoords(3, :) / varargin{2}(1);
                newTypeCoords(4, :) = newTypeCoords(4, :) / varargin{2}(2);
                newTypeCoords(5, :) = newTypeCoords(5, :) / varargin{2}(3);
            else
                error('Invalid input of coordinate type or lattice constants');
            end
        else
            error('Invalid input of coordinate type or lattice constants');
        end
    otherwise
        error('Invalid input');
end

% the possible symmetric atoms in a unit are such that have fractional 
% coordinates of more than or equal to one dimension who are equal to 0 or 1.

% remove the non-diagonal extra symmetric atoms:
RemoveNonDiagonalAtoms(3);
RemoveNonDiagonalAtoms(4);
RemoveNonDiagonalAtoms(5);

% remove the diagonal extra symmetric atoms:
% now symmetric atoms are on either the top or bottom surface along each
% dimension:
% move the atoms on the top surface to the bottom surface:
newTypeCoords(3, abs(newTypeCoords(3, :) - 1) < tolerance) = 0.0;
newTypeCoords(4, abs(newTypeCoords(4, :) - 1) < tolerance) = 0.0;
newTypeCoords(5, abs(newTypeCoords(5, :) - 1) < tolerance) = 0.0;

% there may be repeated atoms at the origin, remove the extra atoms
newTypeCoords = (uniquetol(newTypeCoords', tolerance, 'ByRows', true))';

% move the atoms on the bottom surface to the top surface along z dimension:
newTypeCoords(5, abs(newTypeCoords(5, :)) < tolerance) = 1;

if strcmp(coordType, 'cart')
    newTypeCoords(3, :) = newTypeCoords(3, :) * lattConsts(1);
    newTypeCoords(4, :) = newTypeCoords(4, :) * lattConsts(2);
    newTypeCoords(5, :) = newTypeCoords(5, :) * lattConsts(3);
end

% nested functions:
    function RemoveNonDiagonalAtoms(dim)
        borderAtoms = newTypeCoords(:, abs(newTypeCoords(dim, :)) < tolerance);
        borderAtomNum = size(borderAtoms, 2);
        [dim1, dim2] = OtherDims(dim);
        for borderAtomIdx = 1 : borderAtomNum
            removeIndices = find((abs(newTypeCoords(dim, :) - 1) < tolerance) &...
                (abs(newTypeCoords(dim1, :) - borderAtoms(dim1, borderAtomIdx)) < tolerance) &...
                (abs(newTypeCoords(dim2, :) - borderAtoms(dim2, borderAtomIdx)) < tolerance) &...
                (abs(newTypeCoords(1, :) - borderAtoms(1, borderAtomIdx)) < tolerance));
            newTypeCoords(:, removeIndices) = [];
        end
    end

    function [dim1, dim2] = OtherDims(dim)
        switch dim
            case 3
                dim1 = 4;
                dim2 = 5;
            case 4
                dim1 = 3;
                dim2 = 5;
            case 5
                dim1 = 3;
                dim2 = 4;
            otherwise
                error('Wrong input dimension!');
        end
    end

end

