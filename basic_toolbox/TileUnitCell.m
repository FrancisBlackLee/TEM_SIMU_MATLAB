function [superCell] = TileUnitCell(unitCell, tiles, tolerance, varargin)
%TileUnitCell.m tiles the input unit cell into a super cell.
%   unitCell -- input unit cell;
%   tiles -- tiling along three basis vectors, syntax: 
%       [tileA, tileB, tileC];
%   tolerance -- tolerance for distinguishing two different atoms;
%   varargin -- (coordType, lattConsts), syntax:
%       no input
%       ('frac') -- fractional coordinates:
%       (convMat) -- coordType = 'cart', and the conversion matrix;
%       ('cart', convMat) -- cartesian coordinates, and the conversion
%           matrix.

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

switch nargin
    case 2
        tolerance = 1.0e-8;
        baseUnitCell = RemoveSymmetricAtoms(unitCell, tolerance);
        bases = ConvMatToBases(eye(3));
    case 3
        baseUnitCell = RemoveSymmetricAtoms(unitCell, tolerance);
        bases = ConvMatToBases(eye(3));
    case 4
        if isnumeric(varargin{1})
            if size(varargin{1}, 1) == 3 && size(varargin{1}, 2) == 3
                baseUnitCell = RemoveSymmetricAtoms(unitCell, tolerance,...
                    varargin{:});
                convMat = varargin{1};
                bases = ConvMatToBases(convMat);
            else
                error('Invalid input of conversion matrix');
            end
        elseif ischar(varargin{1})
            if strcmp(varargin{1}, 'frac')
                baseUnitCell = RemoveSymmetricAtoms(unitCell, tolerance);
                bases = ConvMatToBases(eye(3));
            else
                error('Invalid input of coordinate type');
            end
        else
            error('Invalid input');
        end
    case 5
        if ischar(varargin{1}) && isnumeric(varargin{2})
            if strcmp(varargin{1}, 'cart') && size(varargin{2}, 1) == 3 &&...
                    size(varargin{2}, 2) == 3
                baseUnitCell = RemoveSymmetricAtoms(unitCell, tolerance,...
                    varargin{:});
                convMat = varargin{2};
                bases = ConvMatToBases(convMat);
            else
                error('Invalid input of coordinate type or conversion matrix');
            end
        else
            error('Invalid input of coordinate type or conversion matrix');
        end
    otherwise
        error('Invalid input');
end

na = tiles(1);
nb = tiles(2);
nc = tiles(3);

cellNum = na * nb * nc;
atomNumPerUnitCell = size(baseUnitCell, 2);
superCell = zeros(5, cellNum * atomNumPerUnitCell);

if any(mod(tiles, 1))
    error('Invalid tiles');
end

% tile
for tc = 0 : nc - 1
    for tb = 0 : nb - 1
        for ta = 0 : na - 1
            tmpCell = baseUnitCell;
            tmpCell(3 : 5, :) = tmpCell(3 : 5, :) + ta * bases.a' +...
                tb * bases.b' + tc * bases.c';
            head = (tc * na * nb + tb * na + ta) * atomNumPerUnitCell + 1;
            rear = head + atomNumPerUnitCell - 1;
            superCell(:, head : rear) = tmpCell;
        end
    end
end

end