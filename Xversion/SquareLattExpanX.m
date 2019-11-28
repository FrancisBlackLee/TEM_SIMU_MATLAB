%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019  Francis Black Lee and Li Xian

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outllok.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SliceOut] = SquareLattExpanX(SliceIn, LattConst, CellNum, DistError, MassNum, DebyeTemp)
%SquareLattExpan.m expands the sliced unit cell lattice by the input
%CellNum [Nx, Ny, Nz] and output the coordinates. Notice that the input SliceIn
%is the proportional coordinates of atoms in a unit cell or basic square
%lattice.
%   SliceIn -- fractional atomic coordinates (scaled to 0~1) and atomic 
%       types [T1, ..., TN; P1, ..., PN; x1,..., xN; y1,..., yN; z1,..., zN],
%       where T denotes atomic type, represented by the atomic numbers, P
%       denotes the elemental proportion and like before, x, y, z denote
%       the atomic coordinates;
%   LattConst -- planar lattice constants [a, b];
%   SliceDist -- distances between each two slices [D1,..., DN];
%   CellNum -- expansion numbers by which the SliceIn is expanded [Nx, Ny, Nz];
%   NOTE: X denotes an experimental version!
if nargin == 4
    % Transpose the matrix to improve speed
    SliceOut = SliceIn';
    CoordShift = CellNum / 2 + 1; % center the lattice
    % Expand along x
    SliceBase = SliceOut;
    for i = 1: CellNum(1) + 1
        SliceOut = [SliceOut; [SliceBase(:, 1), SliceBase(:, 2), SliceBase(:, 3) + i, SliceBase(:, 4), SliceBase(:, 5)]];
    end
    % Expand along y
    SliceBase = SliceOut;
    for i = 1: CellNum(2) + 1
        SliceOut = [SliceOut; [SliceBase(:, 1), SliceBase(:, 2), SliceBase(:, 3), SliceBase(:, 4) + i, SliceBase(:, 5)]];
    end
    % Expand along z
    SliceBase = SliceOut;
    for i = 1: CellNum(3) + 1
        SliceOut = [SliceOut; [SliceBase(:, 1), SliceBase(:, 2), SliceBase(:, 3), SliceBase(:, 4), SliceBase(:, 5) + i]];
    end
    SliceOut(:, 3) = SliceOut(:, 3) - CoordShift(1);
    SliceOut(:, 4) = SliceOut(:, 4) - CoordShift(2);
    SliceOut(:, 5) = SliceOut(:, 5) - CoordShift(3);
    [row, column] = find((SliceOut(:, 3) <= CellNum(1) / 2 + DistError) & (SliceOut(:, 3) >= -CellNum(1) / 2 - DistError) ...
                       & (SliceOut(:, 4) <= CellNum(2) / 2 + DistError) & (SliceOut(:, 4) >= -CellNum(2) / 2 - DistError) ...
                       & (SliceOut(:, 5) <= CellNum(3) / 2 + DistError) & (SliceOut(:, 5) >= -CellNum(3) / 2 - DistError));
    SliceOut = SliceOut(row, :);
    SliceOut = uniquetol(SliceOut, DistError, 'ByRows', true);
    SliceOut(:, 3) = LattConst(1) * SliceOut(:, 3);
    SliceOut(:, 4) = LattConst(2) * SliceOut(:, 4);
    SliceOut(:, 5) = LattConst(3) * SliceOut(:, 5);
    SliceOut = SliceOut';
else
    % Combine the SliceIn, MassNum and DebyeTemp:
    SliceIn = [SliceIn; MassNum; DebyeTemp];
    % Transpose the matrix to improve speed
    SliceOut = SliceIn';
    CoordShift = CellNum / 2 + 1; % center the lattice
    DistError = 1e-2;
    % Expand along x
    SliceBase = SliceOut;
    for i = 1: CellNum(1) + 1
        SliceOut = [SliceOut; [SliceBase(:, 1), SliceBase(:, 2), SliceBase(:, 3) + i, SliceBase(:, 4), SliceBase(:, 5), SliceBase(:, 6), SliceBase(:, 7)]];
    end
    % Expand along y
    SliceBase = SliceOut;
    for i = 1: CellNum(2) + 1
        SliceOut = [SliceOut; [SliceBase(:, 1), SliceBase(:, 2), SliceBase(:, 3), SliceBase(:, 4) + i, SliceBase(:, 5), SliceBase(:, 6), SliceBase(:, 7)]];
    end
    % Expand along z
    SliceBase = SliceOut;
    for i = 1: CellNum(3) + 1
        SliceOut = [SliceOut; [SliceBase(:, 1), SliceBase(:, 2), SliceBase(:, 3), SliceBase(:, 4), SliceBase(:, 5) + i, SliceBase(:, 6), SliceBase(:, 7)]];
    end
    SliceOut(:, 3) = SliceOut(:, 3) - CoordShift(1);
    SliceOut(:, 4) = SliceOut(:, 4) - CoordShift(2);
    SliceOut(:, 5) = SliceOut(:, 5) - CoordShift(3);
    [row, column] = find((SliceOut(:, 3) <= CellNum(1) / 2 + DistError) & (SliceOut(:, 3) >= -CellNum(1) / 2 - DistError) ...
                       & (SliceOut(:, 4) <= CellNum(2) / 2 + DistError) & (SliceOut(:, 4) >= -CellNum(2) / 2 - DistError) ...
                       & (SliceOut(:, 5) <= CellNum(3) / 2 + DistError) & (SliceOut(:, 5) >= -CellNum(3) / 2 - DistError));
    SliceOut = SliceOut(row, :);
    SliceOut = uniquetol(SliceOut, DistError, 'ByRows', true);
    SliceOut(:, 3) = LattConst(1) * SliceOut(:, 3);
    SliceOut(:, 4) = LattConst(2) * SliceOut(:, 4);
    SliceOut(:, 5) = LattConst(3) * SliceOut(:, 5);
    SliceOut = SliceOut';
end


end
